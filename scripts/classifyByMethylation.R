library(data.table)
library(ranger)
library(RColorBrewer)
library(glmnet)

### load methylation calls

METH_RDATA <- snakemake@input[["meth"]]

load(METH_RDATA)
case <- data.frame(t(sapply(case, function(x) as.numeric(as.character(x)))))

### generate methylation training set

library(rhdf5)
fh5 = snakemake@input[["trainingset_meth"]]

# dump HDF5 training set content
h5ls(fh5)

Dx <- as.factor(h5read(fh5,"Dx"))
sampleIDs <- h5read(fh5,"sampleIDs")
trainingProbes <- h5read(fh5,"probeIDs")

probes <- intersect(colnames(case), trainingProbes)
idxs <- match(probes, trainingProbes)

message(paste(length(probes)," overlapping CpG sites between sample and reference set. Reading training set now...",sep=""))

ts <- data.frame(Dx, (as.matrix(h5read(fh5, "betaValues")) > 0.6)[,idxs] * 1)
colnames(ts) <- c("Dx", trainingProbes[idxs])

#### remove cases without Dx
ts <- subset(ts, Dx != "NA")
ts$Dx <- droplevels(ts$Dx)

message(paste(length(probes)," overlapping CpG sites read from training set.",sep=""))

#### feature selection (top half by SD)

max_CpG <- snakemake@params[["max_CpG"]]
library(matrixStats)
sds <- colSds(as.matrix(ts[,-1]), na.rm=F)
maxSDs <- head(order(sds,decreasing=T),n = min(ncol(ts)-1, max_CpG))
ts <- ts[,c(1,maxSDs + 1)]

message(paste(length(maxSDs)," features selected for training random forest.",sep=""))

### classify

cols <- intersect(colnames(ts),colnames(case))

ts$Dx <- as.factor(ts$Dx)

Dx_fractions <- min(summary(ts$Dx,maxsum=10000)) / summary(ts$Dx,maxsum=10000) # calculate the fraction of samples of each Dx corresponding to number of cases in smallest Dx class

rf <- ranger(dependent.variable.name = "Dx", data = ts[,c("Dx",cols)], num.trees=20000, probability = T, sample.fraction = Dx_fractions)

message(print(rf))

### Recalibration by training GLMs on the output scores

# Get score for all train samples
probs <- predict(rf, ts, predict.all = F)$predictions
scores <- unlist(lapply(1:dim(probs)[1], function(i){probs[i,which(probs[i,] == max(probs[i,]))]}))
pred <- attr(scores, "name")
classes <- levels(ts$Dx)

# A GLM is trained for each class individually (1 vs all)
glm_models <- lapply(classes, function(type){
    scores <- probs[, type]
    scale_this <- data.frame(class = ifelse(pred == type & pred == ts$Dx, 1, 0), score = scores)
    scaled_scores <- glm(class ~ score, scale_this, family = binomial)
    return(scaled_scores)
})

### predict case

#votes <- predict(rf, case[,cols], type="vote")

# x is now an array of numbers, the predicted class is now called x_pred
x_probs <- predict(rf, rbind(case[,cols],case[,cols]), predict.all=F)$predictions[1,]
x_score <- x_probs[which(x_probs == max(x_probs))]
x_pred <- attr(x_score, "name")


votes <- data.frame(x_probs)
colnames(votes) <- c("Freq")

votes$Freq <- votes$Freq / sum(votes$Freq) * 100

# Apply recalibration to case score
x_scaled <- lapply(levels(ts$Dx), function(type){
    x_scores <- x_probs[type]
    scaled_scores <- glm_models[[which(classes == type)]]
    x_scaled <- predict(scaled_scores, newdata = data.frame(score = x_scores), type = "response")
    return(x_scaled)
})
x_mean_this <- unlist(x_scaled)
x_calibrated_scores <- x_mean_this/sum(x_mean_this)
x_calibrated_score <- x_calibrated_scores[x_pred]

votes$cal_Freq <- x_calibrated_scores
votes$cal_Freq <- votes$cal_Freq / sum(votes$cal_Freq) * 100

votes <- votes[order(votes$Freq),, drop = FALSE] 

### Save calibration report

report <- paste(paste0("Number of features: ", rf$num.independent.variables), 
                paste0("Predicted Class: ", x_pred),
                paste0("Initial Score: ", x_score),
                paste0("Calibrated Score: ", x_calibrated_score),
                sep = "\n")
write.table(report, file = snakemake@output[["calibration"]], row.names = F, col.names = F, quote = F)

### save data for tSNE plot

m <- rbind(ts[,c("Dx",cols)],
	   data.frame(Dx="unknown",case[,cols]))

save(m, file=snakemake@output[["matrix"]])

# remove forest from rf object to save disk space
rf$forest <- NULL

save(rf, votes, x_pred, file=snakemake@output[["votes"]])

### plot to PDF

pdf(snakemake@output[["pdf"]])

pie(votes$Freq, labels = paste(rownames(votes)," (",as.integer(votes$Freq)," %)",sep=""), 
    init.angle=90, 
    clockwise=F, 
    main = snakemake@wildcards[["sample"]],
    col=brewer.pal(ncol(votes),"Set1"))

pie(votes$cal_Freq, labels = paste(rownames(votes)," (",round(votes$cal_Freq, digits = 1)," %)",sep=""), 
    init.angle=90, 
    clockwise=F, 
    main = snakemake@wildcards[["sample"]],
    col=brewer.pal(ncol(votes),"Set1"))

dev.off()

