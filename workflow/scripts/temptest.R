library(Rtsne)
library(uwot)
library(rhdf5)
library(ggtext)
library(ggplot2)
library(dplyr)
library(data.table)

### create color scale

colorMap <- fread("static/colorMap_Capper_et_al.txt", blank.lines.skip = TRUE) %>%
	  as_tibble() %>%
	    group_by(group) %>%
	      arrange(methylationClass) %>%
	        group_modify(~ add_row(.x,.before=0,color="white")) %>%
		  mutate(colorLabel = ifelse(is.na(methylationClass),paste0("**",group,"**"),methylationClass))

hexCol <- colorMap$color
names(hexCol) <- colorMap$colorLabel

hexCol[is.na(hexCol)] <- "grey"
hexCol["unknown"] <- "red"

### load methylation calls

load("methylation/NDMA37.RData")
case <- data.frame(t(sapply(case, function(x) as.numeric(as.character(x)))))

### load training set

fh5 = "/home/sander/Documents/static/Trainingsets/Capper_et_al.h5"

# dump HDF5 training set content
h5ls(fh5)

Dx <- as.factor(h5read(fh5,"Dx"))

setdiff(Dx,colorMap$methylationClass) # check for missing elements in color map

sampleIDs <- h5read(fh5,"sampleIDs")
trainingProbes <- h5read(fh5,"probeIDs")

probes <- intersect(colnames(case), trainingProbes)
idxs <- match(probes, trainingProbes)

message(paste(length(probes)," overlapping CpG sites between sample and reference set. Reading training set now...",sep=""))

ts <- data.frame(Dx, (as.matrix(h5read(fh5, "betaValues")) > 0.6)[,idxs] * 1)
colnames(ts) <- c("Dx", trainingProbes[idxs])

m <- rbind(ts, data.frame(Dx = "unknown", case[,probes]))

dim(m)
write.csv(m, "NDMA37_m.csv")


### select most variable 50K probes

library(matrixStats)
beta <- as.matrix(m[,-1])
sds <- colSds(beta, na.rm=F)
maxSDs <- head(order(sds,decreasing=T),n=min(ncol(beta),50000))

write.csv(beta[,maxSDs], "NDMA37.csv")
# perform tSNE or UMAP reduction



