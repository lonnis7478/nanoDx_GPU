library(Rtsne)
library(rhdf5)

### load methylation calls

load("/home/sander/nanopore_new/nanoDx_GPU/methylation/test_sample.RData")
case <- data.frame(t(sapply(case, function(x) as.numeric(as.character(x)))))

### load training set

fh5 = "/home/sander/nanopore_new/nanoDx_GPU/static/Trainingsets/Capper_et_al_old.h5"

# dump HDF5 training set content
h5ls(fh5)

Dx <- as.factor(h5read(fh5,"Dx"))
sampleIDs <- h5read(fh5,"sampleIDs")
trainingProbes <- h5read(fh5,"probeIDs")

probes <- intersect(colnames(case), trainingProbes)
idxs <- match(probes, trainingProbes)

print(head(probes, 10))
print(head(idxs, 10))

message(paste(length(probes)," overlapping CpG sites between sample and reference set. Reading training set now...",sep=""))

ts <- data.frame(Dx, (as.matrix(h5read(fh5, "betaValues")) > 0.6)[,idxs] * 1)
colnames(ts) <- c("Dx", trainingProbes[idxs])

print(colnames(ts))

m <- rbind(ts, data.frame(Dx = "unknown", case[,probes]))

dim(m)

# select most variable 50K probes

library(matrixStats)
beta <- as.matrix(m[,-1])
sds <- colSds(beta, na.rm=F)
maxSDs <- head(order(sds,decreasing=T),n=min(ncol(beta),50000))
print(maxSDs)

# set.seed(42)

write.csv(beta[,maxSDs], file = "data.csv")
tsne <- Rtsne(beta[,maxSDs], partial_pca = T, initial_dims = 94, perplexity = 30, theta = 0, max_iter = 2500, check_duplicates = F, verbose = T)

df <- data.frame(Dx = m[,1], tsne$Y)

#save(df, file="tsne_test.RData")

library(ggtext)
library(ggplot2)
library(dplyr)
library(data.table)

### create color scale

colorMap <- fread("../static/colorMap_Capper_et_al.txt", blank.lines.skip = TRUE) %>%
  as_tibble() %>%
  group_by(group) %>% 
  arrange(methylationClass) %>%
  group_modify(~ add_row(.x,.before=0,color="white")) %>%
  mutate(colorLabel = ifelse(is.na(methylationClass),paste0("**",group,"**"),methylationClass))

#colorMap <- merge(data.frame(methylationClass = unique(df$Dx)), colorMap, all.x = T)

hexCol <- colorMap$color
names(hexCol) <- colorMap$colorLabel

hexCol[is.na(hexCol)] <- "grey"
hexCol["unknown"] <- "red"

### reorder Dx factor levels

print(colorMap)
print(c(colorMap$colorLabel, "unknown"))
df$Dx <- factor(df$Dx, levels = c(colorMap$colorLabel, "unknown"))


### plot

p <- ggplot(data = df  %>% arrange(Dx=="unknown"), aes(x = X1, y = X2, color = Dx, shape = Dx=="unknown", size = Dx=="unknown")) + 
  geom_point() + 
  theme_classic() + 
  ggtitle("t-SNE, perplexity = 30") + 
  scale_colour_manual(values = hexCol, labels = names(hexCol), drop = F) +
  scale_shape_manual(values = c(16, 3)) +
  scale_size_manual(values=c(1,4)) + 
  guides(colour = guide_legend(title = "Methylation class", 
                               title.position = "top",
                               override.aes = list(shape = 15, size = 3)
                               )) +
  theme(legend.text = element_markdown(size = 7))

ggsave(plot = p, width = 14, height = 7, filename = snakemake@output[["pdf"]])

### interactive plot

library(plotly)
library(R.utils)

ip <- plot_ly(data = df  %>% arrange(Dx=="unknown"), 
              x = ~X1, y = ~X2, 
              text = ~Dx, 
              color= ~Dx, colors = hexCol,
              symbol=~Dx=="unknown", symbols=c("circle","+"),
              size=~ifelse(Dx=="unknown",1,0), sizes=c(10,500))

htmlwidgets::saveWidget(as_widget(ip), getAbsolutePath(snakemake@output[["html"]]), selfcontained = T, libdir = getAbsolutePath(paste(snakemake@output[["html"]],"_files",sep="")))

