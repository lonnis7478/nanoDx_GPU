library(Rtsne)
library(uwot)
library(rhdf5)
library(ggtext)
library(ggplot2)
library(dplyr)
library(data.table)

### create color scale
colorMap <- fread(snakemake@input[["colorMap"]], blank.lines.skip = TRUE) %>%
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

load(snakemake@input[["case"]])
case <- data.frame(t(sapply(case, function(x) as.numeric(as.character(x)))))

### load training set

fh5 = snakemake@input[["trainingset"]]

# dump HDF5 training set content
h5ls(fh5)

Dx <- as.factor(h5read(fh5,"Dx"))

setdiff(Dx,colorMap$methylationClass) # check for missing elements in color map

sampleIDs <- h5read(fh5,"sampleIDs")
trainingProbes <- h5read(fh5,"probeIDs")

probes <- intersect(colnames(case), trainingProbes)
idxs <- match(probes, trainingProbes)


write(paste0(snakemake@wildcards[["sample"]],"\t",length(probes)), file=snakemake@params[["cpg_file"]], append=TRUE)
message(paste(length(probes)," overlapping CpG sites between sample and reference set. Reading training set now...",sep=""))

ts <- data.frame(Dx, (as.matrix(h5read(fh5, "betaValues")) > 0.6)[,idxs] * 1)
colnames(ts) <- c("Dx", trainingProbes[idxs])


m <- rbind(ts, data.frame(Dx = "unknown", case[,probes]))

dim(m)

### select most variable 50K probes

library(matrixStats)
beta <- as.matrix(m[,-1])
sds <- colSds(beta, na.rm=F)
maxSDs <- head(order(sds,decreasing=T),n=min(ncol(beta),snakemake@params[["max_CpG"]]))


if(snakemake@params[["save_dataframes"]] == "yes") {
    write.csv(beta[,maxSDs], snakemake@output[["beta"]])}

    write.csv(m[,1], snakemake@output[["m"]])
    quit()

# perform tSNE or UMAP reduction

if (snakemake@params[["dim_reduction_method"]] == "tsne") {
  tsne <- Rtsne(beta[,maxSDs],
                partial_pca = T,
                initial_dims = snakemake@params[["tsne_pca_dim"]],
                perplexity = snakemake@params[["tsne_perplexity"]],
                theta = 0,
                max_iter = snakemake@params[["tsne_max_iter"]], 
                check_duplicates = F, verbose = T)
  df <- data.frame(Dx = m[,1], tsne$Y)
} else if (snakemake@params[["dim_reduction_method"]] == "umap") {
  u <- umap(beta[,maxSDs],
            n_neighbors = snakemake@params[["umap_n_neighbours"]],
            min_dist = snakemake@params[["umap_min_dist"]],
            verbose = T)
  df <- data.frame(Dx = m[,1], u)
}

### reorder Dx factor levels

df$Dx <- factor(df$Dx, levels = c(colorMap$colorLabel, "unknown"))

### plot

p <- ggplot(data = df  %>% arrange(Dx=="unknown"), aes(x = X1, y = X2, color = Dx, shape = Dx=="unknown", size = Dx=="unknown")) + 
  geom_point() + 
  theme_classic() + 
  ggtitle(ifelse(snakemake@params[["dim_reduction_method"]]=="tsne",
                 paste0("t-SNE, no. of PCA dimensions = ", snakemake@params[["tsne_pca_dim"]], 
                        ", perplexity = ", snakemake@params[["tsne_perplexity"]],
                        ", max no. of iterations = ", snakemake@params[["tsne_max_iter"]]),
                 ifelse(snakemake@params[["dim_reduction_method"]]=="umap",
                        paste0("UMAP, no. of neighbours = ", snakemake@params[["umap_n_neighbours"]],
                        ", minimum distance = ", snakemake@params[["umap_min_dist"]]),
                        "")
                 )) +
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
