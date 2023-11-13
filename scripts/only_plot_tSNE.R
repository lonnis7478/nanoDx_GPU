library(ggtext)
library(ggplot2)
library(dplyr)
library(data.table)

df <- read.csv(snakemake@input[["tsne"]])

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
