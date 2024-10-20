library(ggtext)
library(ggplot2)
library(dplyr)
library(data.table)

df <- read.csv("temp_tsne_pca.csv")

colorMap <- fread("/home/sander/Documents/static/colorMap_Capper_et_al.txt", blank.lines.skip = TRUE) %>%
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
  scale_colour_manual(values = hexCol, labels = names(hexCol), drop = F) +
  scale_shape_manual(values = c(16, 3)) +
  scale_size_manual(values=c(1,4)) +
  guides(colour = guide_legend(title = "Methylation class",
                               title.position = "top",
                               override.aes = list(shape = 15, size = 3)
                               )) +
  theme(legend.text = element_markdown(size = 7))

ggsave(plot = p, width = 14, height = 7, filename = "tsne_temp_plot_pca.pdf")

### interactive plot

library(plotly)
library(R.utils)

ip <- plot_ly(data = df  %>% arrange(Dx=="unknown"),
              x = ~X1, y = ~X2,
              text = ~Dx,
              color= ~Dx, colors = hexCol,
              symbol=~Dx=="unknown", symbols=c("circle","+"),
              size=~ifelse(Dx=="unknown",1,0), sizes=c(10,500))

htmlwidgets::saveWidget(as_widget(ip), getAbsolutePath("tsne_temp_plot_pca.html"), selfcontained = T, libdir = getAbsolutePath(paste("tsne_temp_plot.html_pca","_files",sep="")))
