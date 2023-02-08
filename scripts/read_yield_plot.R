library(data.table)
library(ggplot2)

#stats <- fread("Z:/nanopore/revision.pipe/QC/noReads.txt")
stats <- fread(snakemake@input[["stats"]])
colnames(stats) <- c("file","reads")

stats$sample <- sub(".*/(\\w+).*","\\1", stats$file)
#stats$quality <- sub(".*/(\\w+)_(fail|pass|all).*","\\2", stats$file)

#runinfo <- fread("Z:/nanopore/pipeline/data/runs.txt")
runinfo <- fread(snakemake@input[["runinfo"]])
runinfo <- subset(runinfo, `library type`=="WGS")

stats <- merge(stats, runinfo, by.x="sample", by.y="run")

ggplot(stats, aes(x = sample, y = reads)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  # ggtitle("yield vs. read quality") + 
  scale_fill_brewer(type="qual", palette="Set1") +
  facet_grid(~ `kit version` + `flow cell version`, scales = "free_x", space = "free_x")

ggsave(snakemake@output[[1]], width=4,height=4)
