library(rmarkdown)
library(dplyr)

CN <- normalizePath(snakemake@input[["CN"]])
DICT_FILE <- normalizePath(snakemake@input[["DICT_FILE"]])
GC <- normalizePath(snakemake@input[["GC"]])
RL <- normalizePath(snakemake@input[["RL"]])
tSNE <- normalizePath(snakemake@input[["tSNE"]])
DEMUX <- normalizePath(snakemake@input[["demux"]])
NANOSTAT <- normalizePath(snakemake@input[["nanostat"]])
meanDepth <- system(paste("awk 'END{print $4}' ", normalizePath(snakemake@input[["mosdepth"]]), sep=""), intern=T)

load(normalizePath(snakemake@input[["length_dist"]]))
load(normalizePath(snakemake@input[["RFvotes"]]))
row.names(votes) <- votes$Dx
load(normalizePath(snakemake@input[["RFinfo"]]))

if (snakemake@config[["tsne_method"]] == "gpu") {
    cudaTSNE <- normalizePath(snakemake@input[["cuda_tSNE"]])

    render("scripts/WGSreportCUDA.Rmd",
           output_format = "pdf_document",
           output_dir=dirname(normalizePath(snakemake@output[[1]])),
           output_file=basename(normalizePath(snakemake@output[[1]]))
           )
}else{
    render("scripts/WGSreport.Rmd",
           output_format = "pdf_document",
           output_dir=dirname(normalizePath(snakemake@output[[1]])),
           output_file=basename(normalizePath(snakemake@output[[1]]))
           )
}
