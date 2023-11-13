library(data.table)
METH_RDATA <- snakemake@input[["meth"]]
print(METH_RDATA)
load(METH_RDATA)
case=t(case)
save(case, file=snakemake@output[["meth_transformed"]])
message('Transformed')
