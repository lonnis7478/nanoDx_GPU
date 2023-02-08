library(data.table)

# read bedMethyl
cpg <- fread(snakemake@input[["bed"]], data.table = F, header = F)
colnames(cpg) <- c("chromosome","start","end","name","score","strand","start_thick","end_thick","rgb","coverage","methylated_frequency")
cpg$methylated_frequency <- cpg$methylated_frequency / 100 # Scale MAF

probes450K <- fread(snakemake@input[["mapping"]], data.table = F)
colnames(probes450K) <- c("chromosome","start","end","probeID","X","XX")

overlap <- merge(cpg, probes450K, by=c("chromosome","start"))

# remove chrX + chrY probes
CpGcalls <- subset(overlap, !(chromosome %in% c("chrX","chrY")))

CpGcalls <- CpGcalls[!duplicated(CpGcalls$probeID),]

CpGcalls$isMethylated <- (CpGcalls$methylated_frequency > 0.6) * 1
rownames(CpGcalls) <- CpGcalls$probeID

case <- data.frame(t(CpGcalls))
case <- case["isMethylated",]

save(case, file = snakemake@output[[1]])
