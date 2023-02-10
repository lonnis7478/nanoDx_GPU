configfile: "config.yaml"

max_fast5 = config["max_fast5"] if 'max_fast5' in config.keys() else -1

include: "basecall.snakefile"
include: "methylation.snakefile"
include: "common.snakefile"
include: "figures.snakefile"
include: "demux.snakefile"
include: "SV.snakefile"

