configfile: "config.yaml"

include: "common.snakefile"
include: "basecall.snakefile"
include: "methylation.snakefile"
include: "figures.snakefile"
include: "demux.snakefile"
include: "SV.snakefile"

# target rule to batch classification

if 'batch_samples' in config.keys():
  rule batch_reports:
    input: expand("reports/{sample}_WGS_report_{{trainingset}}.pdf", sample=config['batch_samples'])
    output: "reports/batch_reports_{trainingset}.zip"
    shell: "zip {output} {input}"

