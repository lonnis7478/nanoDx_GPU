
rule hg19_GC_dist:
    input:
        ref=config["ref_genome"]
    output:
        chopped="intermediate/hg19_1K.fa",
        subsampled="intermediate/hg19_1K_subsample.fa"
    shell:
        "faSplit size {input} 1000 intermediate/hg19_1K -oneFile ; "
        "fasta-subsample {output.chopped} 20000 > {output.subsampled}"


rule make_GC_plot:
    input:
        BAM="bam/{sample}.bam",
        ref_chopped="static/hg19_1K_subsample.fa"
    output:
        GC_plot="figures/{sample}_GC.pdf",
	readlength_plot="figures/{sample}_readlength.pdf",
	length_dist="stats/{sample}-length_dist.RData"
    conda: "envs/basicR.yaml"
    script:
        "scripts/GC_histogram.R"


localrules: PDFreport_WGS, email_report
rule PDFreport_WGS:
    input:
        demux="stats/{sample}_demux_stats.txt",
        nanostat="stats/{sample}.nanostat.txt",
        CN="plots/{sample}-1000-0.05-CNplot.pdf",
        RFvotes="classification/{sample}-votes-RF5xCVrecal-{trainingSet}.RData",
        RFinfo="classification/{sample}-model_info-RF5xCVrecal-{trainingSet}.RData",
        DICT_FILE="static/{trainingSet}_dictionary.txt",
        GC="figures/{sample}_GC.pdf",
        RL="figures/{sample}_readlength.pdf",
        tSNE="plots/{sample}-tSNE-{trainingSet}.pdf",
        length_dist="stats/{sample}-length_dist.RData",
        mosdepth="stats/{sample}.mosdepth.summary.txt"
    output:
        "reports/{sample}_WGS_report_{trainingSet}.pdf"
    resources: pdfReport=1
    conda: "envs/PDFreport.yaml"
    script:
        "scripts/makePDFreport_WGS.R"


rule PDFreport_WGS_CUDA:
    input:
        demux="stats/{sample}_demux_stats.txt",
        nanostat="stats/{sample}.nanostat.txt",
        CN="plots/{sample}-1000-0.05-CNplot.pdf",
        DICT_FILE="static/{trainingSet}_dictionary.txt",
        GC="figures/{sample}_GC.pdf",
        RL="figures/{sample}_readlength.pdf",
        cuda_tSNE="plots/{sample}-tSNE-CUDA-{trainingSet}.pdf",
        length_dist="stats/{sample}-length_dist.RData",
        mosdepth="stats/{sample}.mosdepth.summary.txt",
        CUDAvotes="classification/{sample}-votes-CudaClassifier-{trainingSet}.RData",
        cuda_model_info="classification/{sample}-model_info-CudaClassifier-{trainingSet}.RData"
    output:
        "reports/{sample}_WGS_report.CUDA_{trainingSet}.pdf"
    resources: pdfReport=1
    conda: "envs/PDFreport.yaml"
    script:
        "scripts/makePDFreport_WGS.R"

rule PDFreport_WGS_CUDA_Fast:
    input:
        demux = "stats/{sample}_demux_stats.txt",
        nanostat = "stats/{sample}.nanostat.txt",
        CN = "plots/{sample}-1000-0.05-CNplot.pdf",
        GC="figures/{sample}_GC.pdf",
        RL="figures/{sample}_readlength.pdf",
        length_dist="stats/{sample}-length_dist.RData",
        mosdepth="stats/{sample}.mosdepth.summary.txt",
        RFvotes="classification/{sample}-votes-RF5xCVrecal-{trainingSet}.RData",
        RFinfo="classification/{sample}-model_info-RF5xCVrecal-{trainingSet}.RData",
        DICT_FILE="static/{trainingSet}_dictionary.txt",
        cuda_tSNE="plots/{sample}-tSNE-CUDA-{trainingSet}.pdf",
        CUDAvotes="classification/{sample}-votes-CudaClassifier-{trainingSet}.RData",
        cuda_model_info="classification/{sample}-model_info-RF5xCVrecal-{trainingSet}.RData"
    output:
        "reports/{sample}_WGS_report.CUDA.Fast_{trainingSet}.pdf"
    resources: pdfReport=1
    conda: "envs/PDFreport.yaml"
    script:
        "scripts/makePDFreport_WGS.R"
    


rule email_report:
    input: 
        report="reports/{sample}_WGS_report_{trainingSet}.pdf",
	    cnv="igv/{sample}-1000-0.05.seg"
    output: 
        touch("tmp/{sample}_WGS_report_{trainingSet}.sent")
    params:
        email=config["email"]
    shell: 
        "echo 'WGS report attached...' | mail -s 'WGS report {wildcards.sample}' -a {input.report} -a {input.cnv} {params.email}"

