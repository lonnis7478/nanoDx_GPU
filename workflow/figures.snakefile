
rule hg19_GC_dist:
    input:
        ref=config["ref_genome"]
    output:
        chopped="intermediate/hg19_1K-{sample}.fa",
        subsampled="intermediate/hg19_1K_subsample-{sample}.fa"
    benchmark: "results/benchmarks/{sample}/hg19_GC_dist.txt"
    shell:
        "faSplit size {input} 1000 intermediate/hg19_1K -oneFile ; "
        "fasta-subsample {output.chopped} 20000 > {output.subsampled}"


rule make_GC_plot:
    input:
        BAM="results/bam/{sample}.bam",
        ref_chopped="static/hg19_1K_subsample.fa"
    output:
        GC_plot="results/figures/{sample}_GC.pdf",
	readlength_plot="results/figures/{sample}_readlength.pdf",
	length_dist="results/stats/{sample}-length_dist.RData"
    conda: "envs/basicR.yaml"
    benchmark: "results/benchmarks/{sample}/make_GC_plot.txt"
    script:
        "scripts/GC_histogram.R"


localrules: PDFreport_WGS, email_report
rule PDFreport_WGS:
    input:
        demux="results/stats/{sample}_demux_stats.txt",
        nanostat="results/stats/{sample}.nanostat.txt",
        CN="results/plots/{sample}-1000-0.05-CNplot.pdf",
        RFvotes="results/classification/{sample}-votes-CudaClassifier-{trainingSet}.RData",
        RFinfo="results/classification/{sample}-model_info-CudaClassifier-{trainingSet}.RData",
        DICT_FILE="static/{trainingSet}_dictionary.txt",
        GC="results/figures/{sample}_GC.pdf",
        RL="results/figures/{sample}_readlength.pdf",
        tSNE="results/plots/{sample}-tSNE-CUDA-{trainingSet}.pdf",
        length_dist="results/stats/{sample}-length_dist.RData",
        mosdepth="results/stats/{sample}.mosdepth.summary.txt"
    output:
        "results/reports/{sample}_WGS_report_{trainingSet}.pdf"
    resources: pdfReport=1
    conda: "envs/PDFreport.yaml"
    benchmark: "results/benchmarks/{sample}/PDFreport_WGS.{trainingSet}.txt"
    script:
        "scripts/makePDFreport_WGS.R"


rule generate_statistics:
    input:
        report="results/reports/{sample}_WGS_report_Capper_et_al.pdf"
    output:
        stats="results/benchmarks/{sample}/stats.csv"
    conda: "envs/stats.yaml"
    script:
        "scripts/gen_stats.py"

rule email_report:
    input: 
        report="results/reports/{sample}_WGS_report_{trainingSet}.pdf",
	    cnv="results/igv/{sample}-1000-0.05.seg"
    output: 
        touch("tmp/{sample}_WGS_report_{trainingSet}.sent")
    params:
        email=config["email"]
    shell: 
        "echo 'WGS report attached...' | mail -s 'WGS report {wildcards.sample}' -a {input.report} -a {input.cnv} {params.email}"

