
rule minimap2_align:
    input:
        fq="fastq/{sample}.fq",
        ref=config["ref_genome"]
    output:
        "bam/{sample}.bam"
    benchmark:
        "benchmarks/{sample}.minimap.align.benchmark.txt"
    threads: 12
    conda: "envs/minimap.yaml"
    shell:
        "minimap2 --MD -L -t 9 -ax map-ont -y {input.ref} {input.fq} | samtools view -bS - | samtools sort -T tmp/{wildcards.sample} - > {output} "

rule indexBAM:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    conda: "envs/minimap.yaml"
    benchmark: "benchmarks/{sample}.indexBAM.benchmark.txt"
    shell:
        "samtools index {input}"

rule coverage_wig:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai"
    output:
        "igv/{sample}.wig"
    threads: 1
    conda: "envs/igvtools.yaml"
    benchmark: "benchmarks/{sample}.coverage_wig.benchmark.txt"
    shell:
        "igvtools count -w 10000 {input.bam} {output} hg19"

rule coverage_mosdepth:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai"
    output:
        "stats/{sample}.mosdepth.summary.txt",
        "stats/{sample}.regions.bed.gz"
    threads: 1
    conda: "envs/mosdepth.yaml"
    benchmark: "benchmarks/{sample}.coverage_mosdepth.benchmark.txt"
    shell:
        "mosdepth -n --fast-mode --by 10000 stats/{wildcards.sample} {input.bam}"

rule QC_nanostat:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai"
    output: "stats/{sample}.nanostat.txt"
    threads: 12
    conda: "envs/qc.yaml"
    benchmark: "benchmarks/{sample}.QC_nanostat.benchmark.txt"
    shell: "NanoStat -t {threads} --bam {input.bam} --no_supplementary > {output}"

rule CN_profile:
    input:
        bam="bam/{sample}.bam",
        bai="bam/{sample}.bam.bai",
        normal_bam="static/FAF04090.bam"
    output:
        pdf="plots/{sample}-{binsize}-{alpha}-CNplot.pdf",
        rplot="plots/{sample}-CN-{binsize}-{alpha}.RData",
        bed="igv/{sample}-{binsize}-{alpha}.seg",
        bed1="igv/{sample}-{binsize}-{alpha}.bed"
    conda: "envs/createCNprofile.yaml"
    benchmark: "benchmarks/{sample}.{binsize}.{alpha}.CN_profile.benchmark.txt"
    script:
        "scripts/createCNprofile.R"

localrules: mirrorBAM
rule mirrorBAM:
    input:
        bam="bam/{sample}.bam",
        bai="bam/{sample}.bam.bai"
    output:
        dir=temp(directory("tmp/ACE/{sample}"))
    benchmark: "benchmarks/{sample}.mirrorBAM.benchmark.txt"
    shell: "mkdir {output} && cp {input} {output}"

rule ACE:
    input:
        bamDir="tmp/ACE/{sample}"
    output:
        dir=directory("ACE/{sample}")
    conda: "envs/ACE.yaml"
    benchmark: "benchmarks/{sample}.ACE.benchmark.txt"
    script: "scripts/CNV_ACE.R"
