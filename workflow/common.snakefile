
rule minimap2_align:
    input:
        fq="results/fastq/{sample}.fq",
        ref=config["ref_genome"]
    output:
        "results/bam/{sample}.bam"
    benchmark:
        "results/benchmarks/{sample}/minimap2_align.txt"
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
    benchmark: "results/benchmarks/{sample}/indexBAM.txt"
    shell:
        "samtools index {input}"

rule coverage_wig:
    input:
        bam = "results/bam/{sample}.bam",
        bai = "results/bam/{sample}.bam.bai"
    output:
        "results/igv/{sample}.wig"
    threads: 1
    conda: "envs/igvtools.yaml"
    benchmark: "results/benchmarks/{sample}/coverage_wig.txt"
    shell:
        "igvtools count -w 10000 {input.bam} {output} hg19"

rule coverage_mosdepth:
    input:
        bam = "results/bam/{sample}.bam",
        bai = "results/bam/{sample}.bam.bai"
    output:
        "results/stats/{sample}.mosdepth.summary.txt",
        "results/stats/{sample}.regions.bed.gz"
    threads: 1
    conda: "envs/mosdepth.yaml"
    benchmark: "results/benchmarks/{sample}/coverage_mosdepth.txt"
    shell:
        "mosdepth -n --fast-mode --by 10000 results/stats/{wildcards.sample} {input.bam}"

rule QC_nanostat:
    input:
        bam = "results/bam/{sample}.bam",
        bai = "results/bam/{sample}.bam.bai"
    output: "results/stats/{sample}.nanostat.txt"
    threads: 12
    conda: "envs/qc.yaml"
    benchmark: "results/benchmarks/{sample}/QC_nanostat.txt"
    shell: "NanoStat -t {threads} --bam {input.bam} --no_supplementary > {output}"

rule CN_profile:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bam.bai",
        normal_bam="static/FAF04090.bam"
    output:
        pdf="results/plots/{sample}-{binsize}-{alpha}-CNplot.pdf",
        rplot="results/plots/{sample}-CN-{binsize}-{alpha}.RData",
        bed="results/igv/{sample}-{binsize}-{alpha}.seg",
        bed1="results/igv/{sample}-{binsize}-{alpha}.bed"
    conda: "envs/createCNprofile.yaml"
    benchmark: "results/benchmarks/{sample}/CN_profile.{binsize}-{alpha}.txt"
    script:
        "scripts/createCNprofile.R"

localrules: mirrorBAM
rule mirrorBAM:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bam.bai"
    output:
        dir=temp(directory("tmp/ACE/{sample}"))
    benchmark: "results/benchmarks/{sample}/mirrorBAM.txt"
    shell: "mkdir {output} && cp {input} {output}"

rule ACE:
    input:
        bamDir="tmp/ACE/{sample}"
    output:
        dir=directory("ACE/{sample}")
    conda: "envs/ACE.yaml"
    benchmark: "results/benchmarks/{sample}/ACE.txt"
    script: "scripts/CNV_ACE.R"
