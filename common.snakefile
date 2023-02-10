

rule filterMAPQ:
    input:
        "megalodon/{sample}.bam"
    output:
        "bam_mapq/{sample}.bam"
    params:
        minMAPQ="20"
    conda: "envs/minimap.yaml"
    shell:
        "samtools view -b -q {params.minMAPQ} {input} > {output}"

rule coverage_wig:
    input:
        bam = "megalodon/{sample}/mappings.sorted.bam",
        bai = "megalodon/{sample}/mappings.sorted.bam.bai"
    output:
        "igv/{sample}.wig"
    threads: 1
    conda: "envs/igvtools.yaml"
    shell:
        "igvtools count -w 10000 {input.bam} {output} hg19"

rule coverage_mosdepth:
    input:
        bam = "megalodon/{sample}/mappings.sorted.bam",
        bai = "megalodon/{sample}/mappings.sorted.bam.bai"
    output:
        "stats/{sample}.mosdepth.summary.txt",
        "stats/{sample}.regions.bed.gz"
    threads: 1
    conda: "envs/mosdepth.yaml"
    shell:
        "mosdepth -n --fast-mode --by 10000 stats/{wildcards.sample} {input.bam}"

rule QC_nanostat:
    input:
        bam = "megalodon/{sample}/mappings.sorted.bam",
        bai = "megalodon/{sample}/mappings.sorted.bam.bai"
    output:
        txt = "stats/{sample}.nanostat.txt"
	
    conda: "/envs/qc.yaml"
    shell:
        "NanoStat -t {threads} --bam {input.bam} > {output.txt}"
        

rule CN_profile:
    input:
        bam="megalodon/{sample}/mappings.sorted.bam",
        bai="megalodon/{sample}/mappings.sorted.bam.bai",
        normal_bam=config["bam_ref"]
    output:
        pdf="plots/{sample}-{binsize}-{alpha}-CNplot.pdf",
        rplot="plots/{sample}-CN-{binsize}-{alpha}.RData",
        bed="igv/{sample}-{binsize}-{alpha}.seg",
        bed1="igv/{sample}-{binsize}-{alpha}.bed"
    conda: "envs/createCNprofile.yaml"
    script:
        "scripts/createCNprofile.R"

localrules: mirrorBAM
rule mirrorBAM:
    input:
        bam="megalodon/{sample}/mappings.sorted.bam",
        bai="megalodon/{sample}/mappings.sorted.bam.bai"
    output:
        dir=temp(directory("tmp/ACE/{sample}"))
    shell: "mkdir {output} && cp {input} {output}"

rule ACE:
    input:
        bamDir="tmp/ACE/{sample}"
    output:
        dir=directory("ACE/{sample}")
    conda: "envs/ACE.yaml"
    script: "scripts/CNV_ACE.R"
