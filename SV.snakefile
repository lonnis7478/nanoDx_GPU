
rule alignNGMLR:
    input:
        fq="fastq/{sample}.fq",
        ref=config["ref_genome"]
    output:
        "ngmlr/{sample}.bam"
    threads: 24
    conda: "envs/SV.yaml"
    shell:
        "ngmlr -t {threads} -r {input.ref} -q {input.fq} -x ont | samtools view -bS - | samtools sort - > {output}"

rule callSVs:
    input:
        bam="bam/{sample}.bam",
        bai="bam/{sample}.bam.bai",
    output:
        "sv/{sample}.{format}"
    wildcard_constraints:
        format="vcf|bedpe"
    conda: "envs/SV.yaml"
    shell:
        "sniffles -m {input.bam} --{wildcards.format} {output} --min_support 5"

rule circosPlot:
    input: "sv/{sample}.bedpe"
    output: "plots/{sample}.SV.circos.pdf"
    conda: "envs/SV.yaml"
    script: "scripts/plotSV.R"
