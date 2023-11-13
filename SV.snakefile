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
