
rule noReads:
    input:
        expand("fasta/{sample}_{quality}.fa",sample=SAMPLES,quality="pass fail all".split())
    output:
        "results/QC/noReads.txt"
    run:
        for f in input:
            shell("cat {f} | grep '>' | wc -l | echo \"{f} $(cat -)\" >> {output} || true")



rule fastQC:
    input:
        "results/fastq/{sample}.fq"
    output:
        directory("results/QC/{sample}/")
    threads: 12
    shell:
        "mkdir {output} ; fastqc -t {threads} {input} --outdir=results/QC/{wildcards.sample} --extract"