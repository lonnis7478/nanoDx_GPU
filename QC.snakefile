
rule noReads:
    input:
        expand("fasta/{sample}_{quality}.fa",sample=SAMPLES,quality="pass fail all".split())
    output:
        "QC/noReads.txt"
    run:
        for f in input:
            shell("cat {f} | grep '>' | wc -l | echo \"{f} $(cat -)\" >> {output} || true")



rule fastQC:
    input:
        "fastq/{sample}.fq"
    output:
        directory("QC/{sample}/")
    threads: 12
    shell:
        "mkdir {output} ; fastqc -t {threads} {input} --outdir=QC/{wildcards.sample} --extract"