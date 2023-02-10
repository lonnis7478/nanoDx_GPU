import subprocess

rule demux_QC:
    input: "megalodon/{sample}/basecalls.fastq"
    output:
        txt="stats/{sample}_demux_stats.txt",
        tmp=temp(directory("tmp/{sample}_qcat"))
    conda: "envs/demux.yaml"
    shell: 
        "qcat -f {input} -b {output.tmp} 2> {output.txt}"
