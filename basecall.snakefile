import subprocess

def listFAST5(wildcards):
    dirs = subprocess.check_output("cd " + config["FAST5_basedir"] + "/" + wildcards.sample + "/ ; find . -name \'*.fast5\' -type f -printf '%T+ %p\n' | sort | cut -f2 -d' '", shell=True).decode().split('\n')[:max_fast5]
    dirs = [i.replace('./', 'basecalling_parallel/' + wildcards.sample + '/') for i in dirs]
    return [i.replace('.fast5', '') for i in dirs]

if config["perform_basecalling"]:	
  rule basecall_nobarcode:
    input:
        config["FAST5_basedir"] + "/{sample}/{file}.fast5"
    output:
        out = temporary(directory("basecalling_parallel/{sample}/{file}")),
        tmp = temporary(directory("tmp/bc/{sample}/{file}"))
    threads: 4
    shell:
        "mkdir {output.tmp} ; cp {input} {output.tmp} ; "
        "guppy_basecaller --num_callers {threads} -c dna_r9.4.1_450bps_fast.cfg --input_path {output.tmp}/ --save_path {output.out}/ --fast5_out"
#        "guppy_basecaller --num_callers {threads} --flowcell FLO-MIN106 --kit SQK-RBK004 --input_path {output.tmp}/ --save_path {output.out}/ --fast5_out"
else:
  rule mirror_fast5:
    input:
      fast5 = config["FAST5_basedir"] + "/{sample}/{file}.fast5"
    wildcard_constraints: sample="[\w\-. ]+"
    output:
      directory("basecalling_parallel/{sample}/{file}")
    threads: 1
    conda: "envs/demux.yaml"
    shell: 
      "mkdir {output} && cp {input.fast5} {output} && python scripts/fast5_to_fastq.py {input.fast5} > {output}/reads.fastq"

rule merge_fastq:
    input:
        listFAST5
    output:
        "fastq/{sample}.fq"
    shell:
        "find basecalling_parallel/{wildcards.sample}/ -name '*.fastq' | xargs -i cat {{}} > {output}"

rule demux_QC:
    input: "fastq/{sample}.fq"
    output:
        txt="stats/{sample}_demux_stats.txt",
        tmp=temp(directory("tmp/{sample}_qcat"))
    conda: "envs/demux.yaml"
    shell: 
        "qcat -f {input} -b {output.tmp} 2> {output.txt}"
