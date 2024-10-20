import glob
import os

def list_fast5_pod5(wildcards):
  fast5 = glob.glob(config["FAST5_basedir"] + "/" + wildcards.sample + "/**/*.fast5", recursive=True)
  pod5 = glob.glob(config["FAST5_basedir"] + "/" + wildcards.sample +  "/**/*.pod5", recursive=True)
  print("Sample " + wildcards.sample + ": " + str(len(fast5)) + " FAST5, " + str(len(pod5)) + " POD5 files detected.")
  return fast5 + pod5

def list_basecall_targets(wildcards):
  files = list_fast5_pod5(wildcards)
  files.sort(key=os.path.getmtime)
  files = [i.replace(config["FAST5_basedir"], 'basecalling_parallel', 1) for i in files]
  return files[:config["max_fast5"]] if 'max_fast5' in config.keys() else files

rule setup_dorado:
  output: "../resources/tools/dorado-0.3.4-linux-x64/bin/dorado"
  shell: "mkdir -p resources/tools && wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.4-linux-x64.tar.gz -O - | tar -xz -C resources/tools"

rule setup_guppy_cpu:
  output: directory("../resources/tools/ont-guppy-cpu")
  shell: "mkdir -p ../resources/tools && wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.4.6_linux64.tar.gz -O - | tar -xz -C resources/tools"

rule setup_guppy_gpu:
  output: directory("../resources/tools/ont-guppy")
  shell: "mkdir -p ../resources/tools && wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_6.4.6_linux64.tar.gz -O - | tar -xz -C resources/tools"

rule download_model:
  input: "resources/tools/dorado-0.3.4-linux-x64/bin/dorado"
  output: directory("../resources/dorado_models/{model}")
  shell: "../resources/tools/dorado-0.3.4-linux-x64/bin/dorado download --model {wildcards.model} --directory resources/dorado_models"

if config["basecalling_mode"]=="dorado":
  localrules: basecall_dorado
  rule basecall_dorado:
    input:
      fast5 = config["FAST5_basedir"] + "/{sample}/{file}",
      model = "../resources/dorado_models/" + config["dorado_model"]
    output:
      out = temporary(directory("basecalling_parallel/{sample}/{file}")),
      tmp = temporary(directory("tmp/bc/{sample}/{file}"))
    threads: 4
    params: options = config["dorado_options"] if 'dorado_options' in config.keys() else "--modified-bases 5mCG_5hmCG" # by default, call 5mC + 5hmC because currently only combined models are available for R10.4.1 chemistry
    conda: "envs/minimap.yaml"
    resources: gpu=1
    shell:
      "mkdir {output.tmp} ; cp {input.fast5} {output.tmp} ; mkdir -p {output.out} ; "
      "../resources/tools/dorado-0.3.4-linux-x64/bin/dorado basecaller {input.model} {output.tmp}/ {params.options} | samtools view -bS - > {output.out}/{wildcards.file}.bam"

elif config["basecalling_mode"]=="guppy_cpu": 
  rule basecall_guppy:
    input:
      "../resources/tools/ont-guppy-cpu",
      fast5 = config["FAST5_basedir"] + "/{sample}/{file}"
    output:
      out = temporary(directory("../basecalling_parallel/{sample}/{file}")),
      tmp = temporary(directory("tmp/bc/{sample}/{file}"))
    params:
      model = config["guppy_model"],
      options = config["guppy_options"] if 'guppy_options' in config.keys() else ""
    threads: 4
    benchmark: "results/benchmarks/{sample}.{file}.basecall_guppy_cpu.benchmark.txt"
    shell:
      "mkdir {output.tmp} ; cp {input.fast5} {output.tmp} ; mkdir -p {output.out} ; "
      "../resources/tools/ont-guppy-cpu/bin/guppy_basecaller -i {output.tmp} -s {output.out} -c {params.model} --bam_out --disable_pings {params.options}"

elif config["basecalling_mode"]=="guppy_gpu":
  rule basecall_guppy:
    input:
      "../resources/tools/ont-guppy",
      fast5 = config["FAST5_basedir"] + "/{sample}/{file}"
    output:
      out = temporary(directory("../basecalling_parallel/{sample}/{file}")),
      tmp = temporary(directory("tmp/bc/{sample}/{file}"))
    params:
      model = config["guppy_model"],
      options = config["guppy_options"] if 'guppy_options' in config.keys() else "--device cuda:0"
    threads: 1
    resources: gpu=1
    benchmark: "results/benchmarks/{sample}.{file}.basecall_guppy_gpu.benchmark.txt"
    shell:
      "mkdir {output.tmp} ; cp {input.fast5} {output.tmp} ; mkdir -p {output.out} ; "
      "../resources/tools/ont-guppy/bin/guppy_basecaller -i {output.tmp} -s {output.out} -c {params.model} --bam_out --disable_pings {params.options}"


if config["basecalling_mode"]=="none_fastq":
  rule mirror_fastq:
    input: config["FAST5_basedir"] + "/{sample}"
    output: "results/fastq/{sample}.fq"
    conda: "results/envs/minimap.yaml"
    shell: "find {input} -name '*.fq' | xargs -n 1 cat > {output}"

elif config["basecalling_mode"]=="none_bam":
  rule mirror_modbam:
    input: config["FAST5_basedir"] + "/{sample}"
    output: "results/fastq/{sample}.fq"
    conda: "envs/minimap.yaml"
    benchmark: "results/benchmarks/{sample}/mirror_modbam.txt"
    shell: "find {input} -name '*.bam' | xargs -n 1 samtools fastq -T MM,ML > {output}"

else:
  rule merge_fastq:
    input: list_basecall_targets
    output: "results/fastq/{sample}.fq"
    conda: "envs/minimap.yaml"
    benchmark: "results/benchmarks/{sample}.merge_fastq.benchmark.txt"
    shell: "find results/basecalling_parallel/{wildcards.sample}/ -name '*.bam' | xargs -n 1 samtools fastq -T MM,ML > {output}"

rule demux_QC:
    input: "results/fastq/{sample}.fq"
    output:
        txt="results/stats/{sample}_demux_stats.txt",
        tmp=temp(directory("tmp/{sample}_qcat"))
    conda: "envs/demux.yaml"
    benchmark: "results/benchmarks/{sample}/demux_QC.txt"
    shell:
      "qcat -f {input} -b {output.tmp} 2> {output.txt}"
