rule demux_qcat:
    input: "fastq/{run}.fq"
    output: directory("tmpDemux/{run}")
    conda: "envs/demux.yaml"
    benchmark: "benchmarks/{run}.demux_qcat.benchmark.txt"
    threads: 12
    shell: "qcat -f {input} -b tmpDemux/{wildcards.run} --min-score 40 -t {threads} --trim"

rule extract_fast5:
    input: 
      fast5 = config["FAST5_basedir"] + "/{run}/",
      fqDemux = "tmpDemux/{run}"
    output: 
      fast5dir = directory("demux/{run}_{barcode}"),
      ids = "tmpDemux/{run}_{barcode}.ids"
    conda: "envs/demux.yaml"
    shell: "awk '{{if(NR%4==1) print $1}}' {input.fqDemux}/{wildcards.barcode}.fastq | sed -e \"s/^@//\" > {output.ids} ; fast5_subset -i {input.fast5} -s {output.fast5dir} -l {output.ids} -r"
