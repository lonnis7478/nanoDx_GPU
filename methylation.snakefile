rule make_minimap_idx:
    input:
        ref=config["ref_genome"]
    output:
        idx="static/ref_genome.mmi"
    threads: 8
    conda: "envs/minimap.yaml"
    shell:
        "minimap2 -t {threads} -ax map-ont -d {output.idx} {input.ref}"


rule nanopolish_per_chunk:
    input:
        dir="basecalling_parallel/{dir}",
        ref=config["ref_genome"],
        idx="static/ref_genome.mmi"
    output:
        fq=temporary("tmp/{dir}/pass.fq"),
        bam=temporary("tmp/{dir}/pass.bam"),
        bai=temporary("tmp/{dir}/pass.bam.bai"),
        calls=temporary("methylation/{dir}/pass_calls.tsv")
    threads: 4
    conda: "envs/nanopolish.yaml"
    shell:
        "find basecalling_parallel/{wildcards.dir}/ -name '*.fastq' | xargs -i cat {{}} > {output.fq} || touch {output.fq} ; "
        "nanopolish index -vvv -d basecalling_parallel/{wildcards.dir} {output.fq} ; "
        "minimap2 -L -t 9 -ax map-ont {input.idx} {output.fq} | samtools view -bS - | samtools sort - > {output.bam} ; "
        "samtools index {output.bam} ; "
        "echo 'Calling nanopolish...' ; "
        "nanopolish call-methylation -vvv --progress --threads {threads} --reads {output.fq} --bam {output.bam} --genome {input.ref} --methylation cpg > {output.calls}"


def assemble_chunk_list(wildcards):
    dirs = subprocess.check_output("cd " + config["FAST5_basedir"] + "/" + wildcards.sample + "/ ; find . -name \"*.fast5\" -printf '%T+ %p\n' | sort | cut -f2 -d' '", shell=True).decode().split('\n')[:max_fast5]
    dirs = [i.replace('./', 'methylation/' + wildcards.sample + '/') for i in dirs]
    dirs = [i.replace('.fast5', '') for i in dirs]
    return [i + "/pass_calls.tsv" for i in dirs]


localrules: merge_calls
rule merge_calls:
    input:
         assemble_chunk_list
    output:
         "methylation/{sample}.calls"
    shell:
         "head -1 {input[0]} > {output} ; "
         "find methylation/{wildcards.sample} -name 'pass_calls.tsv' | xargs awk 'FNR>1' >> {output}"


rule pileup_CpG:
    input:
        "methylation/{sample}.calls"
    output:
        "methylation/{sample}.pileup.txt"
    shell:
        "python scripts/calculate_methylation_frequency.py -c 2.0 -s {input} > {output}"


localrules: nanopolish_to_bedMethyl
rule nanopolish_to_bedMethyl:
        input: "methylation/{sample}.pileup.txt"
        output: "bedMethyl/{sample}.CpG.bed"
        shell: "tail -n +2 {input} | awk '{{print $1 \"\\t\" $2 \"\\t\" ($3 + 1) \"\\t.\\t\" $5 \"\\t+\\t\" $2  \"\\t\" ($3 + 1) \"\\t\" int($7*255) \",0,\" int((1-$7)*255) \"\\t\" $5  \"\\t\" int($7*100)}}' > {output}"

rule megalodon_per_chunk:
    input:
        ref=config["ref_genome"],
    output:
        bed="megalodon/{sample}/modified_bases.5mC.bed"
    conda:"envs/megalodon.yaml"
    shell:
        "megalodon basecalling_parallel/{wildcards.sample} --guppy-server-path=/usr/bin/guppy_basecall_server --mappings-format bam"
        " --reference {input.ref} --devices 0 --processes 12  --overwrite"
        " --outputs  basecalls mappings mod_mappings mods --guppy-config=dna_r9.4.1_450bps_hac.cfg"
        " --remora-modified-bases dna_r9.4.1_e8 hac 0.0.0 5mc CG 0 --sort-mappings --mod-map-base-conv C T" 
        " --mod-map-base-conv m C --output-directory megalodon/{wildcards.sample}"

rule fix_bedMethyl:
    input:
        bed="megalodon/{sample}/modified_bases.5mC.bed"
    output:
        fixed="megalodon/{sample}/modified_bases_fixed.5mC.bed"
    conda:"envs/megalodon.yaml"
    shell:
        "python scripts/fix_bedMethyl.py {input.bed} {output.fixed}"

rule readCpGs:
    input:
        bed="megalodon/{sample}/modified_bases_fixed.5mC.bed",
        mapping="static/450K_hg19.bed"
    output:
        "methylation/{sample}.RData"
    conda: "envs/readCpGs.yaml"
    script:
        "scripts/readCpGs.R"

rule plot_tSNE:
    input:
        case="methylation/{sample}.RData",
        trainingset=config["trainingset_dir"] + "/{trainingset}.h5",
        colorMap="static/colorMap_{trainingset}.txt"
    output:
        pdf="plots/{sample}-tSNE-{trainingset}.pdf",
        html="plots/{sample}-tSNE-{trainingset}.html"
    benchmark: "benchmarks/{sample}.tSNE.{trainingset}.benchmark.txt"
    conda: "envs/tSNE.yaml"
    threads: 4
    script: "scripts/plot_tSNE.R"

rule transform_Rdata:
    input:
        meth="methylation/{sample}.RData"
    output:
        meth_transformed="transformed_rdata/{sample}-transformed.RData"
    conda: "envs/basicR.yaml"
    threads: 1
    script: "scripts/transform_Rdata.R"

rule Feature_importance_calculation:
    input:
        trainingset_meth=config["trainingset_dir"] + "/{trainingset}.h5"
    params:
        method = 'std'
    output:
        FeatureImportance="static/FeatureImportance-{trainingset}.p",
        FilteredFeature="static/FilteredFeature-{trainingset}.p"
    conda: "envs/pyClassifier.yaml"
    threads: 12
    script: "scripts/feature_importance.py"

rule Feature_selection_tfidf:
    input:
        meth="transformed_rdata/{sample}-transformed.RData",
        trainingset_meth=config["trainingset_dir"] + "/{trainingset}.h5",
        FeatureImportance = "static/FeatureImportance-{trainingset}.p",
        FilteredFeature = 'static/FilteredFeature-{trainingset}.p'
    params:
        max_CpG = config["max_CpG"]
    output:
        trainning_data="training/{sample}-FeatureSelection_idf-{trainingset}.p"
    conda: "envs/pyClassifier.yaml"
    threads: 12
    script: "scripts/feature_selection_tfidf.py"

rule RF5xCVrecal:
    input:
        data="training/{sample}-FeatureSelection_idf-{trainingset}.p",
        trainingset_meth=config["trainingset_dir"] + "/{trainingset}.h5",
        meth="transformed_rdata/{sample}-transformed.RData"
    output:
        pdf="plots/{sample}-RF5xCVrecal-{trainingset}.pdf",
        txt="classification/{sample}-votes-RF5xCVrecal-{trainingset}.txt",
        votes="classification/{sample}-votes-RF5xCVrecal-{trainingset}.RData",
        model_info="classification/{sample}-model_info-RF5xCVrecal-{trainingset}.RData"
    benchmark: "benchmarks/{sample}.RF5xCVrecal.{trainingset}.benchmark.txt"
    conda: "envs/pyClassifier.yaml"
    threads: 12
    script: "scripts/pyRF5xCVrecal.py"
