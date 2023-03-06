rule make_minimap_idx:
    input:
        ref=config["ref_genome"]
    output:
        idx="static/ref_genome.mmi"
    threads: 8
    conda: "envs/minimap.yaml"
    shell:
        "minimap2 -t {threads} -ax map-ont -d {output.idx} {input.ref}"

rule megalodon_per_chunk:
    input:
        ref=config["ref_genome"],
        basedir=config["FAST5_basedir"],
        guppy_path=config["guppy_basecall_server_path"]
        
    output:
        "megalodon/{sample}/modified_bases.5mC.bed",
        "megalodon/{sample}/mappings.sorted.bam",
        "megalodon/{sample}/mappings.bam",
        "megalodon/{sample}/mappings.sorted.bam.bai",
        "megalodon/{sample}/basecalls.fastq"
    conda:"envs/megalodon.yaml"
    params:
        guppy_config = config["guppy_config"],
        remora = config["remora_model"]
    shell:
        "megalodon {input.basedir}/{wildcards.sample} \
        --guppy-server-path {input.guppy_path} \
        --mappings-format bam \
        --reference {input.ref} \
        --devices 0 \
        --processes 16 \
        --overwrite \
        --outputs  basecalls mappings mod_mappings mods \
        --guppy-config {params.guppy_config} \
        --remora-modified-bases {params.remora} \
        --sort-mappings \
        --mod-map-base-conv m C \
        --output-directory megalodon/{wildcards.sample}"

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
