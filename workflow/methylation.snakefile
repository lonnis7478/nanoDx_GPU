rule modbam_to_bedMethyl:
  input:
    bam = "results/bam/{sample}.bam",
    bai = "results/bam/{sample}.bam.bai",
    ref = config["ref_genome"]
  output: "results/bedMethyl/{sample}.CpG.bed"
  conda: "envs/minimap.yaml"
  benchmark:"results/benchmarks/{sample}/modbam_to_bedMethyl.txt"
  threads: 8
  shell: "modkit pileup {input.bam} {output} --ref {input.ref} --preset traditional --only-tabs --threads {threads}"

rule readCpGs:
    input:
        bed="results/bedMethyl/{sample}.CpG.bed",
        mapping="static/450K_hg19.bed"
    output:
        "results/methylation/{sample}.RData"
    conda: "envs/readCpGs.yaml"
    benchmark:"results/benchmarks/{sample}/readCpGs.txt"
    script:
        "scripts/readCpGs.R"


rule get_tSNE_data:
    input:
        case="results/methylation/{sample}.RData",
        trainingset=config["trainingset_dir"] + "/{trainingset}.h5",
        colorMap="static/colorMap_{trainingset}.txt",
    output:
        beta="results/betaValues/{sample}.{trainingset}.beta.csv",
        m="results/betaValues/{sample}.{trainingset}.m.csv"
    params:
        save_dataframes = "yes",
        cpg_file="results/benchmarks/CpGs_benchmark.txt",
        max_CpG=config["max_CpG"]
    conda: "envs/cudaTSNE.yaml"
    benchmark: "results/benchmarks/{sample}/get_tSNE_data.{trainingset}.txt"
    threads:4,
    script:"scripts/plot_tSNE.R"

rule tSNE_CUDA:
    input:
        beta="results/betaValues/{sample}.{trainingset}.beta.csv",
        m="results/betaValues/{sample}.{trainingset}.m.csv"
    output:
        tsne="results/cudaTSNE/{sample}.{trainingset}.tsne.csv"
    conda: "envs/cuml.yaml"
    benchmark: "results/benchmarks/{sample}/cuml_tsne.{trainingset}.txt"
    script: "scripts/cuda_tsne.py"

rule plot_tSNE_CUDA:
    input:
        tsne="results/cudaTSNE/{sample}.{trainingset}.tsne.csv",
        colorMap="static/colorMap_{trainingset}.txt"
    output:
        pdf = "results/plots/{sample}-tSNE-CUDA-{trainingset}.pdf",
        html = "results/plots/{sample}-tSNE-CUDA-{trainingset}.html"
    conda:"envs/tSNE.yaml"
    benchmark: "results/benchmarks/{sample}/plot_tSNE_CUDA.{trainingset}.txt"
    script:"scripts/only_plot_tSNE.R"

rule transform_Rdata:
    input:
        meth="results/methylation/{sample}.RData"
    output:
        meth_transformed="results/transformed_rdata/{sample}-transformed.RData"
    conda: "envs/basicR.yaml"
    benchmark: "results/benchmarks/{sample}/transform_Rdata.txt"
    threads: 1
    script: "scripts/transform_Rdata.R"

rule Feature_importance_calculation:
    input:
        trainingset_meth=config["trainingset_dir"] + "/{trainingset}.h5"
    params:
        method = 'std'
    output:
        FeatureImportance="static/FeatureImportance-{sample}-{trainingset}.p",
        FilteredFeature="static/FilteredFeature-{sample}-{trainingset}.p"
    conda: "envs/pyClassifier.yaml"
    benchmark: "results/benchmarks/{sample}/Feature_importance_calculation.{trainingset}.txt"
    threads: 12
    script: "scripts/feature_importance.py"

rule Feature_selection_tfidf:
    input:
        meth="results/transformed_rdata/{sample}-transformed.RData",
        trainingset_meth=config["trainingset_dir"] + "/{trainingset}.h5",
        FeatureImportance = "static/FeatureImportance-{sample}-{trainingset}.p",
        FilteredFeature = 'static/FilteredFeature-{sample}-{trainingset}.p'
    params:
        max_CpG = config["max_CpG"]
    output:
        trainning_data="results/training/{sample}-FeatureSelection_idf-{trainingset}.p"
    conda: "envs/pyClassifier.yaml"
    benchmark: "results/benchmarks/{sample}/Feature_selection_tfidf.{trainingset}.txt"
    threads: 12
    script: "scripts/feature_selection_tfidf.py"

rule CUDA_classifier:
    input:
        data="results/training/{sample}-FeatureSelection_idf-{trainingset}.p",
        trainingset_meth=config["trainingset_dir"] + "/{trainingset}.h5",
        meth="results/transformed_rdata/{sample}-transformed.RData",
    output:
        pdf="results/plots/cuda/{sample}-CudaClassifier-{trainingset}.pdf",
        txt="results/classification/{sample}-votes-CudaClassifier-{trainingset}.txt",
        votes="results/classification/{sample}-votes-CudaClassifier-{trainingset}.RData",
        model_info="results/classification/{sample}-model_info-CudaClassifier-{trainingset}.RData"
    params:
        max_CpG = config["max_CpG"]
    benchmark: "results/benchmarks/{sample}/CUDA_classifier.{trainingset}.txt"
    conda: "envs/pyClassifier.yaml"
    threads: 12
    script: "scripts/cuda_classifier.py"


rule RF5xCVrecal:
    input:
        data="results/training/{sample}-FeatureSelection_idf-{trainingset}.p",
        trainingset_meth=config["trainingset_dir"] + "/{trainingset}.h5",
        meth="results/transformed_rdata/{sample}-transformed.RData"
    output:
        pdf="results/plots/{sample}-RF5xCVrecal-{trainingset}.pdf",
        txt="results/classification/{sample}-votes-RF5xCVrecal-{trainingset}.txt",
        votes="results/classification/{sample}-votes-RF5xCVrecal-{trainingset}.RData",
        model_info="results/classification/{sample}-model_info-RF5xCVrecal-{trainingset}.RData"
    params:
        max_CpG = config["max_CpG"]
    benchmark: "results/benchmarks/{sample}/RF5xCVrecal.{trainingset}.txt"
    conda: "envs/pyClassifier.yaml"
    threads: 12
    script: "scripts/pyRF5xCVrecal.py"


