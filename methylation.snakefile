rule modbam_to_bedMethyl:
  input:
    bam = "bam/{sample}.bam",
    bai = "bam/{sample}.bam.bai",
    ref = config["ref_genome"]
  output: "bedMethyl/{sample}.CpG.bed"
  conda: "envs/minimap.yaml"
  benchmark:"benchmarks/{sample}.modbam_to_bedMethyl.benchmark.txt"
  threads: 8
  shell: "modkit pileup {input.bam} {output} --ref {input.ref} --preset traditional --only-tabs --threads {threads}"

rule readCpGs:
    input:
        bed="bedMethyl/{sample}.CpG.bed",
        mapping="static/450K_hg19.bed"
    output:
        "methylation/{sample}.RData"
    conda: "envs/readCpGs.yaml"
    benchmark:"benchmarks/{sample}.readCpGs.benchmark.txt"
    script:
        "scripts/readCpGs.R"


rule get_tSNE_data:
    input:
        case="methylation/{sample}.RData",
        trainingset=config["trainingset_dir"] + "/{trainingset}.h5",
        colorMap="static/colorMap_{trainingset}.txt",
    output:
        beta="betaValues/{sample}.{trainingset}.beta.csv",
        m="betaValues/{sample}.{trainingset}.m.csv"
    params:
        save_dataframes = "yes",
        cpg_file="benchmarks/CpGs_benchmark.txt"
    conda: "envs/cudaTSNE.yaml"
    benchmark: "benchmarks/{sample}.{trainingset}.get_tSNE_data.benchmark.txt"
    threads:4,
    script:"scripts/plot_tSNE.R"

rule tSNE_CUDA:
    input:
        beta="betaValues/{sample}.{trainingset}.beta.csv",
        m="betaValues/{sample}.{trainingset}.m.csv"
    output:
        tsne="cudaTSNE/{sample}.{trainingset}.tsne.csv"
    conda: "envs/cuml.yaml"
    benchmark: "benchmarks/{sample}.{trainingset}.cuml_tsne.benchmark.txt"
    script: "scripts/cuda_tsne.py"

rule plot_tSNE_CUDA:
    input:
        tsne="cudaTSNE/{sample}.{trainingset}.tsne.csv",
        colorMap="static/colorMap_{trainingset}.txt"
    output:
        pdf = "plots/{sample}-tSNE-CUDA-{trainingset}.pdf",
        html = "plots/{sample}-tSNE-CUDA-{trainingset}.html"
    conda:"envs/tSNE.yaml"
    benchmark: "benchmarks/{sample}.{trainingset}.plot_tSNE_CUDA.benchmark.txt"
    script:"scripts/only_plot_tSNE.R"


rule plot_tSNE:
    input:
        case="methylation/{sample}.RData",
        trainingset=config["trainingset_dir"] + "/{trainingset}.h5",
        colorMap="static/colorMap_{trainingset}.txt"
    output:
        pdf="plots/{sample}-tSNE-{trainingset}.pdf",
        html="plots/{sample}-tSNE-{trainingset}.html"
    params:
        dim_reduction_method = config["dim_reduction_method"] if 'dim_reduction_method' in config.keys() else 'tsne',
        tsne_pca_dim = config["tsne_pca_dim"] if 'tsne_pca_dim' in config.keys() else 94,
        tsne_perplexity = config["tsne_perplexity"] if 'tsne_perplexity' in config.keys() else 30,
        tsne_max_iter = config["tsne_max_iter"] if 'tsne_max_iter' in config.keys() else 2500,
        umap_n_neighbours = config["umap_n_neighbours"] if 'umap_n_neighbours' in config.keys() else 10,
        umap_min_dist = config["umap_min_dist"] if 'umap_min_dist' in config.keys() else 0.5,
        save_dataframes = "no",
        cpg_file="benchmarks/CpGs_benchmark.txt",
    conda: "envs/tSNE.yaml"
    benchmark:"benchmarks/{sample}.{trainingset}.plot_tSNE.benchmark.txt"
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

rule CUDA_classifier:
    input:
        data="training/{sample}-FeatureSelection_idf-{trainingset}.p",
        trainingset_meth=config["trainingset_dir"] + "/{trainingset}.h5",
        meth="transformed_rdata/{sample}-transformed.RData",
    output:
        pdf="plots/cuda/{sample}-CudaClassifier-{trainingset}.pdf",
        txt="classification/{sample}-votes-CudaClassifier-{trainingset}.txt",
        votes="classification/{sample}-votes-CudaClassifier-{trainingset}.RData",
        model_info="classification/{sample}-model_info-CudaClassifier-{trainingset}.RData"
    benchmark: "benchmarks/{sample}.CUDA_classifier.{trainingset}.benchmark.txt"
    conda: "envs/pyClassifier.yaml"
    threads: 12
    script: "scripts/cuda_classifier.py"


rule benchmark:
    input:
        rep="reports/{sample}_WGS_report_CUDA_{trainingSet}.pdf",
    shell:
        "python scripts/benchmarking_pipeline.py {wildcards.sample}"




