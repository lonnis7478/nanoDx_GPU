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

rule plot_tSNE:
    input:
        case="results/methylation/{sample}.RData",
        trainingset=config["trainingset_dir"] + "/{trainingset}.h5",
        colorMap="static/colorMap_{trainingset}.txt"
    output:
        pdf="results/plots/{sample}-tSNE-{trainingset}.pdf",
        html="results/plots/{sample}-tSNE-{trainingset}.html"
    params:
        dim_reduction_method = config["dim_reduction_method"] if 'dim_reduction_method' in config.keys() else 'tsne',
        tsne_pca_dim = config["tsne_pca_dim"] if 'tsne_pca_dim' in config.keys() else 94,
        tsne_perplexity = config["tsne_perplexity"] if 'tsne_perplexity' in config.keys() else 30,
        tsne_max_iter = config["tsne_max_iter"] if 'tsne_max_iter' in config.keys() else 2500,
        umap_n_neighbours = config["umap_n_neighbours"] if 'umap_n_neighbours' in config.keys() else 10,
        umap_min_dist = config["umap_min_dist"] if 'umap_min_dist' in config.keys() else 0.5
    conda: "envs/tSNE.yaml"
    benchmark:"results/benchmarks/{sample}/plot_tSNE.{trainingset}.txt"
    threads: 4
    script: "scripts/plot_tSNE.R"


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


