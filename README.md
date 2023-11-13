# nanoDx pipeline 

nanoDx is an open-source, end-to-end bioinformatics pipeline for DNA methylation-based classification of tumours using nanopore low-pass whole genome sequencing data. 

# News and updates

Please find the release history and changelog here: [https://gitlab.com/pesk/nanoDx/-/releases](https://gitlab.com/pesk/nanoDx/-/releases).

To stay informed about releases and news related to the nanoDx pipeline, you can self-subscribe to the nanodx-users mailing list (no regular postings): [https://mailman.charite.de/mailman/listinfo/nanodx-users](https://mailman.charite.de/mailman/listinfo/nanodx-users)

# Installation and configuration

First, we assume that you have installed the following dependencies:

 - conda (we use [miniconda](https://docs.conda.io/en/latest/miniconda.html))
 - snakemake (>=5.4.0)
 - LaTeX (we use [TinyTeX](https://yihui.name/tinytex/), which can easily be installed from within R: `install.packages('tinytex') ; tinytex::install_tinytex()`)

Then, you can install the pipeline using:

`git clone https://gitlab.com/pesk/nanoDx.git`

Static data (training sets for methylation-based classification, pseudo germline BAM) can be downloaded [here](https://charitede-my.sharepoint.com/:f:/g/personal/dongsheng_yuan_charite_de/EtXtX3enI0JEuYCQjc1_7iQB5EfEptXkdWGy7X4R40nOUw?e=Zry2bT).

Some paths need to be set in the `config.yaml` file (see template `config.EXAMPLE.yaml`) in the pipeline directory to the FAST5 basedir (see Input data below) and the reference genome to be used:

```
ref_genome: /path/to/hg19.fa
FAST5_basedir: /path/to/FAST5
basecalling_mode: guppy_cpu
guppy_model: 
```

## (Modified) base calling

Until version v0.5.1, modified bases (5mC) were called using [nanopolish](https://github.com/jts/nanopolish). However, nanopolish is no longer compatible to recent FAST5/POD5 file formats and sequencing chemistry.
We have therefore implemented base and 5mC calling using ONT's proprietary software (guppy or dorado) which needs to be configured using the `basecalling_mode` option:

 - `guppy_cpu`: basecalling will be performed using a recent guppy version (which will be automatically installed in the `resources/tools` subfolder) in a parallelized fashion (one job per FAST5 file) using CPU computation only.
This option is recommended for high-performance compute cluster with high CPU capacity but no or little GPU resources.
 - `guppy_gpu`: basecalling will be performed using a recent guppy version supporting GPU. Recommended for GPU-equipped workstations. Basecalling is performed for each FAST5 file individually.
 - `dorado`: basecalling using the experimental [dorado](https://github.com/nanoporetech/dorado) basecaller. Basecalling is performed for each FAST5 file individually. This allows incremental basecalling while sequencing without having to re-basecall the entire run. This option is recommended for GPU equipped workstations and near-realtime analysis.
 - `none_bam`: No base or modified base calling is performed. (Unaligned) modified BAM files output by MinKNOW or other pipelines are expected in the input folder (see `FAST5_basedir` config directive). This option requires that modified bases (5mC in CpG context) have been called with the correct model.
 - `none_fastq`: No base or modified base calling is performed. FASTQ containing methylation calls (via MM/ML tags) are expected in the input folder (see `FAST5_basedir` config directive). This option requires that modified bases (5mC in CpG context) have been called with the correct model.

The `perform_basecalling` flag is no longer supported.
Command line options (e.g. to configure CUDA devices) can be passed to basecallers using the `guppy_options` and `dorado_options` options, respectively, in the config file.


## Cluster vs. local work station execution

You should also configure snakemake to work with your cluster. Example configuration files (`cluster_slurm.json` and `cluster_sge.json`) and launch scripts (`run_slurm.sh` and `run_sge.sh`) for use with SGE and SLURM are provided in the repo.

The pipeline is designed to parallelize base calling and methylation calling. 


# Input data

The pipeline takes one top level folder per sample as input. The base directory is set in the `config.yaml` file using the `FAST5_basedir` option.

So for each sample, `<BASEDIR>/<samplename>/` should contain all FAST5 or POD5 files. All subfolders are processed recursively. 
The pipeline can handle multi-read FAST5 or POD5 files. If you have older single-read FAST5 data, you can convert them using [single_to_multi_fast5](https://github.com/nanoporetech/ont_fast5_api) script provided by ONT.

# Output

The pipeline currently performs two main analysis for low-pass WGS samples:

## CNV

A copy number profile is generated using the tumor sample and a pseudo-germline reference that is substracted (currently FAF04090.bam). Copy number profiles can be generated with different window size for binning of reads. To begin with, 1000kbp windows are recommended:

`snakemake plots/<SAMPLE>-CN-<WINDOWSIZE>-<PVAL_SEGMENTATION>-CNplot.pdf --use-conda`

Example:
`snakemake plots/<SAMPLE>-CN-1000-0.05-CNplot.pdf --use-conda`

## Methlyation-based classification

After calling 5mC methylation status for all rates, a random forest classifier is trained using all overlapping CpG sites between the sample and the reference set. Then, the tumor type of the sample is predicted using this classifier. Currently, we are using four training sets:

* `Euskirchen_et_al.h5` A small set of gliomas, MB and common brain metastases from the nanopore pilot study (Euskirchen, Acta Neuropathol 2017)

* `Sturm_et_al.h5` 16 different brain tumor entities with possible histological appearance as "PNET" (Sturm et al., Cell 2016)

* `Capper_et_al.h5` The full reference set of 80+ brain tumor entities from the Heidelberg reference cohort (Capper et al, Nature 2018)

* `TCGA.h5` All 30+ tumor entities from The Cancer Genome Atlas


`snakemake plots/<SAMPLE>-RF-<CLASSIFIER>.pdf --use-conda`

Example:
`snakemake plots/<SAMPLE>-RF-Capper_et_al.pdf --use-conda`

# Recommended use

Typically, we generate both a CN plot and methylation based classification using the Capper et al. training set to be combined in a sequencing report together with some quality metrics. We provide launch example launch scripts for SLURM and SGE. Using SLURM, a typical report can be generated using:

`./run_slurm.sh reports/<SAMPLE>_WGS_report_Capper_et_al.pdf`

This works for other training sets, too:

`./run_slurm.sh reports/<SAMPLE>_WGS_report_<TRAININGSET>.pdf`

#    Batch mode

Sometimes you might want to classify several samples at once. You can use generate a zipped archive of reports for a given training set using the batch_samples option in your main config file (or an additional one). The following command would use an additional config file to define your batch and generate PDF reports for all samples plus a ZIP archive of these reports:

`./run_slurm.sh --configfile my_batch.yaml reports/batch_reports_<TRAININGSET>.zip`

my_batch.yaml should then hold the sample IDs of your batch to be processed:
```
batch_samples:
- sample1
- sample2
```

# Current issues

- Rendering several reports in parallel (by specifying two or more target PDFs) fails due ambigious intermediate/temp file naming. You can work around this by specifing the (imaginary) pdfReport ressource when invoking snakemake: `--ressources pdfReport=1`. This option is automatically passed when using the `run_slurm.sh` or `run_sge.sh` wrappers, but you have to explicitly pass this when invoking snakemake directly. 
- ONT's dorado basecaller currently does not respect system-wide proxy settings (like https_proxy). If you are behind a proxy and choose dorado for basecalling, please set dorado-specific environment variables `dorado_proxy` and `dorado_proxy_port` to allow automatic downloading of models, e.g. 
```
export dorado_proxy='my.proxy.org'
export dorado_proxy_port=8080
```


# Citation

If you use the nanoDx pipeline in your research, *please consider citing the following papers* in your work:

[Kuschel LP, Hench J, Frank S, Hench IB, Girard E, Blanluet M, Masliah-Planchon J, Misch M, Onken J, Czabanka M, Yuan D, Lukassen S, Karau P, Ishaque N, Hain EG, Heppner F, Idbaih A, Behr N, Harms C, Capper D, Euskirchen P. Robust methylation-based classification of brain tumours using nanopore sequencing. Neuropathol Appl Neurobiol. 2023 Feb;49(1):e12856. doi: 10.1111/nan.12856.](https://doi.org/10.1111/nan.12856)

[Euskirchen P, Bielle F, Labreche K, Kloosterman WP, Rosenberg S, Daniau M, Schmitt C, Masliah-Planchon J, Bourdeaut F, Dehais C, Marie Y, Delattre JY, Idbaih A. Same-day genomic and epigenomic diagnosis of brain tumors using real-time nanopore sequencing. Acta Neuropathol. 2017 Nov;134(5):691-703. doi: 10.1007/s00401-017-1743-5.](https://doi.org/10.1007/s00401-017-1743-5)

# Disclaimer

Methylation-based classification using nanopore whole genome sequening is a research tool currently under development.
Interpretation and implementation of the results in a clinical setting is in the sole responsibility of the treating physician.
