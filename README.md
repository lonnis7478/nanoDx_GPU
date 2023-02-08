# nanoDx_GPU pipeline

# Installation and configuration

First, we assume that you have installed the following dependencies:

 - conda (we use [miniconda](https://docs.conda.io/en/latest/miniconda.html))
 - snakemake (>=5.4.0) with DRMAA support
 - LaTeX (we use [TinyTeX](https://yihui.name/tinytex/), which can easily be installed from within R: `install.packages('tinytex') ; tinytex::install_tinytex()`)
 - guppy v3.1.5, ONT's proprietary basecaller (using the instructions and downloads in the ONT Community)

Then, you can checkout the pipeline:

`git clone https://gitlab.com/pesk/nanoDx.git`

Static data (training sets for methylation-based classification, pseudo germline BAM) can be downloaded here:

`https://file-public.bihealth.org/transient/nanodx/`

Some paths need to be set in the `config.yaml` file (see template `config.EXAMPLE.yaml`) in the pipeline directory to the FAST5 basedir (see Input data below) and the reference genome to be used:

```
ref_genome: /fast/projects/cubit/current/static_data/precomputed/BWA/0.7.15/hg19/ucsc/hg19.fa
ref_genome_masked: /fast/projects/cubit/current/static_data/precomputed/BWA/0.7.15/hg19/ucsc/hg19.fa
FAST5_basedir: /fast/projects/nanodx/runs/WGS
perform_basecalling: yes
```

The `perform_basecalling` flag can be set to `no` if your FAST5 files have been basecalled already (e.g. when sequenced on a GridION device).

You should also configure snakemake to work with your cluster. Example configuration files (`cluster_slurm.json` and `cluster_sge.json`) and launch scripts (`run_slurm.sh` and `run_sge.sh`) for use with SGE and SLURM are provided in the repo.

The pipeline is designed to parallelize base calling and methylation calling. 


# Input data


The pipeline takes one top level folder per sample as input. The base directory is set in the `config.yaml` file.

So for each sample, `<BASEDIR>/<samplename>/` should contain all FAST5 files. All subfolders are processed recursively. 
The current version of the pipeline requires multi-FAST5 files. If you have data with the old format, you can convert them using [single_to_multi_fast5](https://github.com/nanoporetech/ont_fast5_api) script provided by ONT.

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

# Current issues

- Rendering several reports in parallel (by specifying two or more target PDFs) fails due ambigious intermediate/temp file naming. You can work around this by specifing the (imaginary) pdfReport ressource when invoking snakemake: `--ressources pdfReport=1`. This option is automatically passed when using the `run_slurm.sh` or `run_sge.sh` wrappers, but you have to explicitly pass this when invoking snakemake directly. 
