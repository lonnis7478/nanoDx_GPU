#!/bin/bash

snakemake reports/NDMA99_SUB1_WGS_report.CUDA_Capper_et_al.pdf --use-conda --jobs 500 && snakemake reports/NDMA99_SUB2_WGS_report.CUDA_Capper_et_al.pdf --use-conda --jobs 500 && snakemake reports/NDMA29_SUB1_WGS_report.CUDA_Capper_et_al.pdf --use-conda --jobs 500 && snakemake reports/NDMA29_SUB2_WGS_report.CUDA_Capper_et_al.pdf --use-conda --jobs 500 && snakemake reports/NDMA32_SUB1_WGS_report.CUDA_Capper_et_al.pdf --use-conda --jobs 500 && snakemake reports/NDMA32_SUB2_WGS_report.CUDA_Capper_et_al.pdf --use-conda --jobs 500
