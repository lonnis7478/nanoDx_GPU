#!/bin/bash


parent_dir="/home/sander/Documents/NDMA_subsampled"

for folder in "$parent_dir"/*/; do

  sample=$(basename "$folder")
  

  stats="results/benchmarks/${sample}/stats.csv"
  
  # Run the Snakemake command for the current sample
  snakemake --use-conda $stats --cores 12 --unlock
done


