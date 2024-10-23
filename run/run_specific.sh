#!/bin/bash


parent_dir="/home/sander/Documents/NDMA_subsampled_25K"
samples=('NDMA11' 'NDMA12' 'NDMA13' 'NDMA14' 'NDMA15' 'NDMA16' 'NDMA17' 'NDMA24'
 'NDMA27' 'NDMA28' 'NDMA29' 'NDMA30' 'NDMA31' 'NDMA32' 'NDMA38' 'NDMA44'
 'NDMA45' 'NDMA48' 'NDMA52' 'NDMA56' 'NDMA60' 'NDMA61' 'NDMA67' 'NDMA74'
 'NDMA76' 'NDMA82' 'NDMA84' 'NDMA85' 'NDMA87' 'NDMA88' 'NDMA90')



for sample in "${samples[@]}"; do
  echo "$sample"
  stats="results/benchmarks/${sample}/stats.csv"
  
  # Run the Snakemake command for the current sample
  snakemake --use-conda $stats --cores 12
done


