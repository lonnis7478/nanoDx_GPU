#!/bin/bash

PARENT_DIR=$1
TARGET_TYPE=$2

for folder in "$PARENT_DIR"/*; do
  sample=$(basename "$folder")
  if [ $TARGET_TYPE = "report" ]; then
    target="../results/reports/${sample}_WGS_report_Capper_et_al.pdf"
  elif [ $TARGET_TYPE = "stats" ]; then
    target="../results/benchmarks/${sample}/stats.csv"
  fi

  snakemake --use-conda "$target" --cores 12
done





