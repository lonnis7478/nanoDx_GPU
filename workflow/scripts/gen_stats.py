import os
import pandas as pd

sample = snakemake.wildcards.sample

df = pd.DataFrame(columns=['rule', 'time'])

for file in os.listdir("benchmarks/"+sample):
    rule = file.split(".")[0]
    bench_df = pd.read_csv("benchmarks/"+sample+"/"+file, sep="\t")

    df.loc[-1] = [rule, bench_df["s"].values[0]]  # adding a row
    df.index = df.index + 1  # shifting index
    df = df.sort_index()

df.loc[-1] = ["total", df['time'].sum()]  # adding a row
df.index = df.index + 1  # shifting index
df = df.sort_index()

df.to_csv(snakemake.output["stats"], index=False)

