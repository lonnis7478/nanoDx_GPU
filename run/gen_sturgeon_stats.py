
import numpy as np
import pandas as pd
import os

def sort(values):
    return values.str.extract('(\d+)$')[0].astype(int)


path = "/home/sander/Documents/sturgeon/result/"

dfs = []
score_threshold = 0.95

for folder in os.listdir(path):
    if os.path.exists(path+folder+"/merged_probes_methyl_calls_general.csv"):
        file = path+folder+"/merged_probes_methyl_calls_general.csv"
        df = pd.read_csv(file)
        df["sample"] = folder

        benchmark_file = open(path+folder+"/benchmark.txt")
        time = int(benchmark_file.read().split(":")[1])
        benchmark_file.close()
        df["time"] = time
        dfs.append(df)

combined_df = pd.concat(dfs, ignore_index=True)
combined_df.drop(columns=['number_probes'], inplace=True)
combined_df = combined_df.sort_values(by=["sample"],key=sort)
combined_df.set_index("sample", inplace=True)
print(combined_df)

max_columns = combined_df.drop(columns=["time"]).idxmax(axis=1)
print(max_columns)
final = pd.DataFrame(columns=["sample", "prediction", "score", "time"])
counter = 0
for column in max_columns:
    new_row = {"sample":max_columns.index[counter], "prediction":column, "score":combined_df.loc[max_columns.index[counter]][column], "time":combined_df.loc[max_columns.index[counter]]["time"]}
    final = pd.concat([final, pd.DataFrame([new_row])], ignore_index=True)
    counter += 1

final.set_index("sample", inplace=True)
final["time"] /= 1000
final["classifiable"] = final["score"] >= score_threshold
print(final)
final.to_csv("sturgeon.csv")

