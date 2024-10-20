import numpy as np
import pandas as pd
import os

def sort(values):
    return values.str.extract('(\d+)$')[0].astype(int)


path = "/home/sander/Documents/MethylzrFinal/MethyLYZR-main/output/"

dfs = []
score_threshold = 0.8
final = pd.DataFrame(columns=["sample", "cpgs", "prediction", "score", "time", "classifiable"])

for folder in os.listdir(path):
    if os.path.exists(f"{path}{folder}/{folder}_results.csv"):
        file = f"{path}{folder}/{folder}_results.csv"
        df = pd.read_csv(file, delimiter="\t", index_col=0)
        #df = df.loc[df["cpgs"] >= 7500].head(1)
        df = df.tail(1)

        benchmark_file = open(f"{path}{folder}/benchmark.txt")
        time = int(benchmark_file.read().split(":")[1])
        benchmark_file.close()

        new_row = {"sample": folder, "cpgs":df["cpgs"], "prediction":df["1st class"], "score":df["1st pred"], "time":time, "classifiable":df["1st pred"] > score_threshold}


        dfs.append(pd.DataFrame(new_row))



combined_df = pd.concat(dfs, ignore_index=True)
combined_df = combined_df.sort_values(by=["sample"],key=sort)
combined_df.set_index("sample", inplace=True)
combined_df["time"] /= 1000
print(np.mean(combined_df["cpgs"].values))
combined_df.to_csv("methylzr.csv")