import numpy as np
import pandas as pd
import os

def sort(values):
    return values.str.extract('(\d+)$')[0].astype(int)


meth_time = pd.DataFrame(columns=['Sample', 'Classification'])
path = "/home/sander/Documents/MethylzrFinal/MethyLYZR-main/output/"

counter = 0
for folder in os.listdir(path):
    if os.path.exists(f"{path}{folder}/benchmark.txt"):
        file = open(f"{path}{folder}/benchmark.txt")
        time = int(file.read().split(":")[1])/1000
        meth_time.loc[counter] = [folder, time]
        counter += 1

meth_time = meth_time.sort_values(by=["Sample"],key=sort)
print(meth_time)
meth_time.to_csv("methylzr_time.csv", index=False)
exit(0)



paths = ["/home/sander/Documents/nanoDXNN/benchmarks_2/"]
time_df = None

for folder in os.listdir(paths[0]):
    if os.path.exists(f"{paths[0]}{folder}/stats.csv"):
        df = pd.read_csv(f"{paths[0]}{folder}/stats.csv").set_index('rule').T
        df = df.reset_index(drop=True)
        df.index.name = 'Index'
        df["Sample"] = folder
        if time_df is None:
            print(df)
            time_df = pd.DataFrame(columns=df.columns)

        time_df = pd.concat([time_df, df])

time_df = time_df.sort_values(by=["Sample"],key=sort)

print(time_df)
time_df.to_csv("nn_time.csv")
