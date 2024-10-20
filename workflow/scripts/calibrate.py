import pandas as pd

base_dir = "/home/sander/Documents/sturgeon/result/"
file_path = base_dir+"NDMA1/merged_probes_methyl_calls_general.csv"


df = pd.read_csv(file_path)
print(df)
