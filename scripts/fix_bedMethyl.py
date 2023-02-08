import pandas as pd
import sys

file_path = sys.argv[1]
out_path  = sys.argv[2]
print("Fixing bedMethyl file : ", file_path)

columns = ["chromosome","start","end","name","score","strand","start_thick","end_thick","rgb","coverage","methylated_frequency"]

bedMethyl_file = pd.read_csv(file_path, sep='\t', header=None, names=columns)
bedMethyl_file.loc[bedMethyl_file.strand == '-', "start"] = bedMethyl_file['start'] -1
bedMethyl_file.loc[bedMethyl_file.strand == '-', "end"] = bedMethyl_file['end'] -1

bedMethyl_file.to_csv(out_path, header=False, index=False, sep="\t")

