import h5py
import numpy as np
import pandas
from sklearn.manifold import TSNE
from tsnecuda import TSNE as cudaTSNE
from plotnine import *
from time import time

def condition(row) :
    if(row['methylated_frequency'] > 0.6) :
        return 1
    return 0

def mapProbeIds(data) :
    return data

def getExistingID(data, value) :
    for id in data :
        if id == value:
            return id

# ------------------- READ CpGs START -------------------

bed_file_path = "../megalodon/test_sample/modified_bases_fixed.5mC.bed"
cpg_columns = ["chromosome","start","end","name","score","strand","start_thick","end_thick","rgb","coverage","methylated_frequency"]
cpg = pandas.read_csv(bed_file_path, delimiter='\t', header=None, names=cpg_columns)

#Scale MAF
cpg['methylated_frequency'] = cpg['methylated_frequency'] / 100

probes_columns = ["chromosome","start","end","probeID","X","XX"]
probes450K = pandas.read_csv("../static/450K_hg19.bed", header=None, names=probes_columns, delimiter='\t')

overlap = pandas.merge(cpg, probes450K, on=["chromosome", "start"])

cpg_calls = overlap[(overlap['chromosome'] != "chrX") & (overlap["chromosome"] != "chrY")]
cpg_calls.drop_duplicates(subset=["probeID"], inplace=True)

# ------------------- READ CPGs DONE -------------------

cpg_calls['isMethylated'] = cpg_calls.apply(condition, axis=1)
cpg_calls = cpg_calls.set_index(pandas.Index(cpg_calls['probeID']))

case = cpg_calls.T
case = pandas.DataFrame(case.loc["isMethylated"])
case = case.T

filename = "/home/sander/nanopore_new/nanoDx_GPU/static/Trainingsets/Capper_et_al_old.h5"
file = h5py.File(filename, "r")

# Grab data from file
dx = pandas.DataFrame({'Dx': file['Dx']})

# Convert 'Dx' to a factor (categorical) variable
dx['Dx'] = pandas.Categorical(dx['Dx'])
print("DX shape : ", dx.shape)
print(dx)


sample_ids = file["sampleIDs"]
training_probes = np.array(file["probeIDs"]).astype(str)
print("Training probes shape : ", training_probes.shape)

tmp = np.array(case.columns)

# Get intersecting probes
probes = np.intersect1d(tmp, training_probes)
print("PROBES : ",probes.shape)
print(probes)

# # Get matching indicies
# idxs = np.where(np.isin(training_probes, probes))[0]
# print(len(probes)," overlapping CpG sites between sample and reference set. Reading training set now...")
# print("INDEX SHAPE : ",idxs)
#
# beta_values = file["betaValues"][:]
# print("Beta values shape : ",beta_values.shape)
#
#
# Create a DataFrame with Dx and the filtered matrix
# ts = pandas.DataFrame((beta_values > 0.6)[idxs].astype(int))
#
# ts = ts.T
# ts.columns = training_probes[idxs]
#
# constant_df = pandas.DataFrame({'Dx': ['unknown'] * len(probes)})
# m = pandas.concat([ts, constant_df], ignore_index=True)

# ---------------------------- PLOT THE TSNE ----------------------------

df = pandas.read_csv("../scripts/data.csv", index_col=0)
dx['Dx'] = dx['Dx'].astype(str)


print("Running TSNE")
start = time()
#tsne = TSNE(n_components=2, perplexity=30, n_iter=2500, verbose=1).fit_transform(df)
tsne = cudaTSNE(n_components=2, perplexity=30, n_iter=2500).fit_transform(df)
done = time()

print("Total elapsed time : ", done - start)
print(tsne)

print("Y : ",tsne[:, 1])
print("DX : ", dx)

df = pandas.DataFrame({"Dx":dx["Dx"].iloc[0], "X1": tsne[:,0], "X2": tsne[:,1]})

print("Final df : ", df)
print("columns: ", df.columns)


# Read data from a file into a DataFrame
colorMap = pandas.read_csv("../static/colorMap_Capper_et_al.txt", skip_blank_lines=True, delimiter="\t")

# Convert DataFrame to a Pandas DataFrame (equivalent to a tibble)
colorMap = pandas.DataFrame(colorMap)

# Group by 'group' column
grouped = colorMap.groupby('group')

# Sort within each group by 'methylationClass' column
colorMap = colorMap.groupby('group').apply(lambda x: x.sort_values('methylationClass'))

# Add a row with 'color' set to 'white' at the beginning of each group
colorMap = grouped.apply(lambda x: x.append({'color': 'white'}, ignore_index=True))

# Create the 'colorLabel' column based on 'methylationClass' and 'group'
colorMap['colorLabel'] = colorMap.apply(lambda row: f"**{row['group']}**" if pandas.isna(row['methylationClass']) else row['methylationClass'], axis=1)

# Replace missing values in 'color' column with 'grey'
colorMap['color'].fillna('grey', inplace=True)

# Set 'unknown' color to 'red'
colorMap.loc[colorMap['colorLabel'] == 'unknown', 'color'] = 'red'

# Extract the 'color' column to 'hexCol'
hexCol = colorMap['color']

# Rename the elements in 'hexCol' with 'colorLabel' values
hexCol.index = colorMap['colorLabel']

colorMap = colorMap.drop_duplicates(subset="colorLabel")
print("LABELS : ", colorMap['colorLabel'])

df['Dx'] = pandas.Categorical(df['Dx'], categories=colorMap['colorLabel'])


# Create the plot
p = (ggplot(df.query('Dx == "unknown"').sort_values('Dx'),
            aes(x='X1', y='X2', color='Dx', shape='Dx=="unknown", size=Dx=="unknown"')) +
     geom_point() +
     theme_classic() +
     ggtitle("t-SNE, perplexity = 30") +
     scale_color_manual(values=hexCol) +
     scale_shape_manual(values=[16, 3]) +
     scale_size_manual(values=[1, 4]) +
     guides(colour=guide_legend(title="Methylation class",
                                title_position="top",
                                override_aes={'shape': 15, 'size': 3}
                                )) +
     theme(legend_text=element_text(size=7))
)

ggsave(p, filename='python_plot.png', dpi=300)






file.close()

# Probes
# [1] "cg08750554" "cg26702212" "cg24517686" "cg14174502" "cg24709745"
#  [6] "cg08927102" "cg27302190" "cg18025136" "cg24638131" "cg14981312"
#
#  Indexes
#  [1] 125967 145349 143065 132257 143401 126160  10817 136563   9818   5888

# Final dataframe
#    cg13804450 cg14209957 cg14556482 cg13741263 cg13768503 cg04189705 cg14398988
# 1           1          0          1          0          0          0          1
# 2           1          0          0          0          1          0          1
# 3           1          1          0          0          1          0          1
# 4           1          1          1          0          1          0          1
# 5           1          1          1          0          1          0          1
# 6           1          1          1          0          0          0          1
# 7           1          1          1          0          1          0          1
# 8           1          1          0          0          0          0          1
# 9           1          0          1          0          1          0          1
# 10          1          1          1          0          0          0          1



