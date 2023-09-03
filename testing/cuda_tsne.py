import h5py
import numpy as np
import pandas

# ------------------- READ CpGs START -------------------

bed_file_path = "../megalodon/test_sample/modified_bases_fixed.5mC.bed"
cpg_columns = ["chromosome","start","end","name","score","strand","start_thick","end_thick","rgb","coverage","methylated_frequency"]
cpg = pandas.read_csv(bed_file_path, delimiter='\t', header=None, names=cpg_columns)

#Scale MAF
cpg['methylated_frequency'] = cpg['methylated_frequency'] / 100

probes_columns = ["chromosome","start","end","probeID","X","XX"]
probes450K = pandas.read_csv("../static/450K_hg19.bed", header=None, names=probes_columns, delimiter='\t')

overlap = pandas.merge(cpg, probes450K, on=["chromosome", "start"])

cpg_calls = overlap[(overlap['chromosome'] == "chrX") | (overlap["chromosome"] == "chrY")]
cpg_calls.drop_duplicates(subset=["probeID"], inplace=True)


# ------------------- READ CPGs DONE -------------------

def condition(row) :
    if(row['methylated_frequency'] > 0.6) :
        return 1
    return 0

def mapProbeIds(data) :
    return data


cpg_calls['isMethylated'] = cpg_calls.apply(condition, axis=1)

cpg_calls = cpg_calls.set_index(pandas.Index(cpg_calls['probeID']))


case = cpg_calls.T
case = pandas.DataFrame(case.loc["isMethylated"])

case = case.T
print(case.head())


filename = "/home/sander/nanopore_new/nanoDx_GPU/static/Trainingsets/Capper_et_al_old.h5"
file = h5py.File(filename, "r")


# Grab data from file
dx = file['Dx'];
sample_ids = file["sampleIDs"]
training_probes  = file["probeIDs"].asstr()
print(training_probes[29])

exit(0)


# Get intersecting probes
probes = case.columns.intersection(training_probes[:])
print(len(probes)," overlapping CpG sites between sample and reference set. Reading training set now...")
print(probes[9])

exit(0)
# Get matching indicies
idxs = probes.match(training_probes)


file.close()



