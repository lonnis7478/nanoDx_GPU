import h5py
import numpy as np
import pandas
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
#from tsnecuda import TSNE as cudaTSNE
from cuml.manifold import TSNE as cumlTSNE

from time import time
from datetime import datetime



def condition(row) :
    if(row['methylated_frequency'] > 0.6) :
        return 1
    return 0

def read_cpgs():
    bed_file_path = "/home/sander/Documents/nanoDx-master/bedMethyl/NDMA37.CpG.bed"
    cpg_columns = ["chromosome","start","end","name","score","strand","start_thick","end_thick","rgb","coverage","methylated_frequency"]
    cpg = pandas.read_csv(bed_file_path, delimiter='\t', header=None, names=cpg_columns)

    #Scale MAF
    cpg['methylated_frequency'] = cpg['methylated_frequency'] / 100

    probes_columns = ["chromosome","start","end","probeID","X","XX"]
    probes450K = pandas.read_csv("../static/450K_hg19.bed", header=None, names=probes_columns, delimiter='\t')

    print(cpg)
    print("\n------\n")
    print(probes450K)
    overlap = pandas.merge(cpg, probes450K, on=["chromosome", "start"])

    cpg_calls = overlap[(overlap['chromosome'] != "chrX") & (overlap["chromosome"] != "chrY")]
    cpg_calls.drop_duplicates(subset=["probeID"], inplace=True)

    return cpg_calls

def build_dataframe(case):
    filename = "/home/sander/Documents/static/Trainingsets/Capper_et_al.h5"
    file = h5py.File(filename, "r")

    # Grab data from file
    dx = pandas.DataFrame({'Dx': file['Dx']})

    # Convert 'Dx' to a factor (categorical) variable
    dx["Dx"] = dx["Dx"].str.decode("utf-8")
    dx['Dx'] = pandas.Categorical(dx['Dx'])

    sample_ids = file["sampleIDs"]
    training_probes = np.array(file["probeIDs"]).astype(str)

    tmp = np.array(case.columns)

    # Get intersecting probes
    probes = np.intersect1d(tmp, training_probes)


    # Get matching indicies
    idxs = np.where(np.isin(training_probes, probes))[0]
    print(len(probes)," overlapping CpG sites between sample and reference set. Reading training set now...")

    beta_values = file["betaValues"][:]
    beta_values = pandas.DataFrame((beta_values > 0.6)[idxs].astype(int))

    #Create a DataFrame with Dx and the filtered matrix
    ts = pandas.concat([dx["Dx"],beta_values.T], axis=1)

    new_columns = training_probes[idxs]
    new_columns = np.insert(new_columns, 0,"Dx")
    ts.columns = new_columns


    unknown_rows = case.loc[:,probes]
    unknown_rows["Dx"] = "unknown"

    file.close()
    return pandas.concat([ts, unknown_rows], ignore_index=True)


def format_beta_values(beta, scale=False, pca=True, pca_components=94):

    if scale:
        scaler = StandardScaler()
        beta = scaler.fit_transform(beta)

    if pca:
        pca = PCA(n_components=pca_components)
        beta = pca.fit_transform(beta)

    return beta



def plot_tsne(m, format=True, save_beta=False, out="cuda", use_old=False, cuml=False):

    beta = None
    if not use_old:
        standard_deviations = m.std(skipna=False)
        max_sd = standard_deviations[np.argsort(standard_deviations)][-50000:]
        beta = m.iloc[:,1:]
        beta = beta[max_sd.index.to_numpy()]

        if format:
            (beta, filename) = format_beta_values(beta, output_file=out, pca_components=70)
            out=filename
    else:
        beta = pandas.read_csv("../cuda_out/R_beta_data.csv", index_col=0)

        if format:
            (beta, filename) = format_beta_values(beta, output_file=out, pca_components=60)
            out = filename
        out+="_Rdata"

    print("Running TSNE")
    start = time()
    tsne = None
    if cuml:
        tsne = cumlTSNE(n_components=2, perplexity=30, n_iter=2500, output_type="numpy").fit_transform(beta)
        out += "_cuml"
    else:
        pass
        #tsne = TSNE(n_components=2, perplexity=30, n_iter=2500, verbose=1).fit_transform(beta)
        #tsne = cudaTSNE(n_components=2, perplexity=30, n_iter=2500).fit_transform(beta)

    done = time()

    # if tsne == None:
    #     print("Something went wrong, tsne is None")
    #     exit(0)

    print("Total elapsed time : ", done - start)
    print("TSNE : ", tsne)
    df = pandas.DataFrame({"Dx":m["Dx"], "X1": tsne[:,0], "X2": tsne[:,1]})

    timestamp = str(datetime.now()).replace(" ", "_")

    if save_beta:
        if format:
            np.savetxt(out+"_beta_"+timestamp+".csv", beta, delimiter=",")
        else:
            beta.to_csv(out+"_beta_"+timestamp+".csv")

    out += "_"+timestamp+".csv"
    print("Saving csv to "+out)
    df.to_csv(out)

output_file = "/home/sander/Documents/cudaResults/cuda"


m = pd.read_csv(snakemake.input["m"], index_col=0)
print("M: ",m["x"])
beta = pd.read_csv(snakemake.input["beta"], index_col=0)
print(beta)
print(beta.shape)

beta = format_beta_values(beta, pca_components=snakemake.config["tsne_pca_dim"])

start = time()
tsne = None
tsne = cumlTSNE(n_components=2, method='exact', perplexity=snakemake.config["tsne_perplexity"], n_iter=snakemake.config["tsne_max_iter"], output_type="numpy", random_state=42, verbose=4)
output = tsne.fit_transform(beta)
done = time()

print("Total elapsed time : ", done - start)

df = pandas.DataFrame({"Dx": m["x"], "X1": output[:, 0], "X2": output[:, 1]})


df.to_csv(snakemake.output["tsne"])

print("T-SNE iterations : ", tsne.n_iter)