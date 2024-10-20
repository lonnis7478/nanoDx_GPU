import h5py
import numpy as np
import pandas
import pyreadr

import pandas as pd
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
#from tsnecuda import TSNE as cudaTSNE
#from cuml.manifold import TSNE as cumlTSNE

from time import time
from datetime import datetime


##### Reading trainingset data

fh5=h5py.File("/home/sander/Documents/static/Trainingsets/Capper_et_al.h5")

case=pyreadr.read_r("../transformed_rdata/NDMA1-transformed.RData")
case=pd.DataFrame(data=list(case.items())[0][1]).T

print(case)

label=list(fh5['Dx'][:])
label=[str(x).replace('b','').replace('\'','') for x in label]
probeIDs=list(fh5['probeIDs'][:])
probeIDs=[str(x).replace('b','').replace('\'','') for x in probeIDs]
train_data=fh5['betaValues']
Probes = list(set(case.columns).intersection(probeIDs))
idxs = [list(probeIDs).index(p) for p in Probes]



print(len(probeIDs))
print(len(Probes))
print(len(idxs))

Dx = pd.Categorical(pd.Series(fh5['Dx'][:]))

ts = pd.DataFrame({
    'Dx': Dx,
    **{Probes[idx]: (fh5["betaValues"][:, idx]).astype(int) for idx in idxs}
})

print(ts)

exit(0)
X_test=case[X_train.columns.to_list()]




X_train=train_data[probes_idx,:]
X_train = np.where(X_train > 0.6, 1, 0)
X_train=pd.DataFrame(X_train.T,columns=[probeIDs[i] for i in probes_idx])





def condition(row) :
    if(row['methylated_frequency'] > 0.6) :
        return 1
    return 0

def format_beta_values(beta, scale=False, pca=True, pca_components=94):

    if scale:
        scaler = StandardScaler()
        beta = scaler.fit_transform(beta)

    if pca:
        pca = PCA(n_components=pca_components)
        beta = pca.fit_transform(beta)

    return beta



output_file = "/home/sander/Documents/cudaResults/cuda"


m = pd.read_csv("/home/sander/Documents/nanoDxGPU/nanoDx_GPU/betaValues/NDMA1.Capper_et_al.m.csv", index_col=0)
print("M: ",m["x"])
beta = pd.read_csv("/home/sander/Documents/nanoDxGPU/nanoDx_GPU/betaValues/NDMA1.Capper_et_al.beta.csv", index_col=0)
print(beta)
print(beta.shape)

beta = format_beta_values(beta, pca_components=94)

start = time()
tsne = None
tsne = cumlTSNE(n_components=2, perplexity=25, n_iter=1500, output_type="numpy", random_state=42, verbose=4, method='exact')
output = tsne.fit_transform(beta)
done = time()

print("Total elapsed time : ", done - start)

df = pandas.DataFrame({"Dx": m["x"], "X1": output[:, 0], "X2": output[:, 1]})


df.to_csv("temp_tsne_pca.csv")

print("T-SNE iterations : ", tsne.n_iter)