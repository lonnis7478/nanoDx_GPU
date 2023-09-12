import numpy as np  
import pandas as pd
import h5py
import logging
import math
import scipy.stats as stats
import pickle

import warnings

fh5 =snakemake.input["trainingset_meth"]

print(' Reading fh5 trainingset_meth.')


fh5=h5py.File(fh5)
label=list(fh5['Dx'][:])
label=[str(x).replace('b','').replace('\'','') for x in label]
probeIDs=list(fh5['probeIDs'][:])
probeIDs=[str(x).replace('b','').replace('\'','') for x in probeIDs]
train_data=fh5['betaValues']


X= train_data[:,:]
X_binary=np.where(X > 0.6, 1, 0)
X=pd.DataFrame(X_binary.T,columns=probeIDs)

def tfidf(X_,label):
    tfidf=pd.DataFrame(columns=X_.columns)
    for tumor in list(set(label)):
        tumor_idx=[i  for i in range(len(label)) if label[i]==tumor]
        tumor_X=X_.T.filter(tumor_idx)
        TF=tumor_X.sum(axis=1) / len(tumor_idx)
        
        IDF=(len(X_.index) / (X_.sum()+1)).apply(lambda x: math.log(x))
        
        TFIDF= TF * IDF
        tfidf.loc[tumor]=TFIDF.tolist()
    return tfidf

if snakemake.params["method"]=='tfidf':
    tfidf=tfidf(X,label)
    X_tfidf_scale=tfidf.T.select_dtypes(include='number').apply(stats.zscore).T
    pickle.dump( X_tfidf_scale.max().sort_values(ascending=False), open( snakemake.output["FeatureImportance"], "wb" ) )
    print('TF_IDF Done')
elif snakemake.params["method"]=='std':
    sds=X.std()
    pickle.dump( sds, open( snakemake.output["FeatureImportance"], "wb" ) )
sum_X=X.sum(axis=0)
all0=sum_X[sum_X==0].index.tolist()
all1=sum_X[sum_X==1].index.tolist()
pickle.dump( [all0,all1], open( snakemake.output["FilteredFeature"], "wb" ) )