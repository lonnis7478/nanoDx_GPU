import numpy as np  
import pandas as pd
import h5py
import pyreadr
import logging
import collections
import warnings
import pickle
##### reading FilteredFeature
Fea1,Fea2 = pickle.load( open( snakemake.input["FilteredFeature"], "rb" ) )

Fea=Fea1+Fea2

warnings.filterwarnings("ignore")

# generate methylation training set from .h5 file
print(' Reading fh5 trainingset_meth.')

fh5 =snakemake.input["trainingset_meth"]
fh5=h5py.File(fh5)
label=list(fh5['Dx'][:])
label=[str(x).replace('b','').replace('\'','') for x in label]
probeIDs=list(fh5['probeIDs'][:])
probeIDs=[str(x).replace('b','').replace('\'','') for x in probeIDs]
train_data=fh5['betaValues']

print(' Reading METH_RDATA.')

##### load input METH_RDATA
METH_RDATA =snakemake.input["meth"]

case=pyreadr.read_r(METH_RDATA)
case=pd.DataFrame(data=list(case.items())[0][1]).T

##### overlaping CpG features with ref data
def matchfeatures(validdata,probeIDs,data):
    # overlaping CpG features with ref data
    probes=[i for i in range(len(probeIDs)) if probeIDs[i] in validdata.columns]
    probesid = [probeIDs[i] for i in probes]
    
    logging.info(len(probesid),' overlapping CpG sites between sample and reference set. Reading training set now...') 
    print(len(probesid),' overlapping CpG sites between sample and reference set. Reading training set now...') 
    
    test_data=validdata[probesid]
    
    train_data=data[probes,:]
    
    train_data = pd.DataFrame(train_data.T,columns=probesid)
    return train_data,test_data 


X_train,X_test=matchfeatures(case,probeIDs,train_data)
fh5.close()

X_train[X_train<=0.6]=0
X_train[X_train>0.6]=1   

logging.info(X_train.shape[1],' overlapping CpG sites read from training set.') 
print(X_train.shape[1],' overlapping CpG sites read from training set.')


#### feature selection 

max_CpG = snakemake.params["max_CpG"]

FeatureImportance= pickle.load( open( snakemake.input["FeatureImportance"], "rb" ) )
FeatureImportance=FeatureImportance.loc[X_train.columns]

##### maxSDs = max feature number , select max feature <= maxSDs
maxSDs=FeatureImportance.sort_values(ascending=False).index[0:min([max_CpG,X_train.shape[1]])]
X_train=X_train[maxSDs.to_list()]
X_test=X_test[maxSDs.to_list()]

Probes=X_train.columns.to_list()
pickle.dump( [Probes,label], open( snakemake.output["trainning_data"], "wb" ) )

