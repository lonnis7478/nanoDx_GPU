#####  Importing packages
import numpy as np
import pandas as pd
import h5py
import warnings
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import LabelEncoder
import pickle
from sklearn.ensemble import RandomForestClassifier
from cuml.ensemble import RandomForestClassifier as cumlRf
import cudf
warnings.filterwarnings("ignore")

##### Reading trainingset data
Probes,label = pickle.load( open( "/home/sander/Documents/nanoDxGPU/nanoDx_GPU/training/NDMA29_SUB1-FeatureSelection_idf-Capper_et_al.p", "rb" ) )
fh5 ="/home/sander/Documents/static/Trainingsets/Capper_et_al.h5"

fh5=h5py.File(fh5)
label=list(fh5['Dx'][:])
label=[str(x).replace('b','').replace('\'','') for x in label]
probeIDs=list(fh5['probeIDs'][:])
probeIDs=[str(x).replace('b','').replace('\'','') for x in probeIDs]
train_data=fh5['betaValues']

probes_idx= [i for i in range(len(probeIDs)) if probeIDs[i] in Probes]

X_train=train_data[probes_idx,:]
X_train = np.where(X_train > 0.6, 1, 0).astype(np.float32)
X_train=pd.DataFrame(X_train.T,columns=[probeIDs[i] for i in probes_idx])

print(X_train[-800:].shape)
cv = RepeatedStratifiedKFold(n_splits=5,n_repeats=1,random_state=1)


label = np.array(label)
label_encoder = LabelEncoder()
int_labels = label_encoder.fit_transform(label)
rf = cumlRf(split_criterion=1,  min_samples_split=4, min_samples_leaf=1,
                        n_estimators=2000,random_state=42, verbose=1)

cuda_x = cudf.DataFrame.from_pandas(X_train[-800:])
cuda_y = cudf.Series(np.array(int_labels)[-800:])

print("Start training")
for train, test in cv.split(X_train, label):
    rf.fit(cuda_x, cuda_y)
    rf.predict(cuda_x)

