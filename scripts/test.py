#####  Importing packages
import math
from collections import OrderedDict

import numpy as np
import pandas as pd
import h5py
import pyreadr
import warnings
import pickle
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
warnings.filterwarnings("ignore")


sample = "NDMA3"


##### Reading trainingset data
Probes,label = pickle.load( open( "training/"+sample+"-FeatureSelection_idf-Capper_et_al.p", "rb" ) )
fh5 = "/home/sander/Documents/static/Trainingsets/Capper_et_al.h5"

fh5=h5py.File(fh5)
label=list(fh5['Dx'][:])
label=[str(x).replace('b','').replace('\'','') for x in label]
label=[str(x).replace(',','-') for x in label]

#np.savetxt("/home/sander/Documents/TempDatasets/labels.csv",unique_labels, delimiter=",", fmt="%s", newline=",")

mapping = []
for l in label:
    if l not in mapping:
        mapping.append(l)

integer_labels = []
for l in label:
    integer_labels.append(mapping.index(l))

print(len(mapping))
print(mapping[integer_labels[662]])
print(mapping.index("A IDH- HG"))
print(mapping[integer_labels[0]])
print(fh5['Dx'][0])
en = np.bincount(np.array(integer_labels))

print(en)

print(en / len(integer_labels))
print(en / len(integer_labels) * np.log2(en / len(integer_labels)))
print(np.sum(en / len(integer_labels) * np.log2(en / len(integer_labels))))
#np.savetxt("/home/sander/Documents/TempDatasets/labels2.csv",np.array(mapping), delimiter=",", fmt="%s", newline=",")

#np.savetxt("/home/sander/Documents/TempDatasets/"+sample+"_y_train2.csv",np.array(integer_labels), delimiter=",", fmt="%d", newline=",")


probeIDs=list(fh5['probeIDs'][:])
probeIDs=[str(x).replace('b','').replace('\'','') for x in probeIDs]
train_data=fh5['betaValues']

probes_idx= [i for i in range(len(probeIDs)) if probeIDs[i] in Probes]

X_train=train_data[probes_idx,:]
X_train = np.where(X_train > 0.6, 1, 0)
X_train=pd.DataFrame(X_train.T,columns=[probeIDs[i] for i in probes_idx])

print("X_TRAIN SHAPE : ")
print(X_train.shape)



#### Reading input case data
METH_RDATA ="transformed_rdata/"+sample+"-transformed.RData"
case=pyreadr.read_r(METH_RDATA)
case=pd.DataFrame(data=list(case.items())[0][1]).T

X_test=case[X_train.columns.to_list()]

X_test.to_csv("/home/sander/Documents/TempDatasets/"+sample+"_x_test2.csv", header=False, index=False)
exit(0)
print('Train random forest...')
rf = RandomForestClassifier(criterion='entropy',  min_samples_split=4, min_samples_leaf=1,
                       n_estimators=2000,random_state=42,n_jobs= 12, verbose=1, bootstrap=False)
rf.fit(X_train, integer_labels)

pre_rf=rf.predict_proba(X_test)
print("PREDICT : ", pre_rf)
print("PREDICT : ", max(pre_rf[0]))
print("int : ", pre_rf[0].index(max(pre_rf)))
print("LABEL : ", label[pre_rf[0].index(max(pre_rf))])


exit(0)

X_train.to_csv("/home/sander/Documents/TempDatasets/"+sample+"_x_train.csv", header=False, index=False)
