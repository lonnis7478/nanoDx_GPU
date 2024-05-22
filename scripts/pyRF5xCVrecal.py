#####  Importing packages
import numpy as np  
import pandas as pd
import h5py
import pyreadr
from sklearn.ensemble import RandomForestClassifier
import logging
import collections
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import brewer2mpl
import warnings
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.calibration import CalibratedClassifierCV
import pickle
from time import time
warnings.filterwarnings("ignore")

##### Reading trainingset data
Probes,label = pickle.load( open( snakemake.input["data"], "rb" ) )
fh5 =snakemake.input["trainingset_meth"]

fh5=h5py.File(fh5)
label=list(fh5['Dx'][:])
label=[str(x).replace('b','').replace('\'','') for x in label]
probeIDs=list(fh5['probeIDs'][:])
probeIDs=[str(x).replace('b','').replace('\'','') for x in probeIDs]
train_data=fh5['betaValues']

probes_idx= [i for i in range(len(probeIDs)) if probeIDs[i] in Probes]

X_train=train_data[probes_idx,:]
X_train = np.where(X_train > 0.6, 1, 0)
X_train=pd.DataFrame(X_train.T,columns=[probeIDs[i] for i in probes_idx])

#### Reading input case data
METH_RDATA =snakemake.input["meth"]
case=pyreadr.read_r(METH_RDATA)
case=pd.DataFrame(data=list(case.items())[0][1]).T

X_test=case[X_train.columns.to_list()]


##### random forest classifier and get the votes
print('Train random forest...')


rf = RandomForestClassifier(criterion='entropy',  min_samples_split=4, min_samples_leaf=1,
                       n_estimators=2000,random_state=42,n_jobs= 12,oob_score=True, verbose=1, max_depth=8)
rf.fit(X_train, label)


pre_rf=rf.predict(X_test)

x_score=max(rf.predict_proba(X_test)[0])
pre_class=pre_rf[0]

##### predict case and get the votes from each trees

def pre_singletree(est,X_input,classes_):
    res=int(est.predict(X_input))
    return classes_[res],res

predictions_trees=list(map(lambda i:pre_singletree(i, X_test,rf.classes_), rf.estimators_))


tmp = [i[0] for i in predictions_trees]
rownames=[i[0] for i in collections.Counter(tmp).items()]

Freq=[i[1] for i in collections.Counter(tmp).items()]
Var1= [i[1] for i in predictions_trees]

votes = pd.DataFrame(list(zip(Var1, Freq)),
               columns =['Var1', 'Freq'], index= rownames)
votes['Freq']=votes['Freq']/sum(votes['Freq']) *100
votes=votes.sort_values(by=['Freq'])
votes['Dx'] = votes.index # class labels

#####  5-fold CV calibration model to get recalibration scores
print('Running 5-fold CV...')
def get_proba_CalibratedClassifierCV(X_train,y_train):
    cv = RepeatedStratifiedKFold(n_splits=5,n_repeats=1,random_state=1)
    clf = RandomForestClassifier(criterion='entropy',  min_samples_split=4, min_samples_leaf=1,
                       n_estimators=200,random_state=42,n_jobs= 12,oob_score=True, verbose=1)

    clf_sigmoid = CalibratedClassifierCV(clf, cv=cv, method='sigmoid',n_jobs =-1,ensemble=False)

    start = time()
    clf_sigmoid.fit(X_train, y_train)
    done = time()
    print("Calibrated classifier time : ", done - start)

    start = time()
    res_proba=clf_sigmoid.predict_proba(X_test)
    done = time()
    print("Calibrated classifier prediction time : ", done -start)
    proba = pd.DataFrame(data={'cal_Freq': res_proba[0]} , index=clf_sigmoid.classes_ )
    return proba


proba=get_proba_CalibratedClassifierCV(X_train,label)
proba['Dx']=proba.index

df=proba.merge(votes,how='left', left_on='Dx', right_on='Dx').fillna(0)

##### write results to txt file
report_content=['Number of features: '+str(rf.n_features_in_),
                "Predicted Class: "+pre_rf[0],
                "Initial Score: "+ str(x_score),
                "Calibrated Score: "+ str( proba.loc[pre_class].values[0])
               ]
report="\n".join(report_content)

with open(snakemake.output["txt"], 'w') as f:
    f.write(report)

##### output votes and RF model info
pyreadr.write_rdata(snakemake.output["votes"], df, df_name="votes")
pyreadr.write_rdata(snakemake.output["model_info"], pd.DataFrame({'oob.error': [1 - rf.oob_score_], 'num.features': [rf.n_features_in_]}), df_name="model_info")

##### plot piechart
with PdfPages(snakemake.output["pdf"]) as pdf:
    plt.subplot(2, 1, 1)
    plt.title(snakemake.wildcards["sample"])
    plt.pie(votes['Freq'],labels=list(votes.index),startangle=90,colors=brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors)
    plt.subplot(2, 1, 2)
    plt.pie(df['cal_Freq'],labels=list(df.index),startangle=90,colors=brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors)
    pdf.savefig()
    plt.close()
