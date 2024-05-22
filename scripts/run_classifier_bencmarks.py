import os
import sys

import numpy as np
import pandas as pd

import h5py
import warnings
from pyreadr import pyreadr
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
import pickle
from time import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import brewer2mpl

np.set_printoptions(suppress=True)
from cudaclassifier import CudaClassifier


categories = ["GBM, G34","DMG, K27","ATRT, SHH","CONTR, CEBM","GBM, MYCN","LGG, PA PF", "CONTR, REACT","LGG, PA MID","LGG, RGNT","MB, WNT",
              "ATRT, MYC", "ATRT, TYR","LGG, PA/GG ST","LGG, SEGA","MB, G4","MB, SHH INF", "MB, SHH CHL AD","MB, G3","EPN, RELA","PIN T,  PB A","SUBEPN, PF",
              "PTPR, B","SUBEPN, ST","EPN, YAP","EPN, PF A","EPN, PF B","EFT, CIC", "ETMR","CNS NB, FOXR2","LYMPHO","HGNET, BCOR","GBM, MES","GBM, RTK II",
              "GBM, RTK I","LGG, MYB","CONTR, HEMI","HGNET, MN1","GBM, MID","HMB", "EPN, MPE","SCHW","O IDH","A IDH, HG","MNG","A IDH","LGG, GG","CN",
              "SUBEPN, SPINE","PIN T,  PB B","PXA","ANA PA","CONTR, INFLAM", "PIN T, PPT","CPH, ADM","CPH, PAP","ENB, A","PITUI","CONTR, PINEAL", "PGG, nC",
              "LGG, DNT","CHGL","MELAN","PLEX, AD","ENB, B","LIPN", "EPN, SPINE","PTPR, A","SFT HMPC","PLEX, PED A","GBM, RTK III","IHG", "MELCYT","DLGNT",
              "PITAD, ACTH","PITAD, STH DNS B","PITAD, PRL", "PITAD, FSH LH","PLEX, PED B","EWS","SCHW, MEL","CONTR, ADENOPIT", "LGG, DIG/DIA","PITAD, STH SPA",
              "PITAD, STH DNS A","CONTR, WM","PLASMA", "CHORDM","RETB","CONTR, PONS","CONTR, HYPTHAL","PITAD, TSH", "bruh"]


#['LYMPHO', 'GBM, MES', 'MB, SHH CHL AD', 'ENB, B', 'LGG, PA PF', 'PIN T, PPT', 'GBM, RTK II', 'MB, WNT', 'PIN T,  PB B', 'HGNET, BCOR', 'ANA PA', 'PTPR, A', 'EPN, SPINE', 'LGG, DNT', 'MB, G3', 'A IDH', 'MNG', 'GBM, RTK I', 'CNS NB, FOXR2', 'ATRT, MYC', 'CONTR, PINEAL', 'EFT, CIC', 'PITUI', 'PTPR, B', 'LGG, DIG/DIA', 'PITAD, TSH', 'PLEX, PED B', 'ATRT, SHH', 'MB, SHH INF', 'CONTR, HEMI', 'MB, G4', 'MELAN', 'DMG, K27', 'GBM, MYCN', 'ENB, A', 'LGG, PA MID', 'ETMR', 'PLEX, AD', 'LIPN', 'GBM, MID', 'EPN, PF A', 'PLEX, PED A', 'HGNET, MN1', 'EPN, RELA', 'EPN, MPE', 'HMB', 'CPH, ADM', 'PITAD, FSH LH', 'PITAD, STH SPA', 'LGG, GG', 'CN', 'MELCYT', 'CONTR, CEBM', 'PXA', 'SUBEPN, PF', 'RETB', 'A IDH, HG', 'GBM, G34', 'O IDH']


warnings.filterwarnings("ignore")

def load_sample(sample_name, trainingset="Capper_et_al"):
    data = f"../training/{sample_name}-FeatureSelection_idf-{trainingset}.p",
    trainingset_meth = f"/home/sander/Documents/static/Trainingsets/{trainingset}.h5"
    meth = f"../transformed_rdata/{sample_name}-transformed.RData"


    ##### Reading trainingset data
    Probes,label = pickle.load( open( data, "rb" ) )
    fh5 =trainingset_meth

    fh5=h5py.File(fh5)
    label=list(fh5['Dx'][:])
    label=[str(x).replace('b','').replace('\'','').replace('\'', '"') for x in label]

    encoded = []

    for i in range(len(label)):
        encoded.append(categories.index(label[i]))


    probeIDs=list(fh5['probeIDs'][:])
    probeIDs=[str(x).replace('b','').replace('\'','') for x in probeIDs]
    train_data=fh5['betaValues']

    probes_idx= [i for i in range(len(probeIDs)) if probeIDs[i] in Probes]

    X_train=train_data[probes_idx,:]
    X_train = np.where(X_train > 0.6, 1, 0)
    X_train=pd.DataFrame(X_train.T,columns=[probeIDs[i] for i in probes_idx])


    #
    # #### Reading input case data
    METH_RDATA = meth
    case=pyreadr.read_r(METH_RDATA)
    case=pd.DataFrame(data=list(case.items())[0][1]).T

    X_test=case[X_train.columns.to_list()]


    y = np.array(encoded)
    X_test = X_test.values.reshape((X_train.shape[1],))
    X_train = X_train.values


    folder = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=42)

    samples = np.arange(X_train.shape[0])

    folds = np.array([])
    training_folds = np.array([])
    validation_folds = np.array([])

    counter = 0
    np.set_printoptions(suppress=True, precision=6)
    bruh = []

    for train_index, test_index in folder.split(X_train, y):
        folds = np.append(folds, samples[train_index])
        training_folds = np.append(training_folds, samples[train_index])

        folds = np.append(folds, samples[test_index])
        validation_folds = np.append(validation_folds, samples[test_index])

    return X_train, y, X_test, folds, training_folds, validation_folds


base_folder = "/home/sander/Documents/NDMA_dorado_basecall"
max_samples = sys.argv[1]

if max_samples is None:
    max_samples = 10



counter = 0
for file in os.listdir(base_folder):

    if counter > max_samples:
        print("Reached maximum samples, terminating...")
        exit(0)

    print(f"Running classification on sample {file}")
    X_train, y, X_test, folds, training_folds, validation_folds = load_sample(file)

    cl = CudaClassifier(n_ensemble=2000)
    cl.print_cuda_properties()

    start = time()
    result = cl.fit(X_train.flatten().astype(np.intc), y.astype(np.intc), X_test.astype(np.intc).flatten(),
                    folds.astype(np.intc), training_folds.astype(np.intc), validation_folds.astype(np.intc),
                    X_train.shape[0], X_train.shape[1], 42)
    done = time()

    print("Time spent : ", done - start, "s")

    votes = pd.DataFrame({"Dx": categories, "Var1": np.arange(len(categories)),
                          "Freq": [x / cl.get_forest_size() for x in result.predictions],
                          "cal_Freq": result.calibrated_proba})
    votes = votes.sort_values(by=['Freq'], ascending=False)

    counter += 1