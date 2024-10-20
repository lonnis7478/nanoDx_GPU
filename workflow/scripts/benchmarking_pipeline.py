import psycopg2
import re
import pyreadr
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

class Benchmarker:

    def __init__(self, input_folder, output_folder):

        self.conn = psycopg2.connect(
            dbname="master",
            user="postgres",
            password="nQgartWZ",
            host="192.168.0.133",
            port=5432
        )

        self.cur = self.conn.cursor()
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.rules = set()

    def random_forest_benchmark(self, filename, sample_name, type="CUDA", n_estimators=2000, max_depth=8, max_features=96):
        X_train, y, label, X_test, folds, training_folds, validation_folds = self.load_sample(filename)
        sample = self.get_existing_sample(sample_name, X_train.shape[1])

        run_config = [n_estimators, max_depth, max_features]

        if type == "CUDA":

            classifier = CudaClassifier(n_ensemble=n_estimators, max_depth=max_depth)
            classifier.print_cuda_properties()

            start = time()
            result = classifier.fit(X_train.flatten().astype(np.intc), y.astype(np.intc), X_test.astype(np.intc).flatten(),
                                    folds.astype(np.intc), training_folds.astype(np.intc), validation_folds.astype(np.intc),
                                    X_train.shape[0], X_train.shape[1], 42)
            done = time()

            votes = pd.DataFrame({"Dx": categories, "Var1": np.arange(len(categories)),
                                  "Freq": [x / classifier.get_forest_size() for x in result.predictions],
                                  "cal_Freq": result.calibrated_proba})
            votes = votes.sort_values(by=['Freq'], ascending=False)

            self.insert_classification(sample["id"], type, done - start, votes)

            print(f"Classification for sample {sample['name']} done")
            print(f"[Time : {done-start}s, Version : CUDA]")

        elif type == "SKLEARN":
            rf = RandomForestClassifier(criterion='entropy', min_samples_split=4, min_samples_leaf=1, max_features=max_features,
                                        n_estimators=n_estimators, random_state=42, n_jobs=12, oob_score=True, verbose=1,
                                        max_depth=max_depth)
            start = time()
            rf.fit(X_train, label)
            done = time()

            print(f"Classification for sample {sample['name']} done")
            print(f"[Time : {done - start}s, Version : SKLEARN]")

    def random_forest_benchmark_batch(self, base_folder, type="CUDA", batch_size=20):
        counter = 0
        for file in os.listdir(base_folder):

            if counter > batch_size:
                print("Reached maximum samples, terminating...")
                exit(0)

            if type == "ALL":
                self.random_forest_benchmark(file, "", "CUDA")
                self.random_forest_benchmark(file, "", "SKLEARN")
            else:
                self.random_forest_benchmark(file, "", type)

            counter += 1

    def post_run_benchmark(self, sample_name):

        if len(self.rules) == 0:
            for file in os.listdir("../benchmarks"):
                if "fast5" not in file and sample_name + "." in file:
                    file = file.replace(".Capper_et_al", "")
                    rule = re.findall(f"^{sample_name}.(.*)\.benchmark.txt$", file)
                    if len(rule) > 0:
                        self.rules.add(rule[0])

        model_info = pyreadr.read_r(f"../classification/{sample_name}-model_info-RF5xCVrecal-Capper_et_al.RData")[
            "model_info"]

        cpgs = model_info["num.features"].values[0]
        sample = self.get_sample(sample_name, cpgs)

        print(sample["name"].values[0])

        votes = pyreadr.read_r(f"../classification/{sample['name'].values[0]}-votes-RF5xCVrecal-Capper_et_al.RData")[
            "votes"]
        cuda_votes = \
        pyreadr.read_r(f"../classification/{sample['name'].values[0]}-votes-CudaClassifier-Capper_et_al.RData")["votes"]

        rule_df = self.build_rule_dataframe(sample_name, self.rules)

        self.insert_classification(sample['id'].values[0], "CUDA", rule_df[rule_df["name"] == "CUDA_classifier"]["time"].values[0], cuda_votes)
        self.insert_classification(sample['id'].values[0], "SKLEARN", rule_df[rule_df["name"] == "RF5xCVrecal"]["time"].values[0], votes)
        self.update_pipeline_benchmarks(rule_df, sample)


    def get_existing_sample(self, sample_name, cpgs):
        return pd.read_sql_query(f"SELECT * FROM samples WHERE name = '{sample_name}' AND cpgs = {cpgs};", self.conn)


    def insert_classification(self, sample_id, type, time, votes):
        votes = votes.sort_values(by="Freq", ascending=False)

        query = f"INSERT INTO classifications (sample_id, type, time, p1, p2, p3, p4, p5) VALUES ({sample_id}, '{type}', {time}, {categories.index(votes['Dx'][0])}, {categories.index(votes['Dx'][1])}, {categories.index(votes['Dx'][2])}, {categories.index(votes['Dx'][3])}, {categories.index(votes['Dx'][4])});"
        self.cur.execute(query)
        self.conn.commit()

    def update_pipeline_benchmarks(self, rules, sample):
        benchmarks = pd.read_sql("SELECT * FROM pipeline_benchmarks;", self.conn)

        filtered_df = rules[rules["name"].isin(benchmarks.columns)]
        unpacked = ','.join('"'+filtered_df["name"]+'"')
        unpacked_filtered = ",".join(filtered_df["time"].astype(str))


        query = f"INSERT INTO pipeline_benchmarks({unpacked}, sample_id) VALUES({unpacked_filtered}, {sample['id'].values[0]});"
        self.cur.execute(query)
        self.conn.commit()

    def build_rule_dataframe(self, sample, rules):
        df = pd.DataFrame(columns=["name", "time"])

        for rule in rules:
            for file in os.listdir("../benchmarks"):
                if rule +"." in file and sample +"." in file:
                    print(file)
                    stats = pd.read_csv("../benchmarks/" + file, sep="\t")
                    df.loc[len(df)] = [rule, stats["s"][0]]

        return df

    def get_sample(self, sample_name, cpgs):
        sample = self.get_existing_sample(self.conn, sample_name, cpgs)

        if sample.empty:
            self.cur.execute(f"INSERT INTO samples(name, cpgs) VALUES('{sample_name}', {cpgs}) RETURNING *;")
            self.conn.commit()
            id, name, cpgs = self.cur.fetchone()
            sample = self.get_existing_sample(self.conn, name, cpgs)

        return sample

    def load_sample(self, sample_name, trainingset="Capper_et_al"):
        data = f"../training/{sample_name}-FeatureSelection_idf-{trainingset}.p",
        trainingset_meth = f"/home/sander/Documents/static/Trainingsets/{trainingset}.h5"
        meth = f"../transformed_rdata/{sample_name}-transformed.RData"

        ##### Reading trainingset data
        Probes, label = pickle.load(open(data, "rb"))
        fh5 = trainingset_meth

        fh5 = h5py.File(fh5)
        label = list(fh5['Dx'][:])
        label = [str(x).replace('b', '').replace('\'', '').replace('\'', '"') for x in label]

        encoded = []

        for i in range(len(label)):
            encoded.append(categories.index(label[i]))

        probeIDs = list(fh5['probeIDs'][:])
        probeIDs = [str(x).replace('b', '').replace('\'', '') for x in probeIDs]
        train_data = fh5['betaValues']

        probes_idx = [i for i in range(len(probeIDs)) if probeIDs[i] in Probes]

        X_train = train_data[probes_idx, :]
        X_train = np.where(X_train > 0.6, 1, 0)
        X_train = pd.DataFrame(X_train.T, columns=[probeIDs[i] for i in probes_idx])

        #
        # #### Reading input case data
        METH_RDATA = meth
        case = pyreadr.read_r(METH_RDATA)
        case = pd.DataFrame(data=list(case.items())[0][1]).T

        X_test = case[X_train.columns.to_list()]

        y = np.array(encoded)
        X_test = X_test.values.reshape((X_train.shape[1],))
        X_train = X_train.values

        folder = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=42)

        samples = np.arange(X_train.shape[0])

        folds = np.array([])
        training_folds = np.array([])
        validation_folds = np.array([])

        np.set_printoptions(suppress=True, precision=6)

        for train_index, test_index in folder.split(X_train, y):
            folds = np.append(folds, samples[train_index])
            training_folds = np.append(training_folds, samples[train_index])

            folds = np.append(folds, samples[test_index])
            validation_folds = np.append(validation_folds, samples[test_index])

        return X_train, y, label, X_test, folds, training_folds, validation_folds



benchmarker = Benchmarker("","")
benchmarker.random_forest_benchmark_batch("/home/sander/Documents/NDMA_dorado_basecall", "CUDA", 100)











