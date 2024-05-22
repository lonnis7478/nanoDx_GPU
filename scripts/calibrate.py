from cuml.ensemble import RandomForestClassifier
from cuml import LogisticRegression
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import RepeatedStratifiedKFold
#from sklearn.ensemble import RandomForestClassifier
from time import time


iris = datasets.load_breast_cancer()
X = iris.data
y = iris.target
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


n_classes = 200
n_fold = 5
rf = RandomForestClassifier()

start = time()
rf.fit(X_train, y_train)

print("Shape :", X_train.shape)
prediction = rf.predict_proba(X_test)[:1]

calibration_layer = []

for i in range(n_classes):
    calibration_layer.append(LogisticRegression())


for i in range(n_classes):
    for k in range(n_fold):
        calibration_layer[i].fit(X_train, y_train)

end = time()

print("Time used : ", end - start)



print("bruh")
exit(0)



cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=1)
clf = RandomForestClassifier(criterion="entropy", min_samples_split=4, min_samples_leaf=1,
                             n_estimators=2000, random_state=42, verbose=1)
clf.fit(X_train, y_train)

clf_sigmoid = CalibratedClassifierCV(clf, cv=cv, method='sigmoid', ensemble=False)

clf_sigmoid.fit(X_train, y_train)
