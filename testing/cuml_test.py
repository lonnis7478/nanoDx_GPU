from cuml.ensemble import RandomForestClassifier as cuRF
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn.model_selection import cross_val_score
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from time import time

def print_image(row, df):
    temp=df.iloc[row,:].values
    temp = temp.reshape(28,28).astype('uint8')
    plt.imshow(temp)
    plt.show()

train = pd.read_csv("/home/sander/Documents/train.csv")
test = pd.read_csv("/home/sander/Documents/test.csv")

print(train.head(10))

df_x = train.iloc[:100,1:]
df_y = train.iloc[:100,0:]

print(df_x.shape)

X_train, X_test, Y_train, Y_test = train_test_split(df_x, df_y, test_size=0.2, random_state=0)

start = time()
rf = RandomForestClassifier(n_estimators=100, verbose=1, warm_start=True)
rf.fit(X_train, Y_train)
end = time()

prediction = rf.predict(X_test)


score = cross_val_score(rf, df_x, df_y)
print("Cross validation score : ")
print(np.mean(score))

print("Total time used : ", end - start)



