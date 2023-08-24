from time import time

import numpy as np
import pandas as pd


# For plotting
from matplotlib import offsetbox
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import seaborn as sns
import plotly.graph_objects as go
from sklearn.manifold import TSNE
from tsnecuda import TSNE as cudaTSNE
from sklearn.preprocessing import StandardScaler
import timeit

sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})


train = pd.read_csv('../../datasets/MNIST_CSV/mnist_train.csv', header=None)
test = pd.read_csv('../../datasets/MNIST_CSV/mnist_test.csv', header=None)


y = train.loc[:,0].values
x = train.loc[:,1:].values
print(x.shape)

standardized_data = StandardScaler().fit_transform(x)
print(standardized_data.shape)

x_subset = x[0:10000]
y_subset = y[0:10000]

start = time()
#tsne = TSNE(n_components=2, perplexity=40, n_iter=300).fit_transform(x_subset)
tsne = cudaTSNE(n_components=2, perplexity=40, n_iter=300).fit_transform(x_subset)
done = time()

print("Total elapsed time : ", done - start)

plt.scatter(tsne[:, 0], tsne[:, 1], s= 5, c=y_subset, cmap='Spectral')
plt.gca().set_aspect('equal', 'datalim')
plt.colorbar(boundaries=np.arange(11)-0.5).set_ticks(np.arange(10))
plt.title('Visualizing Kannada MNIST through t-SNE', fontsize=24);
plt.show()