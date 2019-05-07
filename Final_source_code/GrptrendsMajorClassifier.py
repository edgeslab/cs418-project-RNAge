import pandas as pd
import numpy as np
import AgeGroupViz as agv
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

grptrends = agv.grp_trends
maj_clf = grptrends[['percentage_per_tissue', 'SMTS']].groupby(['SMTS']).max().reset_index()

mpl.style.use('default')
fig = plt.figure(figsize=(5,3))
ax = fig.add_subplot(111)

model_plot = sns.barplot(x = maj_clf['SMTS'],
            y=maj_clf['percentage_per_tissue'],
            data = maj_clf,
            orient='v',
            estimator=np.mean,
            capsize=0.1)

tick_labels = maj_clf['SMTS'].as_matrix()
ax.set_xticklabels(tick_labels,rotation=30, ha='right')
ax.set(ylim=(0,50))
ax.set(title='Majority classifier accuracy across tissue types')
ax.set(xlabel='Tissue types')
ax.set(ylabel='Accuracy')

plt.show()

