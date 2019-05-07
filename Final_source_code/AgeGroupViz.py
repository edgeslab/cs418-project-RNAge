import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import display, HTML,Image
from pathlib import Path
import sklearn
import sklearn.model_selection
import sklearn.feature_selection
from os import listdir

import Init as init

import warnings
warnings.filterwarnings("ignore")

meta = init.meta

counts=pd.DataFrame(meta['SMTS'].value_counts())
df=meta[meta['SMTS'].isin(counts[counts['SMTS']>1].index)]
df=pd.crosstab(index=df['SMTS'],columns=df['AGE'])
sorted_i=df.sum(axis=1).sort_values(ascending=False).index.values
df=df.loc[sorted_i,:]

import matplotlib as mpl
from matplotlib.lines import Line2D
    
df2=meta[meta['SMTS'].isin(counts[counts['SMTS']>200].index)]
df2_cols = df2.columns.values
df2_cols[0] = 'row_id'
df2.columns = df2_cols
grp_trends = df2[['SMTS', 'AGE', 'row_id']].groupby(['SMTS', 'AGE']).count().reset_index()
grp_trends['percentage_per_tissue'] = 1.0
for index, row in grp_trends.iterrows():
    grp_trends.at[index, 'percentage_per_tissue'] = (row['row_id'] / grp_trends['row_id'][grp_trends['SMTS']==row['SMTS']].sum()) * 100

unique_age = grp_trends['AGE'].unique()
unique_tissue = grp_trends['SMTS'].unique()

percent_matrix = np.zeros(shape=(len(unique_age), len(unique_tissue)))

i = 0
for entry in unique_age:
    percent_matrix[i] = (grp_trends['percentage_per_tissue'][grp_trends['AGE']==entry]).as_matrix()
    i += 1

cmap = mpl.cm.get_cmap('hsv', len(unique_tissue)*2)
colorPalette = []
for i in range(cmap.N):
    rgb = cmap(i)[:3] #returns [rgba], hence extracting [rgb]
    colorPalette.append(mpl.colors.rgb2hex(rgb))

mpl.style.use('default')
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)

configs = percent_matrix[0]
N = len(configs)
ind = np.arange(N)
width = 0.4

i = 0
while i < len(unique_age):
    j = i-1
    res = np.zeros(shape=(1, len(unique_tissue)))
    while j >= 0:
        res += percent_matrix[j]
        j -= 1
    plt.bar(ind, percent_matrix[i], width, bottom=res[0], color=colorPalette[i*4-1], tick_label='placeholder')
    i += 1
    
ax.set_xticklabels(unique_tissue, rotation=90)
ax.set(title='Distribution of age groups per tissue type')
ax.set(xlabel='Tissue types')
ax.set(ylabel='Perentage of data points per tissue type')
ax.set(ylim=(0,100))
ax.set(yticks=np.arange(0,110,10))

leg_ele = [
    Line2D([0], [0], lw=2, color=colorPalette[0*4-1], label=unique_age[0]),
    Line2D([0], [0], lw=2, color=colorPalette[1*4-1], label=unique_age[1]),
    Line2D([0], [0], lw=2, color=colorPalette[2*4-1], label=unique_age[2]),
    Line2D([0], [0], lw=2, color=colorPalette[3*4-1], label=unique_age[3]),
    Line2D([0], [0], lw=2, color=colorPalette[4*4-1], label=unique_age[4]),
    Line2D([0], [0], lw=2, color=colorPalette[5*4-1], label=unique_age[5])
    ]
ax.legend(handles = leg_ele, loc='upper left', bbox_to_anchor=(1,1), title='Age groups')
       
plt.show()