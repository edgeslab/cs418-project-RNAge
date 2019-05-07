#Import necessary packages
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
from sklearn.decomposition import PCA

fig, axs=plt.subplots(1,2,squeeze=False,figsize=(15,5),sharex=True)
accuracy_SVMM=pd.read_csv('SVM_ACCURACY_FILE.csv')

print(fig)
f=sns.barplot(accuracy_SVMM['tissue'], accuracy_SVMM['accuracy'], ax=axs[0][0])#.set_title('Accuracy of Different Tissues for SVM method')
f.set_title('Accuracy of Different Tissues for SVM method')
f.set_xlabel('Tissue Type')
f.set_ylabel('Accuracy')
for tick in axs[0][0].get_xticklabels():
    tick.set_rotation(40)
    tick.set_ha('right')
accuracy_KNNN=pd.read_csv('KNN_ACCURACY_FILE.csv')

sns.barplot(accuracy_KNNN['tissue'], accuracy_KNNN['accuracy'],ax=axs[0][1])
plt.title('Accuracy of Different Tissues for KNN method')
plt.xlabel('Tissue Type')
plt.ylabel('Accuracy')
plt.xticks(rotation=40, ha='right')

plt.show()
