import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import display, HTML
from pathlib import Path
from os import listdir
import sklearn.model_selection
import sklearn.feature_selection
sns.set_style("whitegrid")
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn.model_selection import cross_validate,GridSearchCV
from sklearn import metrics #Additional scklearn functions
from sklearn.externals import joblib
from xgboost import plot_tree

def plotXGB(algs,results):
    plt.figure(figsize=(8, 4))
    ax=sns.barplot(x=results.index.values,y=results['Score'])
    x=ax.set_title("XGB Model Accuracies for All Tissues")
    x=ax.set_xlabel("Tissues")
    x=ax.set_ylabel("Accuracy Score")
    x=ax.set_xticklabels(labels=results.index.values,rotation=38)
    plt.show(ax)
    
if __name__ == "__main__":
    data_dir=Path("data")
    tissue_dir=Path("tissue-specific")
    manifest={"data":"All_Tissue_Site_Details.combined.reads.gct",
              "sample_meta":"GTEx_v7_Annotations_SampleAttributesDS.txt",
              "subject_meta":"GTEx_v7_Annotations_SubjectPhenotypesDS.txt",
               "merged_meta":"merged_meta.tsv"}
    meta=pd.read_csv(data_dir/manifest['merged_meta'],sep="\t",dtype={'SMUBRID':object,'SEX':object,'DTHHRDY':object})
    meta=meta[~(meta['AGE'].isnull())]
    counts=pd.DataFrame(meta['SMTS'].value_counts())
    df=meta[meta['SMTS'].isin(counts[counts['SMTS']>200].index)]
    results=pd.read_csv("models/results.tsv",sep="\t",index_col=0)

    algs=[joblib.load('models/{}_model.pkl'.format(results.index.values[x])) for x in range(len(results))]
    plotXGB(algs,results)