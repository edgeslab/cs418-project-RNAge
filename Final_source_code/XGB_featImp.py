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
    uniprot=pd.read_csv("annotation/uniprot_meta.tsv",sep="\t",names=["Uniprot","HGNC_Species","Full HGNC","Location","Type","Desc"])
    biomart=pd.read_csv("annotation/2018-04-12_biomart_ensemblGene_hgnc_uniprot.tab",sep="\t")
    sel=pd.isnull(biomart['hgnc'])
    biomart.loc[sel,"hgnc"]=biomart.loc[sel,"ensembl"] # fills in missing values of hgnc with ensembl
    biomart = biomart.drop_duplicates(subset='ensembl', keep='first')
    biomart=biomart.set_index(biomart['ensembl'])
    res=[]
    for x in range(len(results.index.values)):
        alg=algs[x]
        TISSUE=results.index.values[x]
        feat_imp = pd.Series(alg.get_booster().get_score(importance_type='weight')).sort_values(ascending=False)[:200]
        feat_imp=feat_imp.to_frame('Score').join(biomart)
        sel=pd.isnull(feat_imp['hgnc'])
        feat_imp.loc[sel,"hgnc"]=feat_imp.index.values[sel]
        feat_imp=feat_imp.set_index(feat_imp['hgnc'])
        feat_imp['Score'][:10].plot(kind='bar', title='Feature Importances for {}'.format(TISSUE))
        plt.savefig(fname='plots/{}_feat_imp.png'.format(TISSUE),bbox_inches='tight')
        plt.close()
    for x in range(len(results.index.values)):
        alg=algs[x]
        TISSUE=results.index.values[x]
        feat_imp = pd.Series(alg.get_booster().get_score(importance_type='weight')).sort_values(ascending=False)[:200]
        feat_imp=feat_imp.to_frame('Score').join(biomart)
        sel=pd.isnull(feat_imp['hgnc'])
        feat_imp.loc[sel,"hgnc"]=feat_imp.index.values[sel]
        feat_imp=feat_imp.set_index(feat_imp['hgnc'])
        feat_imp['Score'][:10].plot(kind='bar', title='Feature Importances for {}'.format(TISSUE))
        plt.show()
        plt.close()
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