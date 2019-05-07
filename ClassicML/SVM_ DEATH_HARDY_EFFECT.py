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
%matplotlib inline
from sklearn.decomposition import PCA

#Setting up the environment and importing the necessary data

sns.set_style("darkgrid")
data_dir=Path("data")
tissue_dir=Path("tissue-specific")

#!mkdir data && cp merged_meta.tsv data #Needed after cloning repo
manifest={
    "data":"All_Tissue_Site_Details.combined.reads.gct",
    "sample_meta":"GTEx_v7_Annotations_SampleAttributesDS.txt",
    "subject_meta":"GTEx_v7_Annotations_SubjectPhenotypesDS.txt",
    "merged_meta":"merged_meta.tsv"}

meta=pd.read_csv(manifest['merged_meta'],sep="\t",dtype={'SMUBRID':object,'SEX':object,'DTHHRDY':object})
meta

meta=meta[~(meta['DTHHRDY'].isnull())] # removes all samples without age
meta=meta[~(meta['DTHHRDY']=='3.0')]
meta=meta[~(meta['DTHHRDY']=='4.0')]

meta['DTHHRDY']=meta['DTHHRDY'].map({'1.0': 'Fast', '2.0':'Fast', '0.0':'Vent'})

df=meta[meta['SMTS'].isin(counts[counts['SMTS']>200].index)]
df=pd.crosstab(index=df['SMTS'],columns=df['DTHHRDY'])
display(df)
Tissue_list=df.index

#SVM Method classification for tissues with more than 200 samples 
import random
infiles=listdir(tissue_dir)
def simpleExpressionFilter(counts,min_count):
        """accepts raw counts and a minimum sum count per gene across all samples
        return a boolean array of all genes, which can be applied to any transformed counts.
        True is associated with passing the test.
        """
        keep=np.sum(counts)>min_count
        filtered_counts=counts.loc[:,(keep)] # similar to how the boolean array would be used on any count matrix
        print("Post",filtered_counts.shape[1])
        return keep

accuracy_HARDY_SVM=[]
for TISSUE in Tissue_list:
    cpm=pd.read_csv(tissue_dir/str(TISSUE+"_cpm.tsv"),sep="\t",index_col=0)
    #lcpm=pd.read_csv(tissue_dir/str(TISSUE+"_lcpm.tsv"),sep="\t",index_col=0)
    cdat=pd.read_csv(tissue_dir/str(TISSUE+"_c.tsv"),sep="\t",index_col=0)
    tissue_meta=meta[(meta['SMTS']==TISSUE)]
    cdat=cdat[cdat.index.isin(tissue_meta['SAMPID'])]
    cpm=cpm[cpm.index.isin(tissue_meta['SAMPID'])]
    #X=X.values
    #y=y.values
    pca=PCA(n_components=len(tissue_meta['DTHHRDY']))
    cdat=pca.fit_transform(cdat)
    cdat_train, cdat_test, y_train, y_test = \
        sklearn.model_selection.train_test_split(cdat, tissue_meta['DTHHRDY'], test_size=.3, random_state=1234)
    init=sklearn.svm.classes.SVC(kernel='linear')
    classifier=init.fit(cdat_train, y_train)
    ev=classifier.predict(cdat_test)
    
    ac=0
    y_test=y_test.values
    for i in range (0,len(y_test)):
        if y_test[i]==ev[i]:
            a=0
        else:
            a=1
        ac=ac+a
    ac=1-ac/len(y_test)
    accuracy_HARDY_SVM.append(ac)
    
    
SVM_ACCURACY_HARDY=pd.DataFrame(columns=['tissue', 'accuracy'])
SVM_ACCURACY_HARDY['tissue']=Tissue_list
SVM_ACCURACY_HARDY['accuracy']=accuracy_HARDY_SVM
SVM_ACCURACY_HARDY.to_csv('SVM_ACCURACY_HARDY.csv')

sns.barplot(Tissue_list, accuracy_HARDY_SVM)
plt.title('Accuracy of Different Tissues for SVM Method in Predicting Death Type')
plt.xlabel('Tissue Type')
plt.ylabel('Accuracy')
plt.xticks(rotation=30, ha='right')
plt.savefig('SVM_Accuracy_HARDY.png')


