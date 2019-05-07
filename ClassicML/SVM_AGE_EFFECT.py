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


sns.set_style("darkgrid")
data_dir=Path("data")
tissue_dir=Path("tissue-specific")

manifest={
    "data":"All_Tissue_Site_Details.combined.reads.gct",
    "sample_meta":"GTEx_v7_Annotations_SampleAttributesDS.txt",
    "subject_meta":"GTEx_v7_Annotations_SubjectPhenotypesDS.txt",
    "merged_meta":"merged_meta.tsv"}

meta=pd.read_csv(manifest['merged_meta'],sep="\t",dtype={'SMUBRID':object,'SEX':object,'DTHHRDY':object})
#meta

meta=meta[~(meta['AGE'].isnull())] # removes all samples without age

counts=pd.DataFrame(meta['SMTS'].value_counts())
#display(counts)


df=meta[meta['SMTS'].isin(counts[counts['SMTS']>200].index)]
df=pd.crosstab(index=df['SMTS'],columns=df['AGE'])
#display(df)
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
        #print("Post",filtered_counts.shape[1])
        return keep

accuracy=[]
for TISSUE in Tissue_list:
    cpm=pd.read_csv(tissue_dir/str(TISSUE+"_cpm.tsv"),sep="\t",index_col=0)
    #lcpm=pd.read_csv(tissue_dir/str(TISSUE+"_lcpm.tsv"),sep="\t",index_col=0)
    cdat=pd.read_csv(tissue_dir/str(TISSUE+"_c.tsv"),sep="\t",index_col=0)
    tissue_meta=meta[(meta['SMTS']==TISSUE)]
    cdat_train, cdat_test, y_train, y_test = \
        sklearn.model_selection.train_test_split(cdat, tissue_meta['AGE'], test_size=.3, random_state=1234)
    cpm_train, cpm_test, y_train, y_test = \
        sklearn.model_selection.train_test_split(cpm, tissue_meta['AGE'], test_size=.3, random_state=1234) # random state guarantees that the same split is made for a given tissue.
    sum(cpm_train.iloc[:,0]) # Confirms that the split is the same each time
    # Biased against samples with a smaller library size
    keep_expr=simpleExpressionFilter(cdat_train,10)
    cpm_train_expression_filter=cpm_train.loc[:,(keep_expr)]
    cpm_test_expression_filter=cpm_test.loc[:,(keep_expr)]
    selector=sklearn.feature_selection.VarianceThreshold(threshold=.1)
    selector.fit(cpm_train_expression_filter)
    var_keep=selector.get_support(indices=True)
    train_final=cpm_train_expression_filter.iloc[:,var_keep]
    test_final=cpm_test_expression_filter.iloc[:,var_keep]
    X=train_final.append(test_final)
    y=y_train.append(y_test)
    X=X.values
    y=y.values
    pca=PCA(n_components=len(y))
    X=pca.fit_transform(X)
    #X.shape
    total_number=len(y)
    train_number=round(0.8*total_number)
    random_vals=random.sample(range(total_number+1), train_number)
    X_train=X[[i for i in range(0,total_number) if i in random_vals]]
    X_test=X[[i for i in range(0,total_number) if i not in random_vals]]
    y_train=y[[i for i in range(0,total_number) if i in random_vals]]
    y_test=y[[i for i in range(0,total_number) if i not in random_vals]]
    #['linear', 'rbf', 'poly', 'sigmoid']
    init=sklearn.svm.classes.SVC(kernel='linear')
    classifier=init.fit(X_train, y_train)
    ev=classifier.predict(X_test)
    ac=0
    for i in range (0,len(y_test)):
        if y_test[i]==ev[i]:
            a=0
        else:
            a=1
        ac=ac+a
    ac=1-ac/len(y_test)
    accuracy.append(ac)    




SVM_ACCURACY_FILE=pd.DataFrame(columns=['tissue', 'accuracy'])
SVM_ACCURACY_FILE['tissue']=Tissue_list
SVM_ACCURACY_FILE['tissue']=Tissue_list
SVM_ACCURACY_FILE.to_csv('SVM_ACCURACY_FILE.csv')





sns.barplot(Tissue_list, accuracy)
plt.title('Accuracy of Different Tissues for SVM method')
plt.xlabel('Tissue Type')
plt.ylabel('Accuracy')
plt.xticks(rotation=30, ha='right')
plt.savefig('SVM_Accuracy.png')

plt.show()


### Exporting accuracy data as csv file
KNN_ACCURACY_FILE=pd.DataFrame(columns=['tissue', 'accuracy'])
KNN_ACCURACY_FILE['tissue']=Tissue_list
KNN_ACCURACY_FILE['accuracy']=accuracy_knn
KNN_ACCURACY_FILE.to_csv('KNN_ACCURACY_FILE.csv')

