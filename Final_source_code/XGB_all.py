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

def predOneTissue(alg, dtrain,y_train, useTrainCV=True, cvFolds=5, early_stopping_rounds=50):
    if useTrainCV:
        xgb_param = alg.get_xgb_params()
        extra = {'num_class': len(y_train.unique())}
        xgbParams = alg.get_xgb_params()
        xgbParams.update(extra)
        xgTrain = xgb.DMatrix(dtrain, label=y_train)
        cvresult = xgb.cv(xgbParams,
                      xgTrain,
                      num_boost_round=alg.get_params()['n_estimators'],
                      nfold=cvFolds,
                      stratified=True,
                      metrics={'mlogloss'},
                      early_stopping_rounds=early_stopping_rounds,
                      seed=0,
                      callbacks=[xgb.callback.print_evaluation(show_stdv=False),
                                 xgb.callback.early_stop(3)])

        #print(cvresult)
        alg.set_params(n_estimators=cvresult.shape[0])
    alg.set_params(num_class=len(y_train.unique()))
    #Fit alg
    print(dtrain.shape,y_train.shape)
    alg.fit(dtrain,y_train,eval_metric='auc')
    return alg
def scoreOneTissue(y_test,y_preds):
    """
    Returns a tuple of raw score and F1.
    """
    return (metrics.accuracy_score(y_test,y_preds))
def splitData(X,y,rng=1234):
    return sklearn.model_selection.train_test_split(X, y, test_size=.3, random_state=rng)
def simpleExpressionFilter(counts,min_count):
    """accepts raw counts and a minimum sum count per gene across all samples
    return a boolean array of all genes, which can be applied to any transformed counts.
    True is associated with passing the test.
    """
    keep=np.sum(counts)>min_count
    #print("Pre",counts.shape[1])
    #filtered_counts=counts.loc[:,(keep)] # similar to how the boolean array would be used on any count matrix
    #print("Post",filtered_counts.shape[1])
    return(keep)
def loadData(TISSUE,data_dir,tissue_dir):
    cpm=pd.read_csv(data_dir/tissue_dir/str(TISSUE+"_cpm.tsv"),sep="\t",index_col=0)
    #lcpm=pd.read_csv(data_dir/tissue_dir/str(TISSUE+"_lcpm.tsv"),sep="\t",index_col=0)
    cdat=pd.read_csv(data_dir/tissue_dir/str(TISSUE+"_c.tsv"),sep="\t",index_col=0)
    return (cpm,cdat)
def mainLoop(df,data_dir,tissue_dir):
    data_dir=Path("data")
    tissue_dir=Path("tissue-specific")
    results=[]
    algs=[]
    TISSUE_list=df['SMTS'].unique()
    for TISSUE in TISSUE_list:
        cpm,cdat=loadData(TISSUE,data_dir,tissue_dir)
        #print(cpm.shape)
        #print(cdat.shape)
        all_age=df.loc[(df['SMTS']==TISSUE),'AGE']
        #print(all_age.shape)
        cpm_train,cpm_test,y_train,y_test=splitData(cpm,all_age)
        c_train,c_test,y_train,y_test=splitData(cdat,all_age)
        #y_train.map({'20-29':0,'30-39':1,'40-49':2,'50-59':3, '60-69':4,'70-79':5})
        #y_test.map({'20-29':0,'30-39':1,'40-49':2,'50-59':3, '60-69':4,'70-79':5})
        y_train.replace({"20-29":0, "30-39":1, "40-49":2, "50-59":3, "60-69":4, "70-79":5}, inplace = True)
        y_test.replace({"20-29":0, "30-39":1, "40-49":2, "50-59":3, "60-69":4, "70-79":5}, inplace = True)
        #print(y_train.shape)
        #print(cpm_train.shape)
        #print(c_train.shape)
        keep=simpleExpressionFilter(c_train,10)
        cpm_train=cpm_train.loc[:,(keep)]
        cpm_test=cpm_test.loc[:,(keep)]
        
        selector=sklearn.feature_selection.VarianceThreshold(threshold=.1)
        selector.fit(cpm_train)
        keep=selector.get_support(indices=True)
        cpm_train=cpm_train.iloc[:,keep]
        cpm_test=cpm_test.iloc[:,keep]
        
        xgb1 = XGBClassifier(
         learning_rate =0.1,
         n_estimators=1000,
         max_depth=5,
         min_child_weight=1,
         gamma=0,
         subsample=0.8,
         colsample_bytree=0.8,
         objective= 'multi:softprob',
         nthread=30,
         scale_pos_weight=1,
         seed=1234)
        xgb1=predOneTissue(xgb1,cpm_train,y_train)
        y_preds=xgb1.predict(cpm_test)
        score=scoreOneTissue(y_test,y_preds)
        #print(y_test,y_preds)
        results.append(score)
        algs.append(xgb1)
    results=pd.Series(results,index=TISSUE_list)
    return (algs,results)
def buildModels():
    algs,results=mainLoop(df,data_dir,tissue_dir)
    results.to_csv("models/results.tsv",sep="\t",header=['Score'])
    [joblib.dump(algs[x],'models/{}_model.pkl'.format(results.index.values[x])) for x in range(len(results))]
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
        x=feat_imp['Score'][:10].plot(kind='bar', title='Feature Importances for {}'.format(TISSUE))
        #x.set_ylabel('Feature Importance Score')
        plt.show(x)
        plt.savefig(fname='plots/{}_feat_imp.png'.format(TISSUE))
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