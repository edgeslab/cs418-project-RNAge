import Init as init

import os
from pathlib import Path
import pandas as pd
import numpy as np
import sklearn

from os import listdir

import warnings
warnings.filterwarnings("ignore")

meta = init.meta

meta=meta[~(meta['AGE'].isnull())] # removes all samples without age

counts=pd.DataFrame(meta['SMTS'].value_counts())

data_dir=Path("data")
tissue_dir=Path("tissue-specific")

#There are many tissues with >200 samples with age recorded. Only tissues with 200 samples or more will be considered for predictive analysis.

df=meta[meta['SMTS'].isin(counts[counts['SMTS']>200].index)]
df=pd.crosstab(index=df['SMTS'],columns=df['AGE'])

#Filter data by expression level with the following function
def filter_by_expr(counts,min_count=None,min_sample=None,grp=None):
    lib_size=np.sum(counts,axis=1)
    MedianLibSize=np.median(lib_size)
    norm_cutoff=min_count/MedianLibSize*1e6
    print(norm_cutoff)
    gene_counts=np.sum(counts)
    
#Filter by row count threshold

TISSUE='Colon'
infiles=listdir(data_dir/tissue_dir)
TISSUE_files=[f for f in infiles if  TISSUE in f]

cpm=pd.read_csv(data_dir/tissue_dir/str(TISSUE+"_cpm.tsv"),sep="\t",index_col=0)
cdat=pd.read_csv(data_dir/tissue_dir/str(TISSUE+"_c.tsv"),sep="\t",index_col=0)

tissue_meta=meta[(meta['SMTS']==TISSUE)]

cdat_train, cdat_test, y_train, y_test = \
        sklearn.model_selection.train_test_split(cdat, tissue_meta['AGE'], test_size=.3, random_state=1234)
cpm_train, cpm_test, y_train, y_test = \
        sklearn.model_selection.train_test_split(cpm, tissue_meta['AGE'], test_size=.3, random_state=1234)

# Biased against samples with a smaller library size
def simpleExpressionFilter(counts,min_count):
    """accepts raw counts and a minimum sum count per gene across all samples
    return a boolean array of all genes, which can be applied to any transformed counts.
    True is associated with passing the test.
    """
    keep=np.sum(counts)>min_count
    filtered_counts=counts.loc[:,(keep)] # similar to how the boolean array would be used on any count matrix
    return(keep)
print("Original Gene Count: "+str(len(cpm_train.columns)))
keep_expr=simpleExpressionFilter(cdat_train,10)
cpm_train_expression_filter=cpm_train.loc[:,(keep_expr)]

cpm_test_expression_filter=cpm_test.loc[:,(keep_expr)]
print("Expression Filter Gene Count: "+str(len(cpm_train_expression_filter.columns)))
#Filter by low variance variance and plot mean distribution of genes after each step explained above.

selector=sklearn.feature_selection.VarianceThreshold(threshold=.1)
selector.fit(cpm_train_expression_filter)
var_keep=selector.get_support(indices=True)
train_final=cpm_train_expression_filter.iloc[:,var_keep]
test_final=cpm_test_expression_filter.iloc[:,var_keep]
print("+Variance Filter Gene Count: "+str(len(train_final.columns)))
display(train_final.head(3))