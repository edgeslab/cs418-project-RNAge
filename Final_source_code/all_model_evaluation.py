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
sns.set_style("whitegrid")
if __name__=="__main__": 
    results_dir=Path("model_accuracies")
    infiles=listdir(results_dir)
    df=pd.read_csv(results_dir/infiles[0],sep='\t')
    for f in infiles[1:]:
        x=pd.read_csv(results_dir/f,sep='\t')
        name=x.columns[0]
        df[name]=x[name]
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    grps = list(set(df.index.values))
    res_df=pd.DataFrame(columns=['statistic','p-value'])
    keys=[]
    for g1 in grps:
        for g2 in grps[grps.index(g1)+1:]:
            if g1 != g2:
                keys.append(str(g1+'_'+g2))
                x = list(df.loc[g1,:])
                y = list(df.loc[g2,:])
                res=mannwhitneyu(x,y,alternative='two-sided')
                res_df=res_df.append({'statistic':res[0],
                                        'p-value':res[1]},ignore_index=True)
    res_df=res_df.set_index(pd.Index(keys,'Tissue'))
    corrected_p_values=multipletests(res_df['p-value'])[1]
    res_df['cor_p-value']=pd.Series(corrected_p_values,index=keys)
    res_df=res_df.sort_values(by='cor_p-value')
    df['Average Accuracy']=df.mean(axis=1)
    df['sdev']=df.std(axis=1)
    plt.figure(figsize=(10, 6))
    ax=sns.barplot(x=df.index.values,y=df['Average Accuracy'],
              yerr=df['sdev']*1, capsize=.2)
    x=ax.set_title("Average Model Accuracies")
    x=ax.set_xlabel("Tissues")
    x=ax.set_ylabel("Average Accuracy - 1 SD")
    x=ax.set_xticklabels(labels=df.index.values,rotation=38)
    fig = ax.get_figure()
    fig.savefig("plots/all_model_accuracy.png",dpi=100,bbox_inches = "tight")
    from IPython.display import display, HTML
    display(res_df.head(10))