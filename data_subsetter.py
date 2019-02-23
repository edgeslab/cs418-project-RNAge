#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Benjamin Imlay
import pandas as pd
def GTEx_shrinker(counts,meta,by_col,n=20):
    by=meta[by_col].unique()
    ans=[]
    for i in by:
        nTissue=len(meta[meta[by_col]==i])
        if nTissue<n:
            nn=nTissue
        else:
            nn=n
        ans.append(meta[meta[by_col]==i].sample(nn))
    print(ans)
    selectedMeta=pd.concat(ans)
    selectedCounts=counts.iloc[selectedMeta.index,:]
    return selectedCounts,selectedMeta
if __name__ == '__main__':
    #filteredData=copy.deepcopy(data)
    x,y=GTEx_shrinker(data.rawCounts,data.sampleMeta,'SMTS',20)
    x.to_csv("filteredData.tsv",sep="\t")
    y.to_csv("filteredMeta.tsv",sep="\t")