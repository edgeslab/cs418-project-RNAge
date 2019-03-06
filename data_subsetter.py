#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Benjamin Imlay
import pandas as pd
from pathlib import Path
def GTEx_sample_shrinker(meta,by_col,n=20):
    by=meta[by_col].unique()
    ans=[]
    for i in by:
        nTissue=len(meta[meta[by_col]==i])
        if nTissue<n:
            nn=nTissue
        else:
            nn=n
        ans.append(meta[meta[by_col]==i].sample(nn)['SAMPID'])
    selectedMeta=pd.concat(ans)
    return selectedMeta
if __name__ == '__main__':
    data_dir=Path("data")
    manifest={"data":"All_Tissue_Site_Details.combined.reads.gct",
              "sample_meta":"GTEx_v7_Annotations_SampleAttributesDS.txt",
              "subject_meta":"GTEx_v7_Annotations_SubjectPhenotypesDS.txt"}
    meta=pd.read_csv(data_dir/manifest['sample_meta'],sep="\t")
    y=GTEx_sample_shrinker(meta,'SMTS',20)
    y.to_csv(data_dir/"filteredSAMPID.tsv",sep="\t",index=False)
