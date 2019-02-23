#!/usr/bin/env python
# encoding: utf-8
# Benjamin Imlay
import pandas as pd
from pathlib import Path
class expressionData():
    def __init__(self,data_dir,manifest):
        """
        The constructor expects the data dir of type Path from pathlib, a manifest dictionary with data, sample_meta, and subject_meta.
        
        """
        self.rawCounts=self.readRawCounts(data_dir/manifest["data"])
        self.subjMeta=self.readMeta(data_dir/manifest["subject_meta"])
        self.sampleMeta=self.readMeta(data_dir/manifest["sample_meta"])
        self.sampleMeta=self.orderMeta()
        self.sampleMeta=self.mergeSubjSamp()
    def orderMeta(self):
        """"
        This function reads, removes, and orders metadata entries to align with the rawcounts matrix
        Input:
            self
        Output:
            DataFrame of ordered sample metadata
        """
        SAMPID_DataFrame=pd.DataFrame({'SAMPID':pd.Series(self.rawCounts.index)})
        meta=SAMPID_DataFrame.merge(self.sampleMeta,how='left',on='SAMPID').set_index('SAMPID')
        print("Trimmed and Ordered metadata to dimensions: ",meta.shape)
        return meta
    def mergeSubjSamp(self):
        """
        This function parses the index of samples to get the subject prefix. Then the subject metadata is merged with the sample metadataw.
        Input:
            self.sampleMeta
            self.subjMeta
        Output:
            Dataframe of ordered sample metadata.
        """
        meta=self.subjMeta
        SAMPmeta=self.sampleMeta
        SAMPmeta['SUBJID']=pd.Series(SAMPmeta.index,index=SAMPmeta.index).str.rsplit('-',n=3).str.get(0)
        SAMPmeta=SAMPmeta.merge(meta,on="SUBJID",how="left")
        return SAMPmeta
    def readMeta(self,data_path):
        """
        This function reads raw GTEx metadata file in a tsv format.
        Input:
            Path (str) to the sample metadata file.
        Output:
            Pandas dataframe of sample metadata.
        """
        data=pd.read_csv(data_path,sep="\t")
        print("Imported metadata of dimensions",data.shape)
        return data
    def readRawCounts(self,data_path):
        """
        This function reads raw counts in the gct format, which has two lines of metadata, and two columns for gene name.
        Input:
            Path (str) to the .gct file.
        Output:
            Pandas dataframe with integer count matrix and multiindex.
        """
        data=pd.read_csv(data_path,sep="\t",skiprows=2,index_col=[0,1])
        print("Imported raw count data of dimensions",data.T.shape)
        return data.T
def main():
    pass
if __name__ == '__main__':
    data_dir=Path("C:\data")
    
    
    manifest={"data":"All_Tissue_Site_Details.combined.reads.gct",
              "sample_meta":"GTEx_v7_Annotations_SampleAttributesDS.txt",
              "subject_meta":"GTEx_v7_Annotations_SubjectPhenotypesDS.txt"}
    data=expressionData(data_dir,manifest)