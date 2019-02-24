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
        self.sampleMeta=self.readMeta(data_dir/manifest["sample_meta"])
        self.sampleMeta=self.inputMeta()
    def inputMeta(self):
        """"
        This function reads, removes, and orders metadata entries to align with the rawcounts matrix
        Input:
            self
        Output:
            Modified self.sampleMeta
        """
        SAMPID_DataFrame=pd.DataFrame({'SAMPID':pd.Series(self.rawCounts.index)})
        meta=SAMPID_DataFrame.merge(self.sampleMeta,how='left',on='SAMPID').set_index('SAMPID')
        return meta
    def readMeta(self,data_path):
        """
        This function reads raw GTEx metadata file in a tsv format.
        Input:
            Path (str) to the sample metadata file.
        Output:
            Pandas dataframe of sample metadata.
        """
        data=pd.read_csv(data_path,sep="\t")
        #self.sampleMeta=data
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
        data=pd.read_csv(data_path,sep="\t",skiprows=2,index_col=[0,1], encoding="cp1252")
        #self.rawCounts=data.T
        print("Imported raw count data of dimensions",data.T.shape)
        return data.T
def main():
    pass
if __name__ == '__main__':
    data_dir=Path("D:\Github\cs418-project-RNAge\data")
    
    
    manifest={"data":"All_Tissue_Site_Details.combined.reads.gct",
              "sample_meta":"GTEx_v7_Annotations_SampleAttributesDS.txt",
              "subject_meta":"GTEx_v7_Annotations_SubjectPhenotypesDS.txt"}
    data=expressionData(data_dir,manifest)

#Testing branch sanity