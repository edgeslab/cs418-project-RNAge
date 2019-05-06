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

import warnings
warnings.filterwarnings("ignore")

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