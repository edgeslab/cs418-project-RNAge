#Anna
import pandas as pd
import numpy as np
import cmapPy.pandasGEXpress.parse_gct as pg

#the cmapPy for GCT files

#parsing 
df=pg.parse('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct')
