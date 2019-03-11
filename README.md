# RNAge

Data should be downloaded and placed into directory data/

Subsequent count files, meta subsets, and count subsets also go into data/.
## Data subsetting
1. Generate a filteredSAMPID.tsv and filteredMeta.tsv files.
`$ python data_subsetter.py`
2. Read any count file *(sample x gene)* and pipe out rows with SAMPID's that match those contained in filteredSAMPID.tsv.
`$ python data_filter.py data/count_file.tsv > data/filtered_counts.tsv`

