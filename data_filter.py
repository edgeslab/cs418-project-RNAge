#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Benjamin Imlay
import sys
from pathlib import Path
## sys.argv[1] should be the count file (sample x  gene) you want to filter.
## filteredSAMPID.tsv is simply a list of SAMPID's to include. In this case, it is from data_subsetter.py
def main(selected_SAMPID_file,counts_file):
	with open(selected_SAMPID) as f:
		selected_samples=[SAMPID.rstrip() for SAMPID in f]
	with open(counts_file) as f:
		print(f.readline().rstrip('\n'))
		for l in f:
			SAMPID=l.split("\t")[0]
			if SAMPID in selected_samples:
				print(l.rstrip('\n'))
if __name__ == '__main__':
	data_dir=Path("data")
	selected_SAMPID="data/filteredSAMPID.tsv"
	main(selected_SAMPID,sys.argv[1])
