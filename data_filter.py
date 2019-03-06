#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Benjamin Imlay
import sys
from pathlib import Path
def main(selected_SAMPID_file,counts_file):
	with open(data_dir/selected_SAMPID) as f:
		selected_samples=[SAMPID.rstrip() for SAMPID in f]
	with open(data_dir/counts_file) as f:
		for l in f:
			SAMPID=l.split("\t")[0]
			print(SAMPID)
if __name__ == '__main__':
	data_dir=Path("data")
	selected_SAMPID="filteredSAMPID.tsv"
	main(selected_SAMPID,sys.argv[1])