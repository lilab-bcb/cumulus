#!/usr/bin/env python

from sys import argv, exit
import pegasusio as io

if len(argv) != 4:
	print("Usage: python extract_barcodes_from_rna.py input_raw.h5 output_barcodes.tsv ngene")
	exit(-1)

data = io.read_input(argv[1])
data.filter_data(min_genes=int(argv[3]))

with open(argv[2], "w") as fout:
	fout.write('\n'.join([x + '-1' for x in data.obs_names]) + '\n')
