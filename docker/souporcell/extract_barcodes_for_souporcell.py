#!/usr/bin/env python

from sys import argv, exit
import pegasusio

if len(argv) != 4:
	print("Usage: python extract_barcodes_for_souporcell.py input_raw.h5 output_barcodes.tsv ngene")
	exit(-1)

data = pegasusio.read_input(argv[1], ngene = int(argv[3]))

with open(argv[2], "w") as fout:
	fout.write('\n'.join([x + '-1' for x in data.obs_names]) + '\n')
