#!/usr/bin/env python

from sys import argv, exit
import pegasus as pg

if len(argv) != 4:
	print("Usage: python extract_barcodes_for_souporcell.py input_raw.h5 output_barcodes.tsv ngene")
	exit(-1)

adata = pg.read_input(argv[1])
idx = adata.X.getnnz(axis = 1) >= int(argv[3])

with open(argv[2], "w") as fout:
	fout.write('\n'.join([x + '-1' for x in adata.obs_names[idx]]) + '\n')
