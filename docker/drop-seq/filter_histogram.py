#!/usr/bin/env python

import argparse
import os

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--histogram')
parser.add_argument('--whitelist')
parser.add_argument('--output')

args = parser.parse_args()
# #INPUT=/cromwell_root/test_aligned_tagged_repaired.bam	TAG=XC	FILTER_PCR_DUPLICATES=false	READ_QUALITY=10
# 3155	TGGCCGAACATC
if args.whitelist != '':
    whitelist = pd.read_csv(args.whitelist, sep='\n', header=None, index_col=0)
    histogram = pd.read_csv(args.histogram, sep='\t', header=None, comment='#')
    histogram = histogram[histogram[1].isin(whitelist.index)]
    os.rename(args.histogram, args.histogram.replace('_tag.txt', '_unfiltered_tag.txt'))
    histogram.to_csv(args.output, sep='\t', index=False, header=False)
