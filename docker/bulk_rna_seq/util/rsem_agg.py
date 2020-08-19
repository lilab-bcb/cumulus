#!/usr/bin/env python

import argparse
import json
import os

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Aggregate RSEM results')
parser.add_argument('--gene', help='Path to gene results', nargs='+')
parser.add_argument('--count', help='Path to count results', nargs='+')
parser.add_argument('--qc_vars', help='Path to qc_vars')
parser.add_argument('--normalize_tpm_by_sequencing_depth', action='store_true')
# parser.add_argument('--gzip', help='gzip gene results matrix', action='store_true')

args = parser.parse_args()
count_results = args.count
gene_results = args.gene
qc_vars_json = args.qc_vars
normalize_tpm_by_sequencing_depth = args.normalize_tpm_by_sequencing_depth
gene_names = None
barcodes = []
cntmat = []
prefix = 'tpm.norm' if normalize_tpm_by_sequencing_depth else 'expected_count'
for result_file in gene_results:
    barcodes.append(os.path.basename(result_file)[:-len('.genes.results')])
    df = pd.read_table(result_file, header=0, index_col=0)
    if gene_names is None:
        gene_names = np.array(['_'.join(x.split('_')[1:]) for x in df.index])
    if normalize_tpm_by_sequencing_depth:
        tot_counts = df['expected_count'].sum()
        counts = df['TPM'].values / 10.0  # convert TPMs into TP100Ks
        denom = counts.sum()
        if denom > 0:
            counts = (counts / denom * tot_counts + 0.5).astype(int)
    else:
        counts = df['expected_count'].values
    cntmat.append(counts)
df_out = pd.DataFrame(data=np.stack(cntmat, axis=1), index=pd.Index(gene_names, name='GENE'), columns=barcodes)
df_out.to_csv('{}.dge.txt.gz'.format(prefix), sep='\t', compression='gzip')

qc_df = None
if qc_vars_json is not None:
    with open(qc_vars_json, 'rt') as f:
        qc_vars = json.load(f)

    totals = df_out.values.sum(axis=0)
    idx_upper = df_out.index.str.split('_').str[-1].str.upper().str
    qc_df = pd.DataFrame(index=df_out.columns)
    for name in qc_vars:
        vals = qc_vars[name].split(',')
        selector = None
        for val in vals:
            val = val.upper().strip()
            sel = idx_upper.startswith(val)
            selector = sel if selector is None else selector | sel
        qc_df[name] = df_out[sel].values.sum(axis=0) / totals

count_df = None
if count_results is not None:
    arr = []
    barcodes = []
    for result_file in count_results:
        barcodes.append(os.path.basename(result_file)[:-len('.cnt')])
        with open(result_file) as fin:
            Ns = [int(x) for x in next(fin).strip().split(' ')]
            align_values = [int(x) for x in next(fin).strip().split(' ')]
            res = [str(Ns[3]), str(round(Ns[1] * 100.0 / Ns[3], 2)) + "%",
                   str(round(align_values[0] * 100.0 / Ns[3], 2)) + "%"]
            arr.append(res)
    count_df = pd.DataFrame(data=np.array(arr), index=barcodes,
        columns=["Total reads", "Alignment rate", "Unique rate"])
    count_df.index.name = "id"
if count_df is not None or qc_df is not None:
    if count_df is not None:
        df = count_df
        if count_df is not None:
            if qc_df is not None:
                df = df.join(qc_df)
    else:
        df = qc_df

    df.to_csv('stats.tsv', sep='\t')
