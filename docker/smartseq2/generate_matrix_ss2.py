#!/usr/bin/env python

import os
from sys import argv, exit
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import pegasusio as pio


if len(argv) != 3:
    print('Usage: ./generate_matrix_ss2.py gene_results output_name')
    exit(-1)


gene_results = argv[1]
output_name = argv[2]

barcodes = []
plates = []
total_reads = []
alignment_rates = []
unique_rates = []

gene_ids = []
gene_names = []

X = []

for result_file in gene_results.split(','):
    dirname, basename = os.path.split(result_file)
    slen = len('.genes.results')
    basename = basename[:-slen]
    barcode, sep, plate = basename.rpartition('.')

    barcodes.append(barcode)
    plates.append(plate)

    df = pd.read_table(result_file, header = 0, index_col = 0)

    if len(gene_ids) == 0:
        for x in df.index:
            items = x.split('_')
            gene_ids.append(items[0])
            gene_names.append('_'.join(items[1:]))

    tot_counts = df['expected_count'].sum()
    counts = df['TPM'].values
    denom = counts.sum()
    if denom > 0:
        counts = (counts / denom * tot_counts + 0.5)
    counts = counts.astype(int)
    X.append(counts)

    cnt_file = f'{dirname}/{basename}.stat/{basename}.cnt'
    with open(cnt_file) as fin:
        Ns = [int(x) for x in next(fin).strip().split(' ')]
        align_values = [int(x) for x in next(fin).strip().split(' ')]
        total_reads.append(Ns[3])
        alignment_rates.append(round(Ns[1] * 100.0 / Ns[3], 2))
        unique_rates.append(round(align_values[0] * 100.0 / Ns[3], 2))

X = csr_matrix(np.array(X, dtype = np.int32))
unidata = pio.UnimodalData({'barcodekey': barcodes, 'plate': plates, 'total_reads': total_reads, 'alignment_rate': alignment_rates, 'unique_rate': unique_rates},
    {'featurekey': gene_names, 'featureid': gene_ids},
    {'X': X},
    {'modality': 'rna', 'genome': 'unknown'}
)
pio.write_output(unidata, output_name, file_type = 'mtx')
os.system(f'mv {output_name}/unknown-rna/*.* {output_name}/')
os.system(f'rm -rf {output_name}/unknown-rna')
