#!/usr/bin/env python

import os
from sys import argv, exit
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import pegasusio as pio


if len(argv) != 4:
    print('Usage: ./generate_matrix_ss2.py gene_results cnt_results output_name')
    exit(-1)


def _map_basename_to_file(results, suffix):
    basenames = []
    basename2file = {}
    for result_file in results.split(','):
        basename = os.path.basename(result_file)[0:-len(suffix)]
        basenames.append(basename)
        basename2file[basename] = result_file
    return basename2file, basenames


b2g, basenames = _map_basename_to_file(argv[1], '.genes.results')
b2c, _ = _map_basename_to_file(argv[2], '.cnt')
output_name = argv[3]

barcodes = []
plates = []
total_reads = []
alignment_rates = []
unique_rates = []

gene_ids = []
gene_names = []

X = []

for basename in basenames:
    barcode, sep, plate = basename.rpartition('.')

    barcodes.append(barcode)
    plates.append(plate)

    assert basename in b2g
    df = pd.read_table(b2g[basename], header = 0, index_col = 0)

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

    assert basename in b2c
    cnt_file = b2c[basename]
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
