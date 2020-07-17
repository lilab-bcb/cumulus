import argparse

import pegasusio as io
import pandas as pd

parser = argparse.ArgumentParser(description='Merge demuxlet result with gene-count matrix.')
parser.add_argument('demux_res', metavar = 'demux_result.best', help = 'Demuxlet demultiplexing results.')
parser.add_argument('raw_mat', metavar = 'raw_feature_bc_matrix.h5', help = 'Raw gene count matrix in 10x format.')
parser.add_argument('out_file', metavar = 'output_result.zarr.zip', help = 'Output zarr file.')
args = parser.parse_args()

demux_type_dict = {'SNG': 'singlet', 'DBL': 'doublet', 'AMB': 'unknown'}

def write_output(assignment_file: str, input_mat_file: str, output_zarr_file: str) -> None:
    df = pd.read_csv(assignment_file, sep = '\t', header = 0, index_col = 'BARCODE')
    df.index = pd.Index([x[:-2] for x in df.index])
    df['demux_type'] = df['DROPLET.TYPE'].apply(lambda s: demux_type_dict[s])
    df['assignment'] = ''
    df.loc[df['demux_type'] == 'singlet', 'assignment'] = df.loc[df['demux_type'] == 'singlet', 'SNG.BEST.GUESS']
    df.loc[df['demux_type'] == 'doublet', 'assignment'] = df.loc[df['demux_type'] == 'doublet', 'DBL.BEST.GUESS'].apply(lambda s: ','.join(s.split(',')[:-1]))

    data = io.read_input(input_mat_file)
    data.obs['demux_type'] = ''
    data.obs['assignment'] = ''

    idx = data.obs_names.isin(df.index)
    barcodes = data.obs_names[idx]
    df_valid = df.loc[barcodes, ['demux_type', 'assignment']]
    data.obs.loc[idx, 'demux_type'] = df_valid['demux_type'].values
    data.obs.loc[idx, 'assignment'] = df_valid['assignment'].values

    io.write_output(data, output_zarr_file)


if __name__ == '__main__':
    write_output(args.demux_res, args.raw_mat, args.out_file)
