import sys
import pandas as pd
import numpy as np
import pegasusio as io

def generate_CNV_scored_zarr(result_zarr, cnv_metrics_csv, sample_anno_csv, out_name):

    data = io.read_input(result_zarr)
    barcodes = data.obs_names.values

    df_scores = pd.read_csv(cnv_metrics_csv, index_col = 0)
    cnv_barcodes = df_scores.index.values

    if barcodes.size != cnv_barcodes.size or np.sum(barcodes != cnv_barcodes) > 0:
        print("Cell barcodes are not consistent with CNV result!")
        sys.exit(1)

    df_scores['corr_cnv'] = df_scores['corr_cnv'].fillna(1.0)

    data.obs['cnv_scores'] = df_scores['cnv_scores'].values
    data.obs['cnv_corr'] = df_scores['corr_cnv'].values

    cnv_prod = df_scores['cnv_scores'].values * (df_scores['corr_cnv'].values ** 2)
    data.obs['sqrt_cnv'] = np.sqrt(cnv_prod)

    df_anno = pd.read_csv(sample_anno_csv, header = None, sep = '\t')
    anno_barcodes = df_anno[0].values

    if barcodes.size != anno_barcodes.size or np.sum(barcodes != anno_barcodes) > 0:
        print("Cell barcodes are not consistent with sample annotation result!")
        sys.exit(1)

    data.obs['infercnv_cell_types'] = df_anno[1].astype('str').values

    io.write_output(data, f"{out_name}.cnv.zarr.zip")

if __name__ == '__main__':
    result_zarr = sys.argv[1]
    cnv_metrics_csv = sys.argv[2]
    sample_anno_csv = sys.argv[3]
    out_name = sys.argv[4]

    generate_CNV_scored_zarr(result_zarr, cnv_metrics_csv, sample_anno_csv, out_name)
