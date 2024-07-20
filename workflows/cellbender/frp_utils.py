import pandas as pd
import pegasusio as io

from pegasusio import MultimodalData, UnimodalData
from typing import Union, Optional


def filter_included_genes(
    data: Union[MultimodalData, UnimodalData],
    output_file: Optional[str] = None,
) -> None:
    df_var = pd.read_csv("Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv", comment="#")
    df_var = df_var.loc[df_var["included"]].copy()
    df_var["gene_name"] = df_var["probe_id"].apply(lambda s: s.split("|")[1])
    assert df_var["gene_id"].nunique() == 18082, "Filtered genes inconsistent!"

    data._inplace_subset_var(data.var['featureid'].isin(df_var["gene_id"].unique()))
    correct_gene_names(data)
    if output_file:
        io.write_output(data, output_file)

def correct_gene_names(data):
    rename_dict = {
        "ENSG00000285053": "GGPS1-TBCE",
        "ENSG00000284770": "TBCE",
        "ENSG00000187522": "HSPA14",
        "ENSG00000284024": "MSANTD7",
        "ENSG00000269226": "TMSB15C",
        "ENSG00000158427": "TMSB15B",
    }
    df_var = data.var.reset_index()
    for gene_id, gene_name in rename_dict.items():
        df_var.loc[df_var['featureid']==gene_id, 'featurekey'] = gene_name
    data.var = df_var.set_index('featurekey')
