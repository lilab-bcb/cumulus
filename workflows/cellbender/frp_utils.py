import pandas as pd
import pegasusio as io

from pegasusio import MultimodalData, UnimodalData
from typing import Union, Optional


def filter_included_genes(
    data: Union[MultimodalData, UnimodalData],
    filt_features_tsv: pd.DataFrame,
    output_file: Optional[str] = None,
) -> None:
    df_var = pd.read_csv(filt_features_tsv, sep='\t', header=None)
    assert df_var[0].nunique() == 18082, "Filtered genes inconsistent!"

    data._inplace_subset_var(data.var['featureid'].isin(df_var[0].values))
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
