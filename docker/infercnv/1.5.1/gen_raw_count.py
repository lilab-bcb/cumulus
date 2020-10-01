import os, sys
import pegasusio as io
import pandas as pd
import numpy as np

all_counts_file = sys.argv[1]
used_counts_file = sys.argv[2]
out_name = sys.argv[3]

adata = io.read_input(all_counts_file)
bdata = io.read_input(used_counts_file)

# Use only cells and genes considered in analysis.
adata = adata[bdata.obs_names, bdata.var_names].copy()

# Get raw count matrix.
df_rc = pd.DataFrame(data=np.transpose(adata.X.toarray().astype('int')), index=adata.var['featureid'], columns=adata.obs_names)
df_rc.to_csv(f"{out_name}.csv", sep='\t')
