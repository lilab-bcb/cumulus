import sys
import numpy as np
import pandas as pd
import pegasusio as io

zarr_file = sys.argv[1]
out_name = sys.argv[2]

data = io.read_input(f"{zarr_file}")

mat_key = None
if 'counts' in data._unidata.matrices:
    mat_key = 'counts'
elif 'raw.X' in data._unidata.matrices:
    mat_key = 'raw.X'
else:
    raise Exception("Cannot find raw counts! Must be specified with key 'counts' or 'raw.X'!")
data.select_matrix(mat_key)

df_raw_count = pd.DataFrame(data=np.transpose(data.X.toarray()), index=data.var['featureid'], columns=data.obs_names)
df_raw_count.to_csv(f"{out_name}.csv", sep='\t')

df_anno = pd.DataFrame({'cell': data.obs_names, 'anno': data.obs['anno'].values})
df_anno.to_csv(f"{out_name}.sample_annotation.csv", sep='\t', header=False, index=False)
