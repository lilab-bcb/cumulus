import sys
import numpy as np
import pandas as pd
import pegasus as pg

zarr_file = sys.argv[1]
out_name = sys.argv[2]

data = pg.read_input(f"{zarr_file}")
data.select_matrix('raw.X')

df_raw_count = pd.DataFrame(data=np.transpose(data.X.toarray()), index=data.var['featureid'], columns=data.obs_names)
df_raw_count.to_csv(f"{out_name}.csv", sep='\t')

df_anno = pd.DataFrame({'cell': data.obs_names, 'anno': data.obs['anno'].values})
df_anno.to_csv(f"{out_name}.sample_annotation.csv", sep='\t', header=False, index=False)
