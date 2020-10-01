import sys
import pandas as pd
import pegasus as pg

zarr_file = sys.argv[1]
out_name = sys.argv[2]

data = pg.read_input(f"{zarr_file}")

df_anno = pd.DataFrame({'cell': data.obs_names, 'anno': data.obs['anno'].values})
df_anno.to_csv(f"{out_name}.sample_annotation.csv", sep='\t', header=False, index=False)
