import os, sys
import pegasusio as io
import pandas as pd
import numpy as np

counts_file = sys.argv[1]
out_name = sys.argv[2]

data = io.read_input(counts_file)
data.select_matrix('raw.X')

# Get raw count matrix.
df_rc = pd.DataFrame(data=np.transpose(data.X.toarray().astype('int')), index=data.var['featureid'], columns=data.obs_names)
df_rc.to_csv(f"{out_name}.csv", sep='\t')
