import anndata as ad
import numpy as np
import scipy as sc
import sys
from sklearn.preprocessing import QuantileTransformer

input_file = sys.argv[1]
output_file = sys.argv[2]

adata = ad.read(input_file)

# Get a copy of raw data 
adata.layers["rankit"] = adata.X.copy()
    
# Perform normalizations
for i in range(adata.n_obs):
    
    non_zero = adata.layers["rankit"][i,:].nonzero()[1]
    x = adata.layers["rankit"][i, non_zero].toarray().transpose()
    
    qt = QuantileTransformer(n_quantiles=len(x), random_state=0, output_distribution='normal')
    x_norm = qt.fit_transform(x, 1000)
    adata.layers["rankit"][i, non_zero] = x_norm.transpose()
    
    
adata.write_h5ad(output_file, compression="gzip")


