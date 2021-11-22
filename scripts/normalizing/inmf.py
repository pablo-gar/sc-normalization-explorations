import anndata as ad
import numpy as np
from scipy import sparse
import scanpy
import sys
import _utils
from sklearn.utils import sparsefuncs

input_file = sys.argv[1]
output_file = sys.argv[2]

adata = ad.read(input_file)

# Get a copy of raw data 
temp_X = adata.X.copy()
    
# Perform normalizations
scanpy.pp.normalize_total(adata, key_added = "normalization_factor_sum")

# Scale and center
mean, rms = _utils._get_mean_var(adata.X, RSM=True)
rms[rms == 0] = 1

if sparse.issparse(adata.X):
    adata.X = (adata.X.toarray() - mean) / rms
    #adata.X = sparse.csr_matrix(adata.X)
else:
    adata.X = (adata.X - mean) / rms 

# Store
adata.layers["inmf"] = adata.X.copy()
adata.X = temp_X

adata.write(output_file, compression="gzip")
