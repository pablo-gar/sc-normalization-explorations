import anndata as ad
import numpy as np
from scipy import sparse
import scanpy
import sys
import _utils
from sklearn.utils import sparsefuncs

mode = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

adata = ad.read(input_file)

# Get a copy of raw data 
temp_X = adata.X.copy()
    
# Perform normalizations
scanpy.pp.normalize_total(adata, key_added = "normalization_factor_sum")

# Scale and center
if mode == "rms":
    mean, rms = _utils._get_mean_var(adata.X, RSM=True)
    layer_name = "inmf"
elif mode == "sd":
    mean, rms = _utils._get_mean_var(adata.X, RSM=False)
    layer_name = "z_score"
else:
    raise ValueError("Only modes available are 'rms' and 'sd'")

rms[rms == 0] = 1

if sparse.issparse(adata.X):
    adata.X = (adata.X.toarray() - mean) / rms
    #adata.X = sparse.csr_matrix(adata.X)
else:
    adata.X = (adata.X - mean) / rms 
    
# Store
adata.layers[layer_name] = adata.X.copy()
adata.X = temp_X

# Set zeroes
mask = np.ones(adata.X.shape, dtype=bool)
if sparse.issparse(adata.X):
    nonzeros = adata.X.nonzero()
else:
    nonzeros = np.nonzero(adata.X)
    
mask[nonzeros] = False
adata.layers[layer_name][mask] = 0 
    
if sparse.issparse(adata.X):
    adata.layers[layer_name] = sparse.csr_matrix(adata.layers[layer_name])

adata.write(output_file, compression="gzip")
