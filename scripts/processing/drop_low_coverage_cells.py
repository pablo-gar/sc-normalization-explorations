import scanpy as sc
import sys
import os
import numpy as np


min_expressed_genes = int(sys.argv[1])
input_file = sys.argv[2]
output_file = sys.argv[3]

adata = sc.read(input_file)

# Make raw the main layer
adata.X = adata.raw.X.copy()
del adata.raw

# Get non_zero values per cell
non_zero = [adata.X[i,:].count_nonzero() for i in range(adata.n_obs)]

# Remove those below threshold 
adata = adata[[i >= min_expressed_genes for i in non_zero],:]

adata.write_h5ad(output_file, compression="gzip")
