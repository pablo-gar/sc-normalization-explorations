import anndata as ad
import numpy as np
import sys
import os
from scipy import io
from scipy import sparse
import gzip

layers = sys.argv[1].split(",")
input_file = sys.argv[2]
prefix = sys.argv[3]

# Reading
adata = ad.read(input_file)

# writing obs and var
adata.obs.to_csv(prefix + "_obs.gz", sep="\t", index = False, header = True)
adata.var.to_csv(prefix + "_var.gz", sep="\t", index = False, header = True)

for layer in layers:
    
    
    if layer == "X":
        x = adata.X
    else:
        if layer in adata.layers:
            x = adata.layers[layer]
        else:
            continue
        
    if sparse.issparse(x):
        out_file = prefix + "_" + layer + "_mtx.gz"
        with gzip.open(out_file, "wb") as out:
            io.mmwrite(out, x)
    else:
        out_file = prefix + "_" + layer + "_txt.gz"
        np.savetxt(out_file, x, fmt="%s")
