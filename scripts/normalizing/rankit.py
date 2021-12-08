import anndata as ad
import numpy as np
import scipy as sc
import sys
from sklearn.preprocessing import QuantileTransformer

input_file = sys.argv[1]
output_file = sys.argv[2]

def Q_Q_Prob(data):
    
    ranks = sc.stats.rankdata(data, method="dense")
    
    n = len(data)
    prob_level = []
    
    for i in ranks:
        prob_level.append(np.round((i+1-0.5)/n,5))
    Standard_normal_quantiles = sc.stats.norm.ppf(prob_level)
    
    return Standard_normal_quantiles

adata = ad.read(input_file)

# Get a copy of raw data 
adata.layers["rankit"] = adata.X.copy()
    
# Perform normalizations
for i in range(adata.n_obs):
    
    non_zero = adata.layers["rankit"][i,:].nonzero()[1]
    x = adata.layers["rankit"][i, non_zero].toarray()[0]
    x_norm = Q_Q_Prob(x)
    adata.layers["rankit"][i, non_zero] = x_norm.transpose()
    
    
adata.write_h5ad(output_file, compression="gzip")


