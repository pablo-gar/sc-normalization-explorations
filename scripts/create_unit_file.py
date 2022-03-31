import numpy as np
import scipy as sc
import anndata as ad
import sys
import json
from sklearn.preprocessing import QuantileTransformer

# Creates a test file that has undergone all 
# the processing steps for where's my gene

def main():
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    input_file = "./tabula_muris_lung.h5ad"
    input_file = "./lung_map.h5ad"
    
    genes = ["MALAT1", "CCL5"]
    min_expressed_genes = 500
    min_reads = 3
    offset = 3
    
    adata = ad.read(input_file)
    
    # Get a copy of raw data 
    adata.layers["rankit"] = adata.raw.X.copy()
    
    adata = remove_low_coverage_cells(adata, min_expressed_genes)
    adata = rankit_normalize(adata, offset)
    adata = set_low_to_na(adata, min_reads)
    
    # Get everything in json file
    output = dict()
    output["n_cells"] = adata.n_obs
    output["n_genes"] = adata.n_var
    output["rankit"] = dict()
    for gene in genes:
        output["rankit"][gene] = dict()
        gene_expression = adata.layers["rankit"][:,adata.var_names == gene]
        non_zero_cells = gene_expression.nonzero()[1]
        
    

def remove_low_coverage_cells(adata, n_genes):

    # Get non_zero values per cell
    non_zero = [adata.raw.X[i,:].count_nonzero() for i in range(adata.n_obs)]
    # Remove those below threshold 
    adata = adata[[i >= n_genes for i in non_zero],:]
    
    return adata

def set_low_to_na(adata, n_counts):
    
    rows,cols = adata.raw.X.nonzero()
    for row,col in zip(rows,cols):
        if adata.raw.X[row,col] <= n_counts:
            adata.layers["rankit"][row,col] = 0
    
    return adata

def rankit_normalize(adata, offset):
    
    for i in range(adata.n_obs):
        
        non_zero = adata.layers["rankit"][i,:].nonzero()[1]
        x = adata.layers["rankit"][i, non_zero].toarray()[0]
        x_norm = Q_Q_Prob(x, offset)
        adata.layers["rankit"][i, non_zero] = x_norm.transpose()
    
    return adata
        

def Q_Q_Prob(data, offset):
    
    ranks = sc.stats.rankdata(data, method="dense")
    
    n = len(data)
    prob_level = []
    
    for i in ranks:
        prob_level.append(np.round((i+1-0.5)/n,5))
        
    Standard_normal_quantiles = sc.stats.norm.ppf(prob_level, loc=offset)
    
    return Standard_normal_quantiles

