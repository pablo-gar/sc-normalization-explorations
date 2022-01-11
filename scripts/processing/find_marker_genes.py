from scipy import sparse
import numpy 
import numpy.ma as ma
import scanpy as sc
import pandas as pd
import anndata as ad
import sys

layer = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

adata = ad.read(input_file)
adata.var_names = adata.var["feature_name"]

# Make data that goes 0-1 logarithmic
if adata.layers[layer].max() <= 1:
    if sparse.issparse(adata.layers[layer]):
        min_val = adata.layers[layer].data.min()
        adata.layers[layer].data = numpy.log10(adata.layers[layer].data * (1/min_val))
    else:
        min_val = ma.masked_array(adata.layers[layer], mask = adata.layers[layer] != 0)
        min_val = min_val.min()
        adata.layers[layer] = numpy.log10(adata.layers[layer] * (1/min_val))

results = []
for cell_type in adata.obs.cell_type.drop_duplicates():
    
    print("processing", cell_type)
    # Create a mask for the cell type
    adata.obs["_mask"] = adata.obs.cell_type == cell_type
    adata.obs["_mask"] = adata.obs["_mask"].astype("category")
    
    # Find top 1000 differentially expressed genes
    sc.tl.rank_genes_groups(adata, groupby="_mask", n_genes=1000, layer=layer, rankby_abs=True)
    
    # Rearrange and convert to pandas
    gene_names = pd.DataFrame(adata.uns["rank_genes_groups"]["names"])
    folds = pd.DataFrame(adata.uns["rank_genes_groups"]["logfoldchanges"])
    
    gene_names.rename(columns={"True": "gene"}, inplace=True)
    folds.rename(columns={"True": "foldchange"}, inplace=True)
    
    gene_folds = pd.concat([gene_names, folds], axis=1)
    
    # Find the top 3
    gene_folds.sort_values("foldchange", axis=0, inplace=True, ascending=False)
    gene_folds = gene_folds.iloc[:2,:] 
    gene_folds["cell_type"]  = cell_type
    results.append(gene_folds.loc[:,["gene", "foldchange", "cell_type"]])

results = pd.concat(results, axis=0)

results.to_csv(output_file, sep="\t", header=True, index=False)
