import anndata as ad
import sys
import os


CAT = "cell_type"
MAX_CELLS = 100000

input_file = sys.argv[1]
output_file = sys.argv[2]

adata = ad.read(input_file)

# Get frac of cells to reach MAX_CELLS
frac = MAX_CELLS / adata.n_obs

if frac < 1:
    uniq_category = adata.obs[CAT].drop_duplicates()

    new_cells = []
    for i in uniq_category:
        new_cells+=list(adata.obs.loc[adata.obs[CAT] == i,:].sample(frac=frac, axis=0).index)
        
    adata = adata[new_cells,]
    
adata.write_h5ad(output_file, compression="gzip")
