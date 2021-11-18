import anndata as ad
import sys
import os

tissue = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

ASSAY = "10x"

adata = ad.read(input_file)
adata = adata[adata.obs["tissue"] == tissue,:]
adata = adata[adata.obs["assay"].str.contains(ASSAY),:]
adata.write(output_file, compression="gzip")
