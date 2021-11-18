import anndata as ad
import sys

cell_type_mapping_table_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

# Read mappings
cell_types = {}
with open(cell_type_mapping_table_file, "r") as cell_type_mapping_table:
    for line in cell_type_mapping_table:
        line = line.rstrip().split("\t")
        cell_types[line[0]] = line[1]
        
# Replace in adata
adata = ad.read(input_file)
adata.obs["cell_type_standard"] = adata.obs["cell_type"]
adata.obs["cell_type_standard"].replace(cell_types, inplace=True)

# Save data 
adata.write_h5ad(output_file, compression="gzip")
