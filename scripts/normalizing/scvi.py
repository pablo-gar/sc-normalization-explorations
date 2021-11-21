import scvi
import scanpy as sc
import os

out_dir = "./data/normalized/scvi/"

out_file_single = os.path.join(out_dir, "lung_map.h5ad")
out_file_concatenated = os.path.join(out_dir, "all.h5ad")

adata_1 = sc.read("./data/normalized/inmf/lung_map.h5ad")
adata_2 = sc.read("./data/normalized/inmf/tabula_sapiens.h5ad")
adata_3 = sc.read("./data/normalized/inmf/zemin_covid.h5ad")

adata_single = adata_1.copy()
adatas = adata_1.concatenate(adata_2)
adatas = adatas.concatenate(adata_3)

shared_columns = list(set(adata_1.obs.columns).intersection(set(adata_2.obs.columns)).intersection(set(adata_3.obs.columns)))

adatas.obs = adatas.obs[shared_columns + ["batch"]]

#------
# Run scvi

scvi.model.SCVI.setup_anndata(adatas, batch_key="batch")
scvi.model.SCVI.setup_anndata(adata_single)

model_single = scvi.model.SCVI(adata_single, n_layers=2, n_latent=30)
model_concatenaded = scvi.model.SCVI(adatas, n_layers=2, n_latent=30)

model_single.train()
model_concatenaded.train()

#------
# Get normalized expression
adata_single.layers["scvi"] = model_single.get_normalized_expression(adata_single, return_numpy=True)
adatas.layers["scvi"] = model_concatenaded.get_normalized_expression(adatas, return_numpy=True)

#------
# Write
adata_single.write(out_file_single, compression="gzip")
adatas.write(out_file_concatenated, compression="gzip")
