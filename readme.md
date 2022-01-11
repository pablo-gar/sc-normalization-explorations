# Download data

The following datasets will be used for testing normalization methods. For each only the data from lung will be used.

1. [LungMAP — Human data from a broad age healthy donor group.](https://cellxgene.cziscience.com/collections/625f6bf4-2f33-4942-962e-35243d284837)
2. [Tabula Sapiens - All Cells](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5)
2. [COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas](https://cellxgene.cziscience.com/collections/0a839c4b-10d0-4d64-9272-684c49a2c8ba)


```
bash scripts/download/download_data.sh
```

# Process data

## Select only lung data and 10X

```
bash ./scripts/processing/subsample_lung.sh
```

## Remove cells with <500 expressed genes and make raw data `adata.X`

```
bash ./scripts/processing/drop_low_coverage_cells.sh  
```

## Sub-sample to 100K max cells

This will subsample cells to a max of 100K cell per dataset, the sampling is done per cell type so that all cell tpyes are preserved with the same original proportions

```
 bash ./scripts/processing/subsample_h5ad.sh
```

# Normalize data 

## Do quantile normalization

```
 bash ./scripts/normalizing/rankit.sh
```

## Do iNMF-like normalization normalization

```
bash ./scripts/normalizing/inmf.sh
```

## Do Z-scores

```
bash ./scripts/normalizing/z_scores.sh
```

## Do scvi normalization

```
mkdir -p ./data/normalized/scvi

python3 ./scripts/normalizing/scvi.py
```

# Get marker genes

```
mkdir -p ./data/marker_genes/
python3 ./scripts/processing/find_marker_genes.py rankit ./data/normalized/scvi/all.h5ad ./data/marker_genes/marker_genes_rankit.tsv
python3 ./scripts/processing/find_marker_genes.py scvi ./data/normalized/scvi/all.h5ad ./data/marker_genes/marker_genes_scvi.tsv
```


# Getting data ready for R notebook


Save scvi normalization of data only from lung map

```
mkdir -p ./data/for_notebook/
python3 ./scripts/processing/save_mtx_obs_var.py all ./data/normalized/scvi/lung_map.h5ad ./data/for_notebook/lung_map
```

Save all normalization methods of all datasets

```
mkdir -p ./data/for_notebook/
python3 ./scripts/processing/save_mtx_obs_var.py all ./data/normalized/scvi/all.h5ad ./data/for_notebook/all_datasets
```


