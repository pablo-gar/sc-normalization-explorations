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

## Do scvi normalization

```
mkdir -p ./data/normalized/scvi

python3 ./scripts/normalizing/scvi.py
```


# Getting data ready for R notebook


- Save individual files for mtx, obs and var

```
in_dir="./data/for_notebook/normalized_standard_cell_types/"
out_dir_prefix="./data/for_notebook/matrix_files/"

for i in $in_dir/*
do
    out_dir=$out_dir_prefix/$(basename $i)/
    mkdir -p $out_dir
    python3 ./scripts/processing/save_mtx_obs_var.py rankit $i $out_dir &
done
```
- Make rds file for notebook

```
working_dir="./data/for_notebook/matrix_files/"
ouf_file="./data/for_notebook/all_data.rds"

Rscript ./scripts/processing/make_rds_notebook.R $working_dir $ouf_file
```
