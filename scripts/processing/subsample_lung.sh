in_dir="./data/raw/"
out_dir="./data/subsampled_lung/"

mkdir -p ${out_dir}

for i in ${in_dir}/*
do
    echo "Working with $(basename $i)"
    python3 ./scripts/processing/subsample_tisuse.py lung $i ${out_dir}/$(basename $i) 
done
