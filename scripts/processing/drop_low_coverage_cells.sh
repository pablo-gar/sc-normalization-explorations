in_dir="./data/subsampled_lung/"
out_dir="./data/subsampled_high_coverage_cells/"
genes="500"

mkdir -p ${out_dir}

for i in ${in_dir}/*
do
    echo "Working with $(basename $i)"
    python3 ./scripts/processing/drop_low_coverage_cells.py $genes $i ${out_dir}/$(basename $i) &
done

# wait for all pids
ecode=0
while true; do
    [ -z "$(jobs)" ] && break
    wait -n
    err="$?"
    [ "$err" != "0" ] && ecode="$err"
done
