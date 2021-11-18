in_dir="./data/subsampled_100k_cells/"
out_dir="./data/normalized/rankit/"

mkdir -p ${out_dir}

for i in ${in_dir}/*
do
    echo "Working with $(basename $i)"
    python3 ./scripts/normalizing/rankit.py $i ${out_dir}/$(basename $i) &
done

# wait for all pids
ecode=0
while true; do
    [ -z "$(jobs)" ] && break
    wait -n
    err="$?"
    [ "$err" != "0" ] && ecode="$err"
done
