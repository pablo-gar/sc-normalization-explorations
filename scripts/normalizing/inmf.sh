in_dir="./data/normalized/rankit/"
out_dir="./data/normalized/inmf/"

mkdir -p ${out_dir}

for i in ${in_dir}/*
do
    echo "Working with $(basename $i)"
    python3 ./scripts/normalizing/inmf.py $i ${out_dir}/$(basename $i) &
done

# wait for all pids
ecode=0
while true; do
    [ -z "$(jobs)" ] && break
    wait -n
    err="$?"
    [ "$err" != "0" ] && ecode="$err"
done
