#!/bin/sh
# Aggregate bootstrapped BN structures and compute per-tissue ECv networks.
#
# Usage:  sh 06run_ecv_all.sh [tag]
#
#   tag : "" (default) processes bs_all + 02data_bbknn_t  → ecv_all/
#         "2"          processes bs_all2 + 02data_bbknn2_t → ecv_all2/
set -e

tag="$1"
bs_dir="bs_all${tag}"
src_dir="02data_bbknn${tag}_t"
out_dir="ecv_all${tag}"
sgn="result_all${tag}.sgn3"
txt="result_all${tag}.txt"

# Drop empty/failed bootstrap files (smaller than 1KB).
find "./${bs_dir}" -type f -size -1024c -print -delete

mkdir -p "$out_dir"

ingor --bs prefix=${bs_dir}/result.ing,type=ing,ed=1000,th=0.05 -o "$sgn"
ingor --bs prefix=${bs_dir}/result.ing,type=ing,ed=1000,th=0.05 -o "$txt"

for f in "${src_dir}"/*.txt; do
    name=$(basename "$f" .txt)
    echo "$name"
    ingor --read file="$sgn" --ec data="${src_dir}/${name}.txt",method=ECv -o "${out_dir}/${name}.txt"
done
