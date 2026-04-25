#!/bin/sh
# Per-tissue BN bootstrap + final structure aggregation.
#
# Usage:  sh run_cmd.sh <data_path>
#   <data_path>: a transposed per-tissue matrix (e.g. 02data_bbknn_t/Aorta.txt)
set -e

data_path="$1"
data=$(basename "$data_path" .txt)
out_dir="bs_${data}"

mkdir -p "$out_dir"

for B in $(seq 1 10); do
    ingor -B "$B" -N 10 --single-file off -o "${out_dir}/result.ing" "$data_path" &
done
wait

ingor --bs prefix=${out_dir}/result.ing,type=ing,ed=1000,th=0.05 -o "result_${data}.txt"
