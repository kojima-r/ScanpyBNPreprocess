#!/bin/sh
# Bootstrap-style BN structure estimation with `ingor`.
#
# Usage:  sh 05run.sh <input.txt> [N] [out_prefix]
#
#   <input.txt>   : merged matrix produced by 04merge.py
#   N             : ingor sample count (default 10; used 100 for "all2")
#   out_prefix    : where to write per-bootstrap result.ing files
#                   (default: bs_<basename of input>/result.ing)
set -e

input="$1"
N="${2:-10}"
data=$(basename "$input" .txt)
out_prefix="${3:-bs_${data}/result.ing}"

mkdir -p "$(dirname "$out_prefix")"

for B in $(seq 1 10); do
    ingor -B "$B" -N "$N" --single-file off -o "$out_prefix" "$input" &
done
wait
