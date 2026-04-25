#!/bin/sh
# Small-scale, discretized BN via ILP (PuLP / CBC).
#
# Reads the canonical binary subset staged by ../06prep_disc.py.
set -e

input="${1:-../03data_bbknn_b/all_disc100.tsv}"
python bn_ilp_pulp.py "$input" \
  --max-parents 2 --score bdeu --ess 1.0 \
  --time-limit 3600 --delimiter $'\t'

