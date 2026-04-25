#!/bin/sh
# Small-scale, discretized BN via pgmpy (HillClimb + PC hybrid).
#
# Reads the canonical binary subset staged by ../06prep_disc.py.
set -e

input="${1:-../03data_bbknn_b_r_batch/all_disc100.tsv}"
python learn_bn_pgmpy.py "$input" \
  --discretize quantile --bins 3 \
  --estimator hybrid --score k2 \
  --na-policy drop-rows --na-as-bin \
  --output-prefix bn_result \
  --seed 42 \
  --max-indegree 3 --max-iter 100 \
  --var-threshold 0.0
