#!/bin/sh
# Learn a small-scale Boolean BN skeleton via MatDNF.
#
# Reads the canonical binary subsets staged by ../06prep_disc.py
# (in ../data04_bbknn_r_batch_disc/) and writes the parent map next to the input.
set -e

input="${1:-../data04_bbknn_p_tissue_disc/all_disc10.tsv}"
python learn_bn_matdnf.py "$input" --h 1000 --max-itr 500
