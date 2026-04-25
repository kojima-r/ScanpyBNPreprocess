#!/bin/sh
# Backward-compatible wrapper: bootstrap BN with N=100 over data03_bbknn/all2.txt.
sh "$(dirname "$0")/run_bs.sh" ../data03_bbknn_p_tissue_t/all.txt 10 3
