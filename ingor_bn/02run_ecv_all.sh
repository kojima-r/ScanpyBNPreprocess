#!/bin/sh
# Backward-compatible wrapper: 06run_ecv_all.sh "2".
sh "$(dirname "$0")/run_ecv_all.sh" "" ../data03_bbknn_p_tissue_t/all.txt
