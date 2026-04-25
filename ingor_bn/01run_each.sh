#!/bin/sh
# Run BN structure estimation for every per-tissue file under data02_bbknn_t/.
set -e
for f in ../data02_bbknn_p_tissue_t/*.txt; do
    echo "$f"
    sh "$(dirname "$0")/run_cmd.sh" "$f"
done
