#!/bin/sh
# Run BN structure estimation for every per-tissue file under 02data_bbknn_t/.
set -e
for f in 02data_bbknn_t/*.txt; do
    echo "$f"
    sh "$(dirname "$0")/run_cmd.sh" "$f"
done
