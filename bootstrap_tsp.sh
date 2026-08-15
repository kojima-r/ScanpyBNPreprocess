#!/bin/sh
set -eu

# Tabula Sapiens (Smart-seq) example: row-oriented discrete pipeline.
python 01preprocess_tsp.py --mode smartseq
python 02resample.py --target ss --input-glob "data01_ss/*.txt" --level tissue -n 10
python 04merge.py --target ss --source r --level tissue
python 05disc.py --target ss --source r --level tissue
python 06prep_disc.py --target ss --source r --level tissue
