#!/bin/sh
set -eu

# Tabula Muris Senis (BBKNN) example: row-oriented discrete pipeline.
python 01preprocess_tms.py --mode bbknn
python 02resample.py --target bbknn --input-glob "data01_bbknn/*.txt" --level tissue -n 10
python 04merge.py --target bbknn --source r --level tissue
python 05disc.py --target bbknn --source r --level tissue
python 06prep_disc.py --target bbknn --source r --level tissue
