"""
Bootstrap-style resampling of per-tissue expression matrices.

For each input file, draw N resamples of size m (either the original
number of cells, "same", or its square root, "root"), average the
sampled rows, and write one row per resample to the output.

The cell identifier in the input is "tissue|age|batch|cell_id".
--level controls how rows are stratified before resampling (same
convention as 02pseudo_bulk.py — the leading N "|"-tokens of @name
form the stratification key):

  --level tissue  -> "<tissue>",                 output "s<i>|<tissue>"
  --level age     -> "<tissue>|<age>",           output "s<i>|<tissue>|<age>"
  --level batch   -> "<tissue>|<age>|<batch>",   output "s<i>|<tissue>|<age>|<batch>"

The default output directory is data02_bbknn_r_<level>/.
"""

import argparse
import glob
import os
from functools import partial
from multiprocessing import Pool

import numpy as np


LEVEL_DEPTH = {"tissue": 1, "age": 2, "batch": 3}


def _row_key(name, level):
    return "|".join(name.split("|")[:LEVEL_DEPTH[level]])


def resample_size(n, mode):
    return int(np.sqrt(n)) if mode == "root" else n


def resample(in_path, out_path, n_resamples, size_mode, level):
    keys = []
    rows = []
    with open(in_path) as fp:
        header = next(fp)
        for line in fp:
            arr = line.split("\t")
            keys.append(_row_key(arr[0], level))
            rows.append(np.array(list(map(float, arr[1:]))))
    keys = np.array(keys)
    data = np.array(rows)
    unique_keys = np.unique(keys)
    m = resample_size(len(keys), size_mode)

    with open(out_path, "w") as ofp:
        ofp.write(header)
        for i in range(n_resamples):
            for key in unique_keys:
                pool = np.where(keys == key)[0]
                idx = np.random.choice(pool, size=m, replace=True)
                v = np.mean(data[idx, :], axis=0)
                ofp.write(f"s{i}|{key}\t" + "\t".join(map(str, v)) + "\n")


def _worker(in_path, out_dir, n_resamples, size_mode, level):
    out_path = os.path.join(out_dir, os.path.basename(in_path))
    print(f"{in_path} >> {out_path}")
    resample(in_path, out_path, n_resamples, size_mode, level)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-glob", default="data01_bbknn/*.txt",
                        help="Glob pattern for input per-tissue files")
    parser.add_argument("--level", choices=list(LEVEL_DEPTH), default="tissue",
                        help="Stratification level: tissue / age / batch (default: tissue)")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory (default: data02_bbknn_r_<level>/)")
    parser.add_argument("-n", "--n-resamples", type=int, default=10,
                        help="Number of resamples per stratum per input file")
    parser.add_argument("--size", choices=["same", "root"], default="same",
                        help="Resample size: same as input or its square root")
    parser.add_argument("--workers", type=int, default=16,
                        help="Process-pool size")
    args = parser.parse_args()

    out_dir = args.out_dir or f"data02_bbknn_r_{args.level}/"
    os.makedirs(out_dir, exist_ok=True)

    files = sorted(glob.glob(args.input_glob))
    job = partial(_worker,
                  out_dir=out_dir,
                  n_resamples=args.n_resamples,
                  size_mode=args.size,
                  level=args.level)
    with Pool(args.workers) as pool:
        pool.map(job, files)


if __name__ == "__main__":
    main()
