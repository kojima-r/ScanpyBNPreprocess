"""
Collapse each per-tissue expression matrix into a single
pseudo-bulk row (the mean over all cells).

One output file is written per input file, with a single sample
named "s<basename>".
"""

import argparse
import glob
import os
from functools import partial
from multiprocessing import Pool

import numpy as np


def pseudo_bulk(in_path, out_path):
    with open(in_path) as fp:
        header = next(fp)
        rows = []
        for line in fp:
            arr = line.split("\t")
            rows.append(np.array(list(map(float, arr[1:]))))
    data = np.array(rows)
    name = os.path.basename(in_path)
    v = np.mean(data, axis=0)
    with open(out_path, "w") as ofp:
        ofp.write(header)
        ofp.write(f"s{name}\t" + "\t".join(map(str, v)) + "\n")


def _worker(in_path, out_dir):
    out_path = os.path.join(out_dir, os.path.basename(in_path))
    print(f"{in_path} >> {out_path}")
    pseudo_bulk(in_path, out_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-glob", default="01data_bbknn/*.txt",
                        help="Glob pattern for input per-tissue files")
    parser.add_argument("--out-dir", default="02data_bbknn2/",
                        help="Output directory")
    parser.add_argument("--workers", type=int, default=16,
                        help="Process-pool size")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    files = sorted(glob.glob(args.input_glob))
    with Pool(args.workers) as pool:
        pool.map(partial(_worker, out_dir=args.out_dir), files)


if __name__ == "__main__":
    main()
