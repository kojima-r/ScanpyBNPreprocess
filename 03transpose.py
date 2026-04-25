"""
Transpose each (cells x genes) matrix into a (genes x cells)
matrix, the standard layout expected by the downstream BN
estimator.

Default I/O is derived from --source / --level so the directory
naming stays consistent with steps 02 and 04+:

  data02_bbknn_<source>_<level>/    →   data02_bbknn_<source>_<level>_t/

Pass --input-glob / --out-dir to override.
"""

import argparse
import glob
import os

import pandas as pd


SOURCES = ("r", "p")          # r = resample, p = pseudo_bulk
LEVELS = ("tissue", "age", "batch")


def transpose_dir(input_glob, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for in_path in sorted(glob.glob(input_glob)):
        out_path = os.path.join(out_dir, os.path.basename(in_path))
        print(f"{in_path} >> {out_path}")
        df = pd.read_csv(in_path, sep="\t").transpose()
        df.to_csv(out_path, sep="\t", header=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=SOURCES,
                        help="r=resample, p=pseudo_bulk (default-path generation)")
    parser.add_argument("--level", choices=LEVELS,
                        help="Stratification level used in step 02")
    parser.add_argument("--input-glob", default=None,
                        help="Override: glob pattern for input files")
    parser.add_argument("--out-dir", default=None,
                        help="Override: output directory")
    args = parser.parse_args()

    if args.input_glob and args.out_dir:
        transpose_dir(args.input_glob, args.out_dir)
        return

    if not (args.source and args.level):
        parser.error("either (--source AND --level) or (--input-glob AND --out-dir) is required")

    in_glob = args.input_glob or f"data02_bbknn_{args.source}_{args.level}/*.txt"
    out_dir = args.out_dir   or f"data02_bbknn_{args.source}_{args.level}_t/"
    transpose_dir(in_glob, out_dir)


if __name__ == "__main__":
    main()
