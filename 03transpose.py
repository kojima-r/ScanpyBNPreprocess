"""
Transpose each (cells x genes) matrix into a (genes x cells)
matrix, the standard layout expected by the downstream BN
estimator.

The script transposes every file matching --input-glob into the
matching basename in --out-dir. Pass it twice if you need to handle
both the resample and pseudo-bulk outputs.
"""

import argparse
import glob
import os

import pandas as pd


def transpose_dir(input_glob, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for in_path in sorted(glob.glob(input_glob)):
        out_path = os.path.join(out_dir, os.path.basename(in_path))
        print(f"{in_path} >> {out_path}")
        df = pd.read_csv(in_path, sep="\t").transpose()
        df.to_csv(out_path, sep="\t", header=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-glob", required=False,
                        help="Glob pattern for input files")
    parser.add_argument("--out-dir", required=False,
                        help="Output directory")
    args = parser.parse_args()

    if args.input_glob and args.out_dir:
        transpose_dir(args.input_glob, args.out_dir)
    else:
        # Backward-compatible default: transpose both pipelines.
        transpose_dir("02data_bbknn/*.txt",  "02data_bbknn_t/")
        transpose_dir("02data_bbknn2/*.txt", "02data_bbknn2_t/")


if __name__ == "__main__":
    main()
