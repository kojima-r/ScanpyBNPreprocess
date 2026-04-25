"""
Merge per-tissue matrices produced by the previous step into one
combined matrix.

Two layouts are supported:

* default: matrices are joined column-wise (samples concatenated
  across tissues). Each column is renamed to "<tissue>_<sample>" so
  the originating tissue is preserved.
* --batched: matrices are joined row-wise. A "tissue" column is
  added to record the originating file.
"""

import argparse
import glob
import os

import pandas as pd


def merge_columns(input_glob, out_path):
    frames = []
    for in_path in sorted(glob.glob(input_glob)):
        print(in_path)
        name = os.path.basename(in_path)
        df = pd.read_csv(in_path, sep="\t", index_col=0)
        df.columns = [f"{name}_{c}" for c in df.columns]
        frames.append(df)
    merged = pd.concat(frames, axis=1, join="inner")
    merged.to_csv(out_path, sep="\t", index=True)
    print(f">> {out_path} ({merged.shape})")


def merge_rows(input_glob, out_path):
    frames = []
    for in_path in sorted(glob.glob(input_glob)):
        print(in_path)
        df = pd.read_csv(in_path, sep="\t", index_col=0)
        name, _ = os.path.splitext(os.path.basename(in_path))
        df["tissue"] = name
        frames.append(df)
    merged = pd.concat(frames, axis=0, join="inner")
    merged.to_csv(out_path, sep="\t", index=True)
    print(f">> {out_path} ({merged.shape})")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-glob", default=None,
                        help="Glob pattern for input files")
    parser.add_argument("--out", default=None,
                        help="Output file path")
    parser.add_argument("--batched", action="store_true",
                        help="Stack rows (with a 'tissue' column) instead of columns")
    args = parser.parse_args()

    if args.input_glob and args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        (merge_rows if args.batched else merge_columns)(args.input_glob, args.out)
        return

    if args.batched:
        os.makedirs("03data_bbknn_b/", exist_ok=True)
        merge_rows("02data_bbknn_b/*.txt", "03data_bbknn_b/all.txt")
    else:
        os.makedirs("03data_bbknn/", exist_ok=True)
        merge_columns("02data_bbknn_t/*.txt",  "03data_bbknn/all.txt")
        merge_columns("02data_bbknn2_t/*.txt", "03data_bbknn/all2.txt")


if __name__ == "__main__":
    main()
