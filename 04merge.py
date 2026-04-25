"""
Merge per-tissue matrices produced by the previous step into one
combined matrix.

Two layouts are supported:

* default: matrices are joined column-wise (samples concatenated
  across tissues). Each column is renamed to "<tissue>_<sample>" so
  the originating tissue is preserved.
* --batched: matrices are joined row-wise. A "tissue" column is
  added to record the originating file.

Default I/O is derived from --source / --level so the directory
naming stays consistent with steps 02 and 03:

  default:    02data_bbknn_<source>_<level>_t/   →  03data_bbknn_<source>_<level>/all.txt
  --batched:  02data_bbknn_<source>_<level>/     →  03data_bbknn_b_<source>_<level>/all.txt

Pass --input-glob / --out to override.
"""

import argparse
import glob
import os

import pandas as pd


SOURCES = ("r", "p")
LEVELS = ("tissue", "age", "batch")


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
    parser.add_argument("--source", choices=SOURCES,
                        help="r=resample, p=pseudo_bulk (default-path generation)")
    parser.add_argument("--level", choices=LEVELS,
                        help="Stratification level used in step 02")
    parser.add_argument("--batched", action="store_true",
                        help="Stack rows (with a 'tissue' column) instead of columns")
    parser.add_argument("--input-glob", default=None,
                        help="Override: glob pattern for input files")
    parser.add_argument("--out", default=None,
                        help="Override: output file path")
    args = parser.parse_args()

    fn = merge_rows if args.batched else merge_columns

    if args.input_glob and args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        fn(args.input_glob, args.out)
        return

    if not (args.source and args.level):
        parser.error("either (--source AND --level) or (--input-glob AND --out) is required")

    if args.batched:
        in_glob = args.input_glob or f"02data_bbknn_{args.source}_{args.level}/*.txt"
        out    = args.out        or f"03data_bbknn_b_{args.source}_{args.level}/all.txt"
    else:
        in_glob = args.input_glob or f"02data_bbknn_{args.source}_{args.level}_t/*.txt"
        out    = args.out        or f"03data_bbknn_{args.source}_{args.level}/all.txt"

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fn(in_glob, out)


if __name__ == "__main__":
    main()
