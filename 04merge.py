"""
Merge per-tissue matrices produced by the previous step into one
combined matrix.

Two layouts are supported:

* default: matrices are joined row-wise (samples stacked across
  tissues). A "tissue" column is added to record the originating
  file. Reads the (cells×genes) outputs of step 02 directly.
* --transposed: matrices are joined column-wise (samples concatenated
  across tissues). Each column is renamed to "<tissue>_<sample>" so
  the originating tissue is preserved. Reads the (genes×cells)
  outputs of step 03transpose.

Default I/O is derived from --target / --source / --level so the
directory naming stays consistent with steps 02 and 03:

  default:      data02_<target>_<source>_<level>/    →  data03_<target>_<source>_<level>/all.txt
  --transposed: data02_<target>_<source>_<level>_t/  →  data03_<target>_<source>_<level>_t/all.txt

--target defaults to bbknn; pass facs / tsp / droplet (or any custom
name) to switch pipelines. Pass --input-glob / --out to override.
"""

import argparse
import glob
import os

import pandas as pd


SOURCES = ("r", "p")
LEVELS = ("tissue", "age", "batch")


def merge_columns(input_glob, out_path, min_var=0.0):
    frames = []
    for in_path in sorted(glob.glob(input_glob)):
        print(in_path)
        name = os.path.basename(in_path)
        df = pd.read_csv(in_path, sep="\t", index_col=0)
        df.columns = [f"{name}_{c}" for c in df.columns]
        frames.append(df)
    merged = pd.concat(frames, axis=1, join="inner")
    if min_var > 0:
        before = merged.shape[0]
        keep = merged.var(axis=1) > min_var
        merged = merged.loc[keep]
        print(f"   variance filter (>{min_var}): {before} -> {merged.shape[0]} rows")
    merged.to_csv(out_path, sep="\t", index=True)
    print(f">> {out_path} ({merged.shape})")


def merge_rows(input_glob, out_path, min_var=0.0):
    frames = []
    for in_path in sorted(glob.glob(input_glob)):
        print(in_path)
        df = pd.read_csv(in_path, sep="\t", index_col=0)
        name, _ = os.path.splitext(os.path.basename(in_path))
        df["tissue"] = name
        frames.append(df)
    merged = pd.concat(frames, axis=0, join="inner")
    if min_var > 0:
        gene_cols = merged.columns.drop("tissue")
        before = len(gene_cols)
        keep = merged[gene_cols].var(axis=0) > min_var
        kept_cols = gene_cols[keep]
        merged = merged[list(kept_cols) + ["tissue"]]
        print(f"   variance filter (>{min_var}): {before} -> {len(kept_cols)} cols")
    merged.to_csv(out_path, sep="\t", index=True)
    print(f">> {out_path} ({merged.shape})")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default="bbknn",
                        help="Dataset slot used in default I/O paths "
                             "(data02_<target>_..., data03_<target>_...). "
                             "Default: bbknn.")
    parser.add_argument("--source", choices=SOURCES,
                        help="r=resample, p=pseudo_bulk (default-path generation)")
    parser.add_argument("--level", choices=LEVELS,
                        help="Stratification level used in step 02")
    parser.add_argument("--transposed", action="store_true",
                        help="Stack rows (with a 'tissue' column) instead of columns")
    parser.add_argument("--input-glob", default=None,
                        help="Override: glob pattern for input files")
    parser.add_argument("--out", default=None,
                        help="Override: output file path")
    parser.add_argument("--min-var", type=float, default=0.0,
                        help="Drop genes whose variance is <= this threshold "
                             "after merging (0 = disabled, default).")
    args = parser.parse_args()

    fn = merge_columns if args.transposed else merge_rows

    if args.input_glob and args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        fn(args.input_glob, args.out, min_var=args.min_var)
        return

    if not (args.source and args.level):
        parser.error("either (--source AND --level) or (--input-glob AND --out) is required")

    if args.transposed:
        in_glob = args.input_glob or f"data02_{args.target}_{args.source}_{args.level}_t/*.txt"
        out    = args.out        or f"data03_{args.target}_{args.source}_{args.level}_t/all.txt"
    else:
        in_glob = args.input_glob or f"data02_{args.target}_{args.source}_{args.level}/*.txt"
        out    = args.out        or f"data03_{args.target}_{args.source}_{args.level}/all.txt"

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fn(in_glob, out, min_var=args.min_var)


if __name__ == "__main__":
    main()
