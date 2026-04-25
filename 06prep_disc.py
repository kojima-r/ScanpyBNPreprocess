"""
Prepare discretized BN-input subsets for the small-scale backends
(MatDNF, ilp_bn, pgmpy) and the discretized FastBN path.

Reads the merged discretized matrices produced by 05disc.py and
strips any "@name" / "tissue" annotation columns, then writes:

  <out-dir>/all_disc.tsv           (full)
  <out-dir>/all_disc{10,100,1000}.tsv  (first N gene columns)
  <out-dir>/tissue/<tissue>.tsv    (per-tissue rows, full columns)

The per-tissue split is only emitted when the input carries a
"tissue" column (i.e. the row-merged branch); the column-merged
(--transposed) branch has no tissue annotation and that step is
skipped.

If a ternary input is also present, the same set of files is written
with a "_tri" suffix (and into a "tissue_tri/" subdirectory).
"""

import argparse
import os

import pandas as pd


SUBSET_SIZES = (10, 100, 1000)


def prep(input_path, out_dir, suffix=""):
    if not os.path.exists(input_path):
        print(f"skip (missing): {input_path}")
        return

    os.makedirs(out_dir, exist_ok=True)
    df_full = pd.read_csv(input_path, sep="\t")
    tissue_col = df_full["tissue"] if "tissue" in df_full.columns else None
    drop_cols = [c for c in ("@name", "tissue") if c in df_full.columns]
    df = df_full.drop(columns=drop_cols)

    # Size-tiered subsets.
    full_path = os.path.join(out_dir, f"all_disc{suffix}.tsv")
    df.to_csv(full_path, index=False, sep="\t")
    print(f">> {full_path}")
    for n in SUBSET_SIZES:
        if n >= df.shape[1]:
            continue
        path = os.path.join(out_dir, f"all_disc{suffix}{n}.tsv")
        df.iloc[:, :n].to_csv(path, index=False, sep="\t")
        print(f">> {path}")

    # Per-tissue full-column tables (only when a tissue column is available).
    if tissue_col is None:
        return
    tissue_dir = os.path.join(out_dir, f"tissue{suffix}")
    os.makedirs(tissue_dir, exist_ok=True)
    for t in tissue_col.unique():
        sub = df[tissue_col.values == t]
        path = os.path.join(tissue_dir, f"{t}.tsv")
        sub.to_csv(path, index=False, sep="\t")
        print(f">> {path}")


SOURCES = ("r", "p")
LEVELS = ("tissue", "age", "batch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=SOURCES, default="r",
                        help="r=resample, p=pseudo_bulk (default-path generation)")
    parser.add_argument("--level", choices=LEVELS, default="batch",
                        help="Stratification level used in step 02 (default: batch)")
    parser.add_argument("--transposed", action="store_true",
                        help="Read 05disc.py outputs from the --transposed branch")
    parser.add_argument("--bin-input", default=None,
                        help="Override: binary discretized merged matrix from 05disc.py")
    parser.add_argument("--tri-input", default=None,
                        help="Override: ternary discretized merged matrix from 05disc.py")
    parser.add_argument("--out-dir", default=None,
                        help="Override: where to write the subsetted .tsv files")
    args = parser.parse_args()

    suffix = "_t" if args.transposed else ""
    base = f"data03_bbknn_{args.source}_{args.level}{suffix}"
    bin_in  = args.bin_input or f"{base}/all_disc.txt"
    tri_in  = args.tri_input or f"{base}/all_disc_tri.txt"
    out_dir = args.out_dir   or f"data04_bbknn_{args.source}_{args.level}{suffix}_disc/"

    prep(bin_in, out_dir, suffix="")
    prep(tri_in, out_dir, suffix="_tri")


if __name__ == "__main__":
    main()
