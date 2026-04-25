"""
Prepare discretized BN-input subsets for the small-scale backends
(MatDNF, ilp_bn, pgmpy) and the discretized FastBN path.

Reads the merged discretized matrices produced by 05disc.py and
strips the "@name" and "tissue" annotation columns, then writes:

  <out-dir>/all_disc.tsv           (full)
  <out-dir>/all_disc{10,100,1000}.tsv  (first N gene columns)
  <out-dir>/tissue/<tissue>.tsv    (per-tissue rows, full columns)

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
    tissue_col = df_full["tissue"]
    df = df_full.drop(columns=["@name", "tissue"])

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

    # Per-tissue full-column tables.
    tissue_dir = os.path.join(out_dir, f"tissue{suffix}")
    os.makedirs(tissue_dir, exist_ok=True)
    for t in tissue_col.unique():
        sub = df[tissue_col.values == t]
        path = os.path.join(tissue_dir, f"{t}.tsv")
        sub.to_csv(path, index=False, sep="\t")
        print(f">> {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-input", default="03data_bbknn_b/all_disc.txt",
                        help="Binary discretized merged matrix from 05disc.py")
    parser.add_argument("--tri-input", default="03data_bbknn_b/all_disc_tri.txt",
                        help="Ternary discretized merged matrix from 05disc.py")
    parser.add_argument("--out-dir", default="03data_bbknn_b/",
                        help="Where to write the subsetted .tsv files")
    args = parser.parse_args()

    prep(args.bin_input, args.out_dir, suffix="")
    prep(args.tri_input, args.out_dir, suffix="_tri")


if __name__ == "__main__":
    main()
