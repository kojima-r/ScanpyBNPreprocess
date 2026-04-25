"""
Discretize a merged expression matrix into binary and ternary
versions for Bayesian-network input.

Accepts either of the two layouts emitted by 04merge.py:

* default (row-merged, samples × genes; the input may carry a
  "tissue" column from 04merge — it is passed through as-is) —
  reads data03_bbknn_<source>_<level>/all.txt.
* --transposed (column-merged, genes × samples) — reads
  data03_bbknn_<source>_<level>_t/all.txt and transposes it to
  samples × genes for discretization. Column headers from the
  merged file become @name values verbatim.

For each gene column, compute the 0.1% and 75% quantiles and bin
values into:

* binary  : 1 if value > 0.1% quantile else 0
* ternary : 0 below 0.1%, 1 between 0.1% and 75%, 2 above 75%
            (collapses to a binary column when the two quantiles tie)

Constant gene columns are dropped. @name is left untouched and no
"tissue" column is added by this step.
"""

import argparse
import os

import pandas as pd

EPS = 1.0e-10

META_COLS = ("@name", "tissue")


def _gene_columns(df):
    return [c for c in df.columns if c not in META_COLS]


def discretize_ternary(df, q_lo=0.001, q_hi=0.75):
    out = df.copy()
    for col in _gene_columns(out):
        q = out[col].quantile([q_lo, q_hi])
        a = q[q_lo] + EPS
        b = q[q_hi] + EPS
        if a == b:
            out[col] = pd.cut(out[col],
                              bins=[-float("inf"), a, float("inf")],
                              labels=[0, 1], right=False, include_lowest=True)
        else:
            out[col] = pd.cut(out[col],
                              bins=[-float("inf"), a, b, float("inf")],
                              labels=[0, 1, 2], right=False, include_lowest=True)
    return out


def to_binary(df_tri):
    """Map ternary {0,1,2} gene columns to binary {0,1} (any non-zero → 1)."""
    out = df_tri.copy()
    for col in _gene_columns(out):
        out[col] = (df_tri[col].astype(int) > 0).astype(int)
    return out


SOURCES = ("r", "p")
LEVELS = ("tissue", "age", "batch")


def load_merged(in_path, transposed):
    """Return a (samples × genes) frame with @name as the first column."""
    if not transposed:
        return pd.read_csv(in_path, sep="\t")

    raw = pd.read_csv(in_path, sep="\t").set_index("@name").T
    raw.index.name = "@name"
    return raw.reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=SOURCES, default="r",
                        help="r=resample, p=pseudo_bulk (default-path generation)")
    parser.add_argument("--level", choices=LEVELS, default="batch",
                        help="Stratification level used in step 02 (default: batch)")
    parser.add_argument("--transposed", action="store_true",
                        help="Read column-merged (genes×samples) output of 04merge.py --transposed")
    parser.add_argument("--input", default=None,
                        help="Override: merged matrix from 04merge.py")
    parser.add_argument("--out-dir", default=None,
                        help="Override: where to write *_disc.txt outputs")
    args = parser.parse_args()

    suffix = "_t" if args.transposed else ""
    base = f"data03_bbknn_{args.source}_{args.level}{suffix}"
    in_path  = args.input   or f"{base}/all.txt"
    out_dir  = args.out_dir or f"{base}/"

    os.makedirs(out_dir, exist_ok=True)
    df = load_merged(in_path, args.transposed)
    df_tri = discretize_ternary(df)

    constant = [c for c in _gene_columns(df_tri) if df_tri[c].nunique() == 1]
    if constant:
        print(f"dropping constant columns: {constant}")
    df_tri = df_tri.drop(columns=constant)

    df_bin = to_binary(df_tri)

    bin_path = os.path.join(out_dir, "all_disc.txt")
    tri_path = os.path.join(out_dir, "all_disc_tri.txt")
    df_bin.to_csv(bin_path, sep="\t", index=False)
    df_tri.to_csv(tri_path, sep="\t", index=False)
    print(f">> {bin_path}")
    print(f">> {tri_path}")


if __name__ == "__main__":
    main()
