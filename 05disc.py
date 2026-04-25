"""
Discretize a merged expression matrix into binary and ternary
versions for Bayesian-network input.

For each gene column, compute the 0.1% and 75% quantiles and bin
values into:

* binary  : 1 if value > 0.1% quantile else 0
* ternary : 0 below 0.1%, 1 between 0.1% and 75%, 2 above 75%
            (collapses to a binary column when the two quantiles tie)

Constant columns are dropped, the "@name" column gets its age token
zero-padded so it sorts correctly, and a "tissue" column is restored
from the input.
"""

import argparse
import os

import pandas as pd

EPS = 1.0e-10


def discretize_ternary(df, q_lo=0.001, q_hi=0.75):
    out = df.copy()
    for col in out.columns[1:]:  # skip @name
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
    """Map ternary {0,1,2} columns to binary {0,1} (any non-zero → 1)."""
    name_col = df_tri[["@name"]]
    body = {col: (df_tri[col] > 0).astype(int) for col in df_tri.columns[1:]}
    return pd.concat([name_col, pd.DataFrame(body, index=df_tri.index)], axis=1)


def _pad_age(name):
    """Zero-pad the first '<n>m' / '<n>d' age token in a '|'-joined name.

    Works for any number of parts, so it accepts both the legacy
    "s<i>|<age>|<batch>" form and the aligned
    "s<i>|<tissue>|<age>|<batch>" form emitted by 02resample.py.
    """
    parts = name.split("|")
    for i, p in enumerate(parts):
        if len(p) >= 2 and p[-1] in ("m", "d") and p[:-1].isdigit():
            parts[i] = f"{int(p[:-1]):02d}{p[-1]}"
            break
    return "|".join(parts)


def postprocess(df, tissue_col):
    constant = [c for c in df.columns if df[c].nunique() == 1]
    if constant:
        print(f"dropping constant columns: {constant}")
    df = df.drop(columns=constant)
    df["tissue"] = tissue_col
    df["@name"] = df["@name"].apply(_pad_age)
    return df.sort_values(by="@name")


SOURCES = ("r", "p")
LEVELS = ("tissue", "age", "batch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=SOURCES, default="r",
                        help="r=resample, p=pseudo_bulk (default-path generation)")
    parser.add_argument("--level", choices=LEVELS, default="batch",
                        help="Stratification level used in step 02 (default: batch)")
    parser.add_argument("--input", default=None,
                        help="Override: merged batched matrix from 04merge.py --batched")
    parser.add_argument("--out-dir", default=None,
                        help="Override: where to write *_disc.txt outputs")
    args = parser.parse_args()

    base = f"03data_bbknn_b_{args.source}_{args.level}"
    in_path  = args.input   or f"{base}/all.txt"
    out_dir  = args.out_dir or f"{base}/"

    os.makedirs(out_dir, exist_ok=True)
    df = pd.read_csv(in_path, sep="\t")

    tissue = df["tissue"]
    df_tri = discretize_ternary(df.drop(columns=["tissue"]))

    constant = [c for c in df_tri.columns if df_tri[c].nunique() == 1]
    if constant:
        print(f"dropping constant columns (pre-binarize): {constant}")
    df_tri = df_tri.drop(columns=constant)

    df_bin = postprocess(to_binary(df_tri), tissue)
    df_tri = postprocess(df_tri, tissue)

    bin_path = os.path.join(out_dir, "all_disc.txt")
    tri_path = os.path.join(out_dir, "all_disc_tri.txt")
    df_bin.to_csv(bin_path, sep="\t", index=False)
    df_tri.to_csv(tri_path, sep="\t", index=False)
    print(f">> {bin_path}")
    print(f">> {tri_path}")


if __name__ == "__main__":
    main()
