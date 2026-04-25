"""
Learn a Boolean Bayesian-network skeleton from a binary-discretized
expression table by fitting one MatDNFClassifier per gene.

Input:  TSV with a header row and 0/1 cells (e.g. all_disc100.tsv
        produced by ../06prep_disc.py).
Output: JSON with the dependency parents discovered for each gene:

    {"<child>": ["<parent_1>", "<parent_2>", ...], ...}

The intent is to mirror the small-scale, discretized-input contract
shared with the ilp_bn and pgmpy backends so all three are
interchangeable from the same prep step.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def learn_parents(df: pd.DataFrame, *, h: int, max_itr: int, device: str,
                  threshold: int) -> dict[str, list[str]]:
    from mat_dnf.learner import MatDNFClassifier  # Python 3.10+ runtime
    cols = list(df.columns)
    X_full = df.values.astype(np.int64)
    parents: dict[str, list[str]] = {}

    for j, target in enumerate(cols):
        feature_idx = [i for i in range(len(cols)) if i != j]
        X = X_full[:, feature_idx]
        y = X_full[:, j]

        if len(np.unique(y)) < 2:
            parents[target] = []
            continue

        clf = MatDNFClassifier(h=h, max_itr=max_itr, device_name=device)
        clf.fit(X, y, feature_names=[cols[i] for i in feature_idx],
                target_name=target)
        # Variables that appear in the learned DNF for this target.
        dep = clf.get_supported_vars(X, threshold=threshold)
        parents[target] = sorted({cols[feature_idx[i]] for i in np.where(dep)[0]})
        print(f"{target} <- {parents[target]}")

    return parents


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Binary-discretized TSV (header row, 0/1 values)")
    parser.add_argument("--out", default="./",
                        help="Output directory of <input stem>_matdnf.json (default: ./)")
    parser.add_argument("--h", dest="h", type=int, default=1000,
                        help="Max number of disjuncts in each DNF")
    parser.add_argument("--max-itr", type=int, default=500,
                        help="Max training iterations per target")
    parser.add_argument("--threshold", type=int, default=0,
                        help="Disjunct support threshold for picking a parent")
    parser.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    args = parser.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.out) + in_path.stem + "_matdnf.json")

    df = pd.read_csv(in_path, sep="\t")
    parents = learn_parents(df, h=args.h, max_itr=args.max_itr,
                            device=args.device, threshold=args.threshold)

    with open(out_path, "w") as fp:
        json.dump(parents, fp, indent=2)
    print(f">> {out_path}")


if __name__ == "__main__":
    main()
