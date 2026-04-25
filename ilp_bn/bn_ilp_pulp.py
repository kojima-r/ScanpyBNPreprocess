#!/usr/bin/env python3
"""
Bayesian Network structure learning from 0/1 table data using Integer Linear Programming (PuLP).

- Reads a CSV of binary variables (0/1) into pandas.
- Precomputes local scores (BIC or BDeu) for each node and each candidate parent set up to max_parents.
- Formulates an ILP with one binary variable y[j,S] for choosing the parent set S for node j.
- Enforces acyclicity via topological-order variables u[j] with big-M constraints.
- Solves with PuLP's CBC solver (bundled).

Usage:
    python bn_ilp_pulp.py path/to/data.csv --max-parents 2 --score bic --time-limit 30
    python bn_ilp_pulp.py path/to/data.csv --max-parents 2 --score bdeu --ess 1.0 --time-limit 60

If a time limit is set and the solver stops early, the script outputs the best incumbent solution found
up to that point (best-so-far), along with the solver status.
"""

import argparse
import json
import math
import sys
from pathlib import Path
from itertools import combinations, product
from typing import Dict, Iterable, List, Sequence, Tuple, FrozenSet

import numpy as np
import pandas as pd
import pulp


def powerset_up_to_k(items: Sequence[str], k: int) -> Iterable[FrozenSet[str]]:
    for r in range(0, k + 1):
        for comb in combinations(items, r):
            yield frozenset(comb)


def bic_local_score_binary(df: pd.DataFrame, child: str, parents: FrozenSet[str]) -> float:
    """
    Compute the local BIC score for a binary child given a (possibly empty) set of binary parents.
    BIC = log-likelihood - 0.5 * (#params) * log(N)
    For binary variables:
        #params = (2^|parents|) * (2 - 1) = 2^{|parents|}
    """
    N = len(df)
    if N == 0:
        raise ValueError("Empty dataframe.")
    if len(parents) == 0:
        p1 = df[child].mean()
        p1 = min(max(p1, 1e-9), 1 - 1e-9)
        ll = (df[child] * math.log(p1) + (1 - df[child]) * math.log(1 - p1)).sum()
        penalty = 0.5 * (1) * math.log(N)
        return float(ll - penalty)

    parent_cols = list(parents)
    ll = 0.0
    for vals in product([0, 1], repeat=len(parent_cols)):
        mask = np.ones(N, dtype=bool)
        for col, v in zip(parent_cols, vals):
            mask &= (df[col].values == v)
        Nj = int(mask.sum())
        if Nj == 0:
            continue
        y = df.loc[mask, child].values
        p1 = float(y.mean())
        p1 = min(max(p1, 1e-9), 1 - 1e-9)
        ll += float((y * math.log(p1) + (1 - y) * math.log(1 - p1)).sum())

    penalty = 0.5 * (2 ** len(parents)) * math.log(N)
    return float(ll - penalty)


def bdeu_local_score_binary(df: pd.DataFrame, child: str, parents: FrozenSet[str], ess: float = 1.0) -> float:
    """
    Compute BDeu local score for a binary child given binary parents.
    Uniform Dirichlet priors with equivalent sample size (ess).
    r = 2 (child states), q = 2^{|parents|} (parent configurations).
    Score = sum_{parent configs} [ logGamma(alpha_q) - sum_x logGamma(alpha_qx)
                                   + sum_x logGamma(alpha_qx + N_qx) - logGamma(alpha_q + N_q) ]
    where alpha_qx = ess / (r * q), alpha_q = ess / q.
    """
    N = len(df)
    if N == 0:
        raise ValueError("Empty dataframe.")
    r = 2
    q = 2 ** len(parents) if len(parents) > 0 else 1
    alpha_q = ess / q
    alpha_qx = ess / (r * q)

    if len(parents) == 0:
        y = df[child].values
        N1 = int(y.sum())
        N0 = int(len(y) - N1)
        return (math.lgamma(alpha_q) - r * math.lgamma(alpha_qx)
                + math.lgamma(alpha_qx + N0) + math.lgamma(alpha_qx + N1)
                - math.lgamma(alpha_q + (N0 + N1)))

    parent_cols = list(parents)
    score = 0.0
    for vals in product([0, 1], repeat=len(parent_cols)):
        mask = np.ones(N, dtype=bool)
        for col, v in zip(parent_cols, vals):
            mask &= (df[col].values == v)
        Nj = int(mask.sum())
        if Nj == 0:
            score += (math.lgamma(alpha_q) - r * math.lgamma(alpha_qx))
            continue
        y = df.loc[mask, child].values
        N1 = int(y.sum())
        N0 = int(Nj - N1)
        score += (math.lgamma(alpha_q) - r * math.lgamma(alpha_qx)
                  + math.lgamma(alpha_qx + N0) + math.lgamma(alpha_qx + N1)
                  - math.lgamma(alpha_q + Nj))
    return float(score)


def build_and_solve_ilp(
    df: pd.DataFrame,
    max_parents: int = 2,
    time_limit: int = None,
    score: str = "bic",
    ess: float = 1.0,
) -> Dict:
    """
    Build and solve the ILP for Bayesian network structure learning.
    Returns a dict describing the (best found) solution; if a time limit is given, this
    may be the best incumbent instead of a proven optimum.
    """
    variables = list(df.columns)
    n = len(variables)
    if n == 0:
        raise ValueError("No variables in dataframe.")
    if max_parents < 0:
        raise ValueError("max_parents must be >= 0")

    # Precompute candidate parent sets and their local scores
    candidate_sets: Dict[str, List[FrozenSet[str]]] = {}
    local_score: Dict[Tuple[str, FrozenSet[str]], float] = {}
    score = score.lower()
    if score not in {"bic", "bdeu"}:
        raise ValueError("score must be either 'bic' or 'bdeu'")
    for child in variables:
        others = [v for v in variables if v != child]
        cands = [S for S in powerset_up_to_k(others, max_parents)]
        candidate_sets[child] = cands
        for S in cands:
            if score == "bic":
                local_score[(child, S)] = bic_local_score_binary(df, child, S)
            else:
                local_score[(child, S)] = bdeu_local_score_binary(df, child, S, ess=ess)

    # ILP model
    prob = pulp.LpProblem("BN_Structure_Learning", pulp.LpMaximize)

    # y_(child,S) variables: choose exactly one parent set for each child
    y_vars: Dict[Tuple[str, FrozenSet[str]], pulp.LpVariable] = {}
    for child in variables:
        for S in candidate_sets[child]:
            y_vars[(child, S)] = pulp.LpVariable(
                f"y__{child}__{'_'.join(sorted(S)) if S else 'empty'}",
                lowBound=0, upBound=1, cat="Binary"
            )

    # Topological order variables u_j in [0, n-1]
    u_vars: Dict[str, pulp.LpVariable] = {
        v: pulp.LpVariable(f"u__{v}", lowBound=0, upBound=n - 1, cat="Continuous") for v in variables
    }

    # Objective: sum of local scores for chosen parent sets
    prob += pulp.lpSum(local_score[(child, S)] * y_vars[(child, S)]
                       for child in variables for S in candidate_sets[child])

    # For each child, select exactly one parent set
    for child in variables:
        prob += pulp.lpSum(y_vars[(child, S)] for S in candidate_sets[child]) == 1, f"one_parent_set__{child}"

    # Acyclicity via order constraints
    M = n  # safe Big-M since u in [0, n-1]
    for child in variables:
        for S in candidate_sets[child]:
            for parent in S:
                prob += u_vars[parent] + 1 <= u_vars[child] + M * (1 - y_vars[(child, S)]), \
                        f"acyc__{parent}_to_{child}__{'_'.join(sorted(S)) or 'empty'}"

    # Solve with (optional) time limit
    solver = pulp.PULP_CBC_CMD(msg=True, timeLimit=time_limit) if time_limit else pulp.PULP_CBC_CMD(msg=True)
    result_status = prob.solve(solver)
    status_str = pulp.LpStatus[result_status]

    # Extract best-so-far incumbent regardless of status
    selected_parent_sets: Dict[str, List[str]] = {}
    for child in variables:
        # get values; if None, treat as 0
        vals = [(S, (y_vars[(child, S)].varValue or 0.0)) for S in candidate_sets[child]]
        # choose argmax
        best_S = max(vals, key=lambda t: t[1])[0]
        selected_parent_sets[child] = sorted(list(best_S))

    edges = [[p, c] for c, Ps in selected_parent_sets.items() for p in Ps]

    # Robust objective computation from incumbent
    obj_value = 0.0
    for child in variables:
        # find the chosen S
        chosen_S = frozenset(selected_parent_sets[child])
        obj_value += local_score[(child, chosen_S)]

    sol = {
        "status": status_str,
        "objective_value": float(obj_value),
        "order": {v: float(u_vars[v].varValue) if u_vars[v].varValue is not None else None for v in variables},
        "selected_parent_sets": selected_parent_sets,
        "edges": edges,
        "score": score,
        "ess": ess if score == "bdeu" else None,
        "time_limit": time_limit,
    }
    return sol


def main():
    parser = argparse.ArgumentParser(description="BN structure learning via ILP (PuLP) for 0/1 data.")
    parser.add_argument("csv_path", help="Path to CSV file of 0/1 variables.")
    parser.add_argument("--max-parents", type=int, default=2, help="Maximum parents per node (default: 2).")
    parser.add_argument("--delimiter", type=str, default=",", help="CSV delimiter (default: ',').")
    parser.add_argument("--header", dest="header", action="store_true", help="CSV has header row (default).")
    parser.add_argument("--no-header", dest="header", action="store_false", help="CSV has no header row.")
    parser.set_defaults(header=True)
    parser.add_argument("--time-limit", type=int, default=None, help="Solver time limit in seconds (optional).")
    parser.add_argument("--score", type=str, default="bic", choices=["bic", "bdeu"],
                        help="Score function to use: bic or bdeu (default: bic).")
    parser.add_argument("--ess", type=float, default=1.0,
                        help="Equivalent sample size for BDeu (default: 1.0).")
    parser.add_argument("--out", type=str, default="./",
                        help="Output directory (default: ./).")

    args = parser.parse_args()

    csv_path = Path(args.csv_path)
    if not csv_path.exists():
        print(f"File not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    if args.header:
        df = pd.read_csv(csv_path, sep=args.delimiter)
    else:
        df = pd.read_csv(csv_path, sep=args.delimiter, header=None)
        df.columns = [f"X{i}" for i in range(df.shape[1])]

    # Validate binary data
    for col in df.columns:
        unique_vals = set(pd.Series(df[col]).dropna().unique().tolist())
        if not unique_vals.issubset({0, 1}):
            raise ValueError(f"Column '{col}' contains non-binary values: {sorted(unique_vals)}")

    sol = build_and_solve_ilp(df,
                              max_parents=args.max_parents,
                              time_limit=args.time_limit,
                              score=args.score,
                              ess=args.ess)

    # Save JSON next to input
    out_name = csv_path.with_suffix("").name
    out_json = args.out + out_name + "_bn_structure.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(sol, f, ensure_ascii=False, indent=2)

    # Pretty print
    print("Status:", sol["status"])
    if args.time_limit:
        print(f"Time limit (sec): {args.time_limit}")
    print("Score:", sol["score"], "(ESS:", sol["ess"], ")" if sol["score"] == "bdeu" else "")
    print("Objective (best-so-far total score):", sol["objective_value"])
    print("\nChosen parent sets:")
    for child, Ps in sol["selected_parent_sets"].items():
        print(f"  {child} <- {{{', '.join(Ps) if Ps else ''}}}")
    print("\nEdges:")
    for p, c in sol["edges"]:
        print(f"  {p} -> {c}")
    print(f"\nStructure JSON saved to: {out_json}")


if __name__ == "__main__":
    main()

