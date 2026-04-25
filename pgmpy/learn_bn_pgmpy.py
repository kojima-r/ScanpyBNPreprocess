import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pgmpy.estimators import HillClimbSearch, PC, BDeu, BIC, K2, ExpertKnowledge
#from pgmpy.models import BayesianNetwork
from pgmpy.models import DiscreteBayesianNetwork
from pgmpy.estimators import BayesianEstimator
from pgmpy.readwrite import BIFWriter

# plotting
#import networkx as nx
#import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Structure learning for a Bayesian Network from a CSV using pgmpy (1.0.0)."
    )
    p.add_argument("csv", type=str, help="Path to input CSV file.")
    p.add_argument("--id-column", type=str, default=None,
                   help="Optional column to drop (e.g., sample ID).")
    p.add_argument("--na-policy", choices=["drop-rows", "fill-mean", "fill-median", "fill-mode"],
                   default="drop-rows", help="How to handle missing values.")
    p.add_argument("--discretize", choices=["none", "quantile", "uniform"], default="quantile",
                   help="Discretization method for continuous variables.")
    p.add_argument("--bins", type=int, default=3, help="Number of bins for discretization.")
    p.add_argument("--na-as-bin", action="store_true",
                   help="When discretizing (quantile/uniform), put NaNs into a dedicated extra bin (code: last bin index).")
    p.add_argument("--var-threshold", type=float, default=0.0,
                   help="Drop columns whose variance is <= this threshold (applied before discretization).")
    p.add_argument("--estimator", choices=["hc", "pc" ,"hybrid"], default="hc",
                   help="Structure learning algorithm: 'hc' (HillClimb) or 'pc' (constraint-based).")
    p.add_argument("--score", choices=["k2", "bdeu", "bic"], default="k2",
                   help="Scoring method for hill climb.")
    p.add_argument("--equivalent-sample-size", type=float, default=10.0,
                   help="Equivalent sample size for BDeu score and Bayesian parameter fitting.")
    p.add_argument("--max-indegree", type=int, default=5,
                   help="Optional max number of parents per node (HillClimb only).")
    p.add_argument("--max-iter", type=int, default=1000,
                   help="Optional max number of iterations (HillClimb only).")
    p.add_argument("--whitelist", type=str, default=None,
                   help="Optional path to a JSON file with list of required edges [['A','B'], ...].")
    p.add_argument("--blacklist", type=str, default=None,
                   help="Optional path to a JSON file with list of forbidden edges [['A','B'], ...].")
    p.add_argument("--temporal-order", type=str, default=None,
                   help="Optional path to JSON list specifying temporal/causal ordering of variables ['T0','T1', ...].")
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    p.add_argument("--output-prefix", type=str, default="bn_result",
                   help="Prefix for output files (PNG, JSON, BIF).")
    p.add_argument("--significance-level", type=float, default=0.01,
                   help="Significance level for PC algorithm (constraint-based).")
    p.add_argument("--enforce-expert", action="store_true",
                   help="If set, algorithms will try to enforce expert edges (where supported).")
    p.add_argument("--n-jobs", type=int, default=8,
                   help="Number of cpu cores to perform pc algorithm (constraint-based).")
    return p.parse_args()


def handle_missing(df: pd.DataFrame, policy: str) -> pd.DataFrame:
    if policy == "drop-rows":
        return df.dropna(axis=0)
    if policy == "fill-mean":
        return df.fillna(df.mean(numeric_only=True))
    if policy == "fill-median":
        return df.fillna(df.median(numeric_only=True))
    if policy == "fill-mode":
        modes = df.mode(dropna=True)
        if not modes.empty:
            mode_vals = modes.iloc[0]
            return df.fillna(mode_vals)
        return df.fillna(0)
    raise ValueError(f"Unknown na-policy: {policy}")


def discretize(df: pd.DataFrame, method: str, bins: int, na_as_bin: bool) -> pd.DataFrame:
    result = df.copy()

    if method == "none":
        for col in result.select_dtypes(include=["object", "category"]).columns:
            codes = result[col].astype("category").cat.codes  # NaN -> -1
            if na_as_bin:
                if (codes == -1).any():
                    max_code = int(codes[codes >= 0].max()) if (codes >= 0).any() else -1
                    codes = codes.replace(-1, max_code + 1)
            result[col] = codes
        return result.astype("int64", errors="ignore")

    numeric_cols = result.select_dtypes(include=[np.number]).columns.tolist()
    for col in numeric_cols:
        series = result[col]
        if series.nunique(dropna=True) <= 1:
            if na_as_bin and series.isna().any():
                result[col] = np.where(series.isna(), 1, 0)
            else:
                result[col] = 0
            continue

        if method == "quantile":
            try:
                binned = pd.qcut(series, q=bins, labels=False, duplicates="drop")
            except ValueError:
                binned = pd.cut(series, bins=bins, labels=False, include_lowest=True)
        elif method == "uniform":
            binned = pd.cut(series, bins=bins, labels=False, include_lowest=True)
        else:
            raise ValueError(f"Unknown discretize method: {method}")

        binned = binned.astype("float64")
        if na_as_bin:
            if np.isnan(binned).all():
                binned = np.where(series.isna(), 0, np.nan)
            max_code = int(np.nanmax(binned)) if not np.isnan(binned).all() else -1
            na_code = max_code + 1
            binned = np.where(series.isna(), na_code, binned)
        else:
            if np.isnan(binned).any():
                vals, counts = np.unique(binned[~np.isnan(binned)], return_counts=True)
                if len(vals) > 0:
                    mode_val = vals[np.argmax(counts)]
                    binned = np.where(np.isnan(binned), mode_val, binned)
                else:
                    binned = np.zeros_like(binned)

        result[col] = binned.astype(int)

    for col in result.select_dtypes(include=["object", "category"]).columns:
        codes = result[col].astype("category").cat.codes  # NaN -> -1
        if na_as_bin and (codes == -1).any():
            max_code = int(codes[codes >= 0].max()) if (codes >= 0).any() else -1
            codes = codes.replace(-1, max_code + 1)
        else:
            if (codes == -1).any():
                nonneg = codes[codes >= 0]
                fill_code = int(nonneg.mode().iloc[0]) if not nonneg.empty else 0
                codes = codes.replace(-1, fill_code)
        result[col] = codes

    return result.astype("int64", errors="ignore")


def load_edge_list(path):
    if path is None:
        return None
    with open(path, "r", encoding="utf-8") as f:
        return [tuple(e) for e in json.load(f)]


def load_temporal_order(path):
    if path is None:
        return None
    with open(path, "r", encoding="utf-8") as f:
        arr = json.load(f)
        if not isinstance(arr, list):
            raise ValueError("--temporal-order must be a JSON list of variable names")
        return arr


def make_expert_knowledge(whitelist, blacklist, temporal):
    ek = ExpertKnowledge()
    if whitelist:
        ek.required_edges = list(map(tuple, whitelist))
    if blacklist:
        ek.forbidden_edges = list(map(tuple, blacklist))
    if temporal:
        ek.temporal_order = list(temporal)
    return ek


def get_score(df: pd.DataFrame, score_name: str, ess: float):
    s = score_name.lower()
    if s == "k2":
        return K2(df)
    if s == "bdeu":
        return BDeu(df, equivalent_sample_size=ess)
    if s == "bic":
        return BIC(df)
    raise ValueError(f"Unknown score: {score_name}")


def learn_structure_hc(df: pd.DataFrame, score_name: str, ess: float,max_iter, max_indegree, expert, seed: int, start_dag=None):
    print("start: hc algorithm")
    score = get_score(df, score_name, ess)
    search = HillClimbSearch(df,use_cache=True)
    dag = search.estimate(
            max_indegree=max_indegree,
            expert_knowledge=expert,
            max_iter=max_iter,
            scoring_method=score,
            start_dag=start_dag,
            )
    return dag


def learn_structure_pc(df: pd.DataFrame, significance_level: float, expert, enforce: bool):
    print("start: pc algorithm")
    pc = PC(data=df)
    # PC supports enforce_expert_knowledge flag (docs 1.0.0)
    dag = pc.estimate(
            significance_level=significance_level,
            expert_knowledge=expert,
            enforce_expert_knowledge=enforce,
            return_type="dag",
            n_jobs=32)
    return dag


def fit_parameters(model: DiscreteBayesianNetwork, df: pd.DataFrame, ess: float):
    model.fit(df, estimator=BayesianEstimator, prior_type="BDeu", equivalent_sample_size=ess)
    return model


def save_outputs(model: DiscreteBayesianNetwork, prefix: str):
    edges = list(model.edges())
    nodes = list(model.nodes())
    print("#edges:",len(edges))
    print("#nodes:",len(nodes))
    with open(f"{prefix}_edges.json", "w", encoding="utf-8") as f:
        json.dump(edges, f, ensure_ascii=False, indent=2)

    #bif_path = f"{prefix}_model.bif"
    #writer = BIFWriter(model)
    #writer.write_bif(filename=bif_path)

    """
    plt.figure(figsize=(max(6, len(nodes) * 0.5), max(6, len(nodes) * 0.5)))
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    pos = nx.spring_layout(G, seed=42)
    nx.draw_networkx_nodes(G, pos)
    nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle='-|>', arrowsize=12)
    nx.draw_networkx_labels(G, pos, font_size=8)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"{prefix}_graph.png", dpi=200)
    plt.close()
    """
    return {
        "edges_json": f"{prefix}_edges.json",
        #"bif": f"{prefix}_model.bif",
        #"graph_png": f"{prefix}_graph.png",
    }


def main():
    args = parse_args()
    np.random.seed(args.seed)

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"CSV not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    df_disc = pd.read_csv(csv_path, sep="\t")
    whitelist = load_edge_list(args.whitelist)
    blacklist = load_edge_list(args.blacklist)
    temporal = load_temporal_order(args.temporal_order)
    expert = make_expert_knowledge(whitelist, blacklist, temporal)
    
    print("input dataset(col x row): ",len(df_disc.columns),"x",len(df_disc))
    if args.estimator == "hc":
        dag = learn_structure_hc(
            df_disc, args.score, args.equivalent_sample_size, args.max_iter,  args.max_indegree, expert, args.seed,
            )
    elif args.estimator == "hybrid":
        dag = learn_structure_pc(df_disc, args.significance_level, expert, args.enforce_expert)
        dag = learn_structure_hc(
            df_disc, args.score, args.equivalent_sample_size, args.max_iter,  args.max_indegree, expert, args.seed,
            start_dag=dag
        )
    else:
        dag = learn_structure_pc(df_disc, args.significance_level, expert, args.enforce_expert)

    # Convert to BayesianNetwork and fit CPDs
    #model = BayesianNetwork(dag.edges())
    model = DiscreteBayesianNetwork(dag.edges())
    model.add_nodes_from(dag.nodes())
    model = fit_parameters(model, df_disc, args.equivalent_sample_size)

    outputs = save_outputs(model, args.output_prefix)

    print("=== Learned structure (edges) ===")
    for u, v in model.edges():
        print(f"{u} -> {v}")
    print("\nSaved files:")
    for k, v in outputs.items():
        print(f"- {k}: {v}")


if __name__ == "__main__":
    main()

