"""
Read a Tabula Muris Senis h5ad file and split it into per-tissue
tab-separated expression matrices (cells x highly variable genes).

Output: <out_dir>/<tissue>.txt with the first column being the cell
identifier "tissue|age|batch|cell_id" (the "batch" component is
omitted for FACS/droplet data which do not have a batch annotation).
"""

import argparse
import os

import scanpy as sc


DEFAULT_INPUTS = {
    "bbknn":   "data/tabula-muris-senis-bbknn-processed-official-annotations.h5ad",
    "facs":    "data/tabula-muris-senis-facs-processed-official-annotations.h5ad",
    "droplet": "data/tabula-muris-senis-droplet-processed-official-annotations.h5ad",
}


def write_tissue(adata, tissue, mode, out_path):
    sub = adata[adata.obs["tissue"] == tissue, adata.var["highly_variable"] == True]
    X = sub.X.todense()

    obs = sub.obs
    if mode == "bbknn":
        names = (obs["tissue"].astype(str) + "|" +
                 obs["age"].astype(str)    + "|" +
                 obs["batch"].astype(str)  + "|" +
                 obs.index.astype(str)).tolist()
    else:
        names = (obs["tissue"].astype(str) + "|" +
                 obs["age"].astype(str)    + "|" +
                 obs.index.astype(str)).tolist()

    genes = sub.var.index.tolist()
    with open(out_path, "w") as fp:
        fp.write("@name\t" + "\t".join(genes) + "\n")
        for i, name in enumerate(names):
            row = X[i].tolist()[0]
            fp.write(name + "\t" + "\t".join(map(str, row)) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=DEFAULT_INPUTS.keys(), default="bbknn",
                        help="Which Tabula Muris Senis dataset to process")
    parser.add_argument("--input", default=None,
                        help="Path to the input .h5ad file (defaults depend on --mode)")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory (default: data01_<mode>/)")
    args = parser.parse_args()

    input_path = args.input or DEFAULT_INPUTS[args.mode]
    out_dir = args.out_dir or f"data01_{args.mode}/"
    os.makedirs(out_dir, exist_ok=True)

    adata = sc.read(input_path)
    tissues = adata.obs["tissue"].unique().tolist()
    print(f"#tissues = {len(tissues)}")
    for tissue in tissues:
        out_path = os.path.join(out_dir, tissue + ".txt")
        print(f">> {out_path}")
        write_tissue(adata, tissue, args.mode, out_path)


if __name__ == "__main__":
    main()
