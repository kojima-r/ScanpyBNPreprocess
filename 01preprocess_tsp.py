"""
Read Tabula Sapiens h5ad files and split them into per-tissue
tab-separated expression matrices (cells x genes).

Output: <out_dir>/<tissue>.txt with the first column being the cell
identifier "tissue|development_stage|cell_type|cell_id".
"""

import argparse
import os
import glob

import scanpy as sc

def write_tissue(adata, mode, tissue, out_path):
    #sub =adata[(adata.obs["method"]=="x10")&(adata.obs["tissue"] == tissue),adata.var["std"]>0.01]
    sub =adata[(adata.obs["method"]==mode)&(adata.obs["tissue"] == tissue),:]
    print(sub.shape)
    X = sub.X.todense()

    obs = sub.obs
    print(obs.keys())
    names = (obs["tissue"].astype(str) + "|" +
             obs["development_stage"].astype(str)    + "|" +
             obs["cell_type"].astype(str)  + "|" +
             obs.index.astype(str)).tolist()

    genes = sub.var.index.tolist()
    with open(out_path, "w") as fp:
        fp.write("@name\t" + "\t".join(genes) + "\n")
        for i, name in enumerate(names):
            row = X[i].tolist()[0]
            fp.write(name + "\t" + "\t".join(map(str, row)) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--mode", choices=["10x", "smartseq"], default="10x",
                        help="Which Tabula Sapiens assay to process")
    parser.add_argument("--input-glob", default="data/Tabula_Sapiens*.h5ad",
                        help="Glob pattern for input .h5ad files")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory (default: data01_ss/ for smartseq, data01_10x/ for 10x)")
    args = parser.parse_args()

    if args.mode == "smartseq":
        out_dir = args.out_dir or f"data01_ss/"
    else: # 10X
        out_dir = args.out_dir or f"data01_10x/"
    os.makedirs(out_dir, exist_ok=True)

    data_mode = "10X" if args.mode == "10x" else "smartseq"
    input_paths = sorted(glob.glob(args.input_glob))
    if not input_paths:
        parser.error(f"no input files matched: {args.input_glob}")
    for input_path in input_paths:
        adata = sc.read(input_path)
        sub =adata[(adata.obs["method"]==data_mode),:]
        tissues = sub.obs["tissue"].unique().tolist()
        print(input_path)
        print(f"#tissues = {tissues}")
        for tissue in tissues:
            tissue_name=tissue.replace(" ","_")
            out_path = os.path.join(out_dir, tissue_name + "." + args.mode + ".txt")
            print(f">> {out_path}")
            write_tissue(adata, data_mode, tissue, out_path)


if __name__ == "__main__":
    main()
