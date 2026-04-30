"""
Read a Tabula Muris Senis h5ad file and split it into per-tissue
tab-separated expression matrices (cells x highly variable genes).

Output: <out_dir>/<tissue>.txt with the first column being the cell
identifier "tissue|age|batch|cell_id" (the "batch" component is
omitted for FACS/droplet data which do not have a batch annotation).
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

    parser.add_argument("--mode", choices=["10X","smartseq"], default="10X",
                        help="Which Tabula Muris Senis dataset to process")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory (default: data01_tps/)")
    args = parser.parse_args()

    out_dir = args.out_dir or f"data01_tps/"
    os.makedirs(out_dir, exist_ok=True)
    
    for input_path in glob.glob("data/Tabula_Sapiens*.h5ad"):
        adata = sc.read(input_path)
        sub =adata[(adata.obs["method"]==args.mode),:]
        tissues = sub.obs["tissue"].unique().tolist()
        print(input_path)
        print(f"#tissues = {tissues}")
        for tissue in tissues:
            tissue_name=tissue.replace(" ","_")
            out_path = os.path.join(out_dir, tissue_name + "."+args.mode+".txt")
            print(f">> {out_path}")
            write_tissue(adata, args.mode, tissue, out_path)


if __name__ == "__main__":
    main()
