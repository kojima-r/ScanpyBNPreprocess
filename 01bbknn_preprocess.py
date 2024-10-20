import numpy as np
import pandas as pd
import scanpy as sc
import os

results_file = 'tabula-muris-senis-bbknn-processed-official-annotations.h5ad'


adata = sc.read(results_file)
os.makedirs("01data_bbknn/",exist_ok=True)

print(len(adata.obs["tissue"].unique().tolist()))
for tissue in adata.obs["tissue"].unique().tolist():
    Z=adata[adata.obs["tissue"]==tissue, adata.var["highly_variable"]==True ]
    a=Z.X.todense()
    obs=(Z.obs["tissue"].astype(str)+"|"+Z.obs["age"].astype(str)+"|"+Z.obs.index.astype(str)).tolist()
    col=Z.var.index.tolist()
    with open("01data_bbknn/"+tissue+".txt","w") as fp:
        h="@name\t"+"\t".join(col)
        fp.write(h)
        fp.write("\n")
        for i, el in enumerate(obs):
            v=a[i].tolist()[0]
            line=el+"\t"+"\t".join(map(str,v))
            fp.write(line)
            fp.write("\n")

