import glob
import numpy as np
import os
import pandas as pd
N=10
os.makedirs("03data_bbknn/",exist_ok=True)

all_df=[]
for filename in glob.glob("02data_bbknn_t/*.txt"):
    print(filename)
    name=os.path.basename(filename)
    df=pd.read_csv(filename,sep="\t",index_col=0)
    print(df.columns)
    cs=[name+"_"+e for e in df.columns]
    df.columns=cs
    print(df)
    all_df.append(df)


m_df=pd.concat(all_df,axis=1,join="inner")
print(m_df)
m_df.to_csv("03data_bbknn/all.txt",sep="\t",index = True)
