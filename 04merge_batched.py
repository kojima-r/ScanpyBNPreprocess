import glob
import numpy as np
import os
import pandas as pd

def run(input_files,out_dir,out_filename):
    os.makedirs(out_dir,exist_ok=True)

    all_df=[]
    for filename in glob.glob(input_files):
        print(filename)
        df=pd.read_csv(filename,sep="\t",index_col=0)
        name, _ = os.path.splitext(os.path.basename(filename))
        df["tissue"]=name
        all_df.append(df)

    m_df=pd.concat(all_df,axis=0,join="inner")
    print(m_df)
    m_df.to_csv(out_dir+out_filename,sep="\t",index = True)

input_files="02data_bbknn_b/*.txt"
out_dir="03data_bbknn_b/"
out_filename="all.txt"
run(input_files,out_dir,out_filename)

