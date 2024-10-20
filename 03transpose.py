import glob
import numpy as np
import os
import pandas as pd


def run(out_dir,input_files):
    os.makedirs(out_dir,exist_ok=True)

    for filename in glob.glob(input_files):
        print(filename)
        name=os.path.basename(filename)
        df=pd.read_csv(filename,sep="\t").transpose()
        df.to_csv(out_dir+name,sep="\t",header = False)


out_dir="02data_bbknn_t/"
input_files="02data_bbknn/*.txt"
run(out_dir,input_files):

out_dir="02data_bbknn2_t/"
input_files="02data_bbknn2/*.txt"
run(out_dir,input_files):

