import glob
import numpy as np
import os
import pandas as pd
os.makedirs("02data_bbknn_t/",exist_ok=True)

for filename in glob.glob("02data_bbknn/*.txt"):
    print(filename)
    name=os.path.basename(filename)
    df=pd.read_csv(filename,sep="\t").transpose()
    df.to_csv("02data_bbknn_t/"+name,sep="\t",header = False)



