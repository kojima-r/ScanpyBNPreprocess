import glob
import numpy as np
import os

out_dir="02data_bbknn2/"
input_files="01data_bbknn/*.txt"

os.makedirs(out_dir,exist_ok=True)

def run(filename):
    print(filename)
    id_list=[]
    with open(filename) as fp:
        head=next(fp)
        data=[]
        for line in fp:
            arr=line.split("\t")
            id_list.append(arr[0])
            v=np.array(list(map(float,arr[1:])))
            data.append(v)
        data=np.array(data)
    
    name=os.path.basename(filename)
    out_name=out_dir+name
    print(">>",out_name)
    with open(out_name,"w") as ofp:
        ofp.write(head)
        new_data=[]
        v=np.mean(data,axis=0)
        s="\t".join(map(str,v))
        ofp.write("s"+str(name)+"\t"+s)
        ofp.write("\n")

from multiprocessing import Pool

f_list=[filename for filename in glob.glob(input_files)]
p = Pool(16)
p.map(run,f_list)

