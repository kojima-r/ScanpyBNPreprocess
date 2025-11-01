import glob
import numpy as np
import os

N=10
out_dir="02data_bbknn/"
input_files="01data_bbknn/*.txt"
#Sample the same number of times as the original sample size.
M="same" #resample 
#Sample only the root of the original sample size.
# M="root"

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
        if M=="root":
            m=int(np.sqrt(len(id_list)))
        else: # M=="same"
            m=len(id_list)
        new_data=[]
        for i in range(N):
            idx=np.random.randint(0,len(id_list),size=m)
            X=data[idx,:]
            v=np.mean(X,axis=0)
            s="\t".join(map(str,v))
            ofp.write("s"+str(i)+"\t"+s)
            ofp.write("\n")

from multiprocessing import Pool

f_list=[filename for filename in glob.glob(input_files)]
p = Pool(16)
p.map(run,f_list)

