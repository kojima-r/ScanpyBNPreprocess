import glob
import numpy as np
import os

N=10
out_dir="02data_bbknn_b/"
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
            id_arr=arr[0].split("|")
            id_list.append(id_arr[1]+"|"+id_arr[2])
            v=np.array(list(map(float,arr[1:])))
            data.append(v)
        data=np.array(data)
    u_id_list=np.unique(id_list)
    id_list=np.array(id_list)
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
            for target_id in u_id_list:
                target_id_list=np.where(id_list==target_id)[0]
                idx=np.random.choice(target_id_list,size=m,replace=True)
                X=data[idx,:]
                v=np.mean(X,axis=0)
                s="\t".join(map(str,v))
                ofp.write("s"+str(i)+"|"+target_id+"\t"+s)
                ofp.write("\n")

from multiprocessing import Pool

f_list=[filename for filename in glob.glob(input_files)]
p = Pool(16)
p.map(run,f_list)

