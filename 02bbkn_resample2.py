import glob
import numpy as np
import os

os.makedirs("02data_bbknn2/",exist_ok=True)


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
    out_name="02data_bbknn2/"+name
    print(">>",out_name)
    with open(out_name,"w") as ofp:
        ofp.write(head)
        #m=int(np.sqrt(len(id_list)))
        m=len(id_list)
        new_data=[]
        v=np.mean(data,axis=0)
        s="\t".join(map(str,v))
        ofp.write("s"+str(name)+"\t"+s)
        ofp.write("\n")
        #new_data.append(v)
        #print(np.array(new_data).shape)

from multiprocessing import Pool

f_list=[filename for filename in glob.glob("01data_bbknn/*.txt")]
p = Pool(16)
p.map(run,f_list)

