import os
import glob
data={}
with open("data/GSE98969_experimental_design_f.txt") as fp:
    for _ in range(18):
        next(fp)
    for line in fp:
        arr=line.strip().split("\t")
        #print(arr)
        key=arr[2]
        batch_id=arr[1]
        label=arr[6]
        data[key]=(batch_id,label)
with open("info.tsv","w") as fp:
    for k,v in data.items():
        arr=[k,v[0],v[1]]
        s="\t".join(map(str,arr))
        fp.write(s)
        fp.write("\n")

quit()
for filename in glob.glob("data/GSM*.txt"):
    print(filename)
