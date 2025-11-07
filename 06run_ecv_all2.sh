find ./bs_all2 -type f -size -1024c | xargs rm
mkdir -p ecv_all2

#ingor --bs prefix=bs_all2/result.ing,type=ing,ed=100,th=0.05 -o result_all2.sgn3
ingor --bs prefix=bs_all2/result.ing,type=ing,ed=1000,th=0.05 -o result_all2.sgn3
ingor --bs prefix=bs_all2/result.ing,type=ing,ed=1000,th=0.05 -o result_all2.txt

for f in `ls 02data_bbknn2_t/*.txt`
do
name=`basename ${f} .txt`
ingor --read file=result_all2.sgn3 --ec data="./02data_bbknn2_t/${name}.txt",method=ECv -o "ecv_all2/${name}.txt"
done



