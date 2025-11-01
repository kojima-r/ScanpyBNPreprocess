mkdir -p ecv_all2

#ingor --bs prefix=bs_all2/result.ing,type=ing,ed=100,th=0.05 -o result_all2.sgn3
ingor --bs prefix=bs_all2/result.ing,type=ing,ed=1000,th=0.05 -o result_all2.sgn3
ingor --bs prefix=bs_all2/result.ing,type=ing,ed=1000,th=0.05 -o result_all2.txt

for f in `ls 02data_bbknn2_t/*`
do
name=`basename ${f}`
ingor --read file=result_all2.sgn3 --ec data="./02data_bbknn2_t/${name}",method=ECv -o "ecv_all2/${name}.txt"
#ingor --read file=result_all2.sgn3 --ec data="./03data_bbknn/all2.txt",method=ECv -o "ecv_all2/${name}.txt"
#sleep 1
done

#wait

#sleep 1

