mkdir -p ecv_all

#ingor --bs prefix=bs_all/result.ing,type=ing,ed=100,th=0.05 -o result_all.sgn3
ingor --bs prefix=bs_all/result.ing,type=ing,ed=1000,th=0.05 -o result_all.sgn3
ingor --bs prefix=bs_all/result.ing,type=ing,ed=1000,th=0.05 -o result_all.txt

for f in `ls 02data_bbknn_t/*`
do
name=`basename ${f}`
echo ${name}
ingor --read file=./result_all.sgn3 --ec data="./02data_bbknn_t/${name}",method=ECv -o "ecv_all/${name}.txt"
#ingor --read file=result_all.sgn3 --ec data="./03data_bbknn/all.txt",method=ECv -o "ecv_all/${name}.txt"
#sleep 1
done

#wait

#sleep 1

