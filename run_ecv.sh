data_path=$1
target_path=$2
data=`basename ${data_path} .txt`

target=`basename ${target_path} .txt`

ingor --bs prefix=bs_${data}/result.ing,type=ing,ed=100,th=0.05 -o result_${data}.sgn3
ingor --read file=result_${data}.sgn3 --ec data=${target_path},method=ECv -o ecv_${data}_${target}.txt


