data_path=$1
path=02data_bbknn
data=`basename ${data_path} .txt`
mkdir -p bs_${data}


ingor -B 1  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 2  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 3  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 4  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 5  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 6  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 7  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 8  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 9  -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 10 -N 10 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
sleep 1
wait
sleep 1


ingor --bs prefix=bs_${data}/result.ing,type=ing,ed=1000,th=0.05 -o result_${data}.txt

#ingor --read file=result.txt --score data=sample002-p10-n20.txt -o result_param.txt

#ingor --bs prefix=bs/result.ing,type=ing,ed=100,th=0.05 -o result.sgn3
#ingor --read file=result.sgn3 --ec data=sample002-p10-n20.txt,method=ECv -o result_ecv.txt
