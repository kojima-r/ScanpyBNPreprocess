#sh run_cmd.sh 03data_bbknn/all.txt

data=all2
path=03data_bbknn

mkdir -p bs_all2

ingor -B 1  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 2  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 3  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 4  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 5  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 6  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 7  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 8  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 9  -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
ingor -B 10 -N 100 --single-file off -o bs_${data}/result.ing ./${path}/${data}.txt &
sleep 1
wait
sleep 1



