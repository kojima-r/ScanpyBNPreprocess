#sh run_cmd.sh 03data_bbknn/all.txt

data=all
path=03data_bbknn

mkdir -p bs_all

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



