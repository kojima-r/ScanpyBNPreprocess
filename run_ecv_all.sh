mkdir -p ecv_all_

for f in `ls data02_bbknn/*`
do
b=`basename ${f}`
ingor --read file=result_all_.sgn3 --ec data=./data02_bbknn/${name},method=ECv -o ecv_all_/${name}.txt
done
