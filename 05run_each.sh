
for f in `ls 02data_bbknn_t/*.txt`
do
	echo $f
	sh run_cmd.sh $f
done
