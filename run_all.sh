#sh run_cmd.sh ./data/inc_cortex_tpm.txt
#sh run_cmd.sh ./data/inc_cortex_tpm_ctrl.txt
#sh run_cmd.sh ./data/inh_cortex_tpm.txt
#sh run_cmd.sh ./data/inh_cortex_tpm_ctrl.txt
#sh run_cmd.sh ./data/inc_hippo_tpm.txt
#sh run_cmd.sh ./data/inc_hippo_tpm_ctrl.txt

#sh run_cmd.sh ./data/inc_cortex_logtpm.large.txt
#sh run_cmd.sh ./data/inc_cortex_logtpm_ctrl.large.txt
#sh run_cmd.sh ./data/inh_cortex_logtpm.large.txt
#sh run_cmd.sh ./data/inh_cortex_logtpm_ctrl.large.txt
#sh run_cmd.sh ./data/inc_hippo_logtpm.large.txt
#sh run_cmd.sh ./data/inc_hippo_logtpm_ctrl.large.txt

#sh run_cmd.sh ./data/full_logtpm.large.txt
for f in `ls 02data_bbknn_t/*.txt`
do
	echo $f
	sh run_cmd.sh $f
done
