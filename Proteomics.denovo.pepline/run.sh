date
# rm -rf .snakemake/
logdir=`readlink -f ./log`
# rm -rf $logdir
mkdir -p $logdir
nohup  /software/miniconda3/miniconda3/envs/snakemake_env/bin/snakemake --latency-wait 600 -j 999 -s DENOVO.pepline.smk --configfile DENOVO.config.yaml    --cluster-config DENOVO.cluster.yaml    --cluster "/opt/gridengine/bin/lx-amd64/qsub -S /bin/bash -wd $logdir -q {cluster.Q} -binding linear:{cluster.BL} -l {cluster.L} -P {cluster.P}" --ri  > ./nohup.out 2>&1 &
#  /software/miniconda3/miniconda3/envs/snakemake_env/bin/snakemake -s DENOVO.pepline.smk --configfile DENOVO.pepline.smk -n -p
date
