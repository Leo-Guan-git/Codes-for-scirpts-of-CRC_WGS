date
## !! dir = config.yaml[outdir]
dir=`readlink -f ./RT.SpecSim.result`
nohup  /software/miniconda3/miniconda3/envs/snakemake_env/bin/snakemake   --latency-wait 600   -j 999   -s   /RT.SpecSim.prediction.pepline.smk   --configfile   /RT.SpecSim.prediction.config.yaml    --cluster-config   /RT.SpecSim.prediction.cluster.yaml    --cluster "/opt/gridengine/bin/lx-amd64/qsub -S /bin/bash -wd $dir/{cluster.dir} -q {cluster.Q} -binding linear:{cluster.BL} -l {cluster.L} -P {cluster.P} -N {cluster.N}" --ri > ./nohup.out 2>&1 &

#  /software/miniconda3/miniconda3/envs/snakemake_env/bin/snakemake -s   /RT.SpecSim.prediction.pepline.smk --configfile   /RT.SpecSim.prediction.config.yaml -n -p
date
