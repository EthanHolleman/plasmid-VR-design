mkdir -p logs
source ~/anaconda3/etc/profile.d/conda.sh
conda activate snakemake
snakemake -j 24 --cluster-config cluster.yml \
--cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} \
-n {cluster.cpus} --mem {cluster.mem} -J {cluster.name} -o {cluster.output} \
-e {cluster.output} --mail-type ALL --mail-user {cluster.email}" \
--conda-frontend=mamba --configfile config.yml \
--latency-wait 60 --verbose --use-conda --rerun-incomplete