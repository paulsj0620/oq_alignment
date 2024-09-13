wkdir="/data/SJ/playground/oq_align"
conffn="/data/SJ/playground/oq_align/dev.config.yaml"

mkdir -p ${wkdir}/logs
cmd="snakemake \
    --cores all \
    --printshellcmds \
    --snakefile /data/SJ/pipeline/onequeue_alignment/Snakefile \
    --config 'ref=/data/SJ/pipeline/onequeue_alignment/refs/ref.yaml' \
    --max-jobs-per-second 5 \
    --jobs 18 \
    --configfile ${conffn} \
    --directory ${wkdir} "
#    --cluster-config /data/SJ/pipeline/onequeue_alignment/envs/cluster.json \
#    --cluster 'sbatch /bin/bash'
echo $cmd
