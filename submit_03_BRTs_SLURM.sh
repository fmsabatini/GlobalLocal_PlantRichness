#!/bin/bash

#SBATCH --job-name=BRTs_direct99

#SBATCH --output=/work/%u/%x-%A-%a.log

#SBATCH --time 0-15:00:00
#SBATCH --mem-per-cpu=30G

#SBATCH --cpus-per-task=1

#SBATCH --array=1-99

module load foss/2019b R/4.0.0

fornonf="$1"
shift

output="/work/$USER/BRT_BiasCorr/$fornonf$SLURM_JOB_NAME-$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID"

Rscript \
    cli_03_BRTs.r \
    --verbose \
    --fornonf "$fornonf" \
    --nrows 99 \
    --index "$SLURM_ARRAY_TASK_ID" \
    /data/splot/_data/Mydata_global_NatCommR1.RData \
    /data/splot/_data/world.over.RData \
    "$output" 

