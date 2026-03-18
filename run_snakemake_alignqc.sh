#!/bin/bash

#SBATCH --job-name=snakemake-alignqc
#SBATCH --mail-type=END
#SBATCH --mail-user=giuseppe.gautieri@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x_%A.log

# Load environment
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate alignqc_env
conda info | grep "active environment"

# This command removes the stale lock file.
snakemake --snakefile alignqc.smk --configfile config.yml --unlock


# Run Snakemake with SLURM
snakemake \
     --snakefile alignqc.smk \
     --configfile config.yml \
     --cores 4 \
     --use-singularity \
     --cluster "sbatch --job-name={rule}.{wildcards.sample} --partition=cpuq --time={resources.time} --cpus-per-task={threads} --mem={resources.mem_mb}M --output=logs/{rule}_%j.log" \
     --jobs 12 \
     --keep-going \
     --rerun-incomplete \
     --latency-wait 60