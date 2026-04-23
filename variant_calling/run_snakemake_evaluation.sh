#!/bin/bash

#SBATCH --job-name=snakemake_benchmark_variant_calling
#SBATCH --mail-type=END
#SBATCH --mail-user=giuseppe.gautieri@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=%x_%A.log

# 1. Load your base Snakemake environment
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate alignqc_env  

# 2. Create a logs directory if it doesn't exist
mkdir -p logs

# 3. Unlock the directory
snakemake --snakefile benchmark_master.smk --configfile config.yml --unlock

# 4. Run Snakemake

snakemake \
     --snakefile benchmark_master.smk \
     --configfile config.yml \
     --cores 50 \
     --jobs 20 \
     --use-conda \
     --conda-frontend mamba \
     --use-singularity \
     --keep-going \
     --rerun-incomplete \
     --latency-wait 120 \
     --cluster "sbatch \
               --job-name=smk.{rule} \
                --partition=cpuq \
                --time={resources.time} \
                --cpus-per-task={threads} \
                --mem={resources.mem_mb}M \
                --output=logs/{rule}_%j.log" \
     