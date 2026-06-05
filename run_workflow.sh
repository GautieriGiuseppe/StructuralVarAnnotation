#!/bin/bash
#SBATCH --job-name=StructuralVarAnnotation
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/StructuralVarAnnotation_%j.out
#SBATCH --error=logs/StructuralVarAnnotation_%j.err

set -euo pipefail

cd "${SVA_REPO_DIR}"

mkdir -p logs
mkdir -p "${SVA_OUTDIR}/logs"

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate alignqc_env

bin/StructuralVarAnnotation \
    --samples "${SVA_SAMPLES}" \
    --config "${SVA_CONFIG}" \
    --outdir "${SVA_OUTDIR}" \
    --jobs "${SVA_JOBS}" \
    --partition "${SVA_PARTITION}" \
    --cluster-jobs