#!/bin/bash

#SBATCH --job-name=StructuralVarAnnotation
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=5-00:00:00
#SBATCH --output=%x_%A.log
#SBATCH --error=logs/StructuralVarAnnotation_%j.err

set -euo pipefail

cd "${SVA_REPO_DIR}"

mkdir -p logs
mkdir -p "${SVA_OUTDIR}/logs"

# ------------------------------------------------------------
# Optional private conda/snakemake environment
# ------------------------------------------------------------

CONDA_BIN="${SVA_CONDA_BIN:-}"

if [[ -n "$CONDA_BIN" ]]; then
    CONDA_BIN="${CONDA_BIN%/}"

    if [[ ! -x "${CONDA_BIN}/conda" ]]; then
        echo "ERROR: SVA_CONDA_BIN does not contain executable conda: ${CONDA_BIN}/conda" >&2
        exit 1
    fi

    export PATH="${CONDA_BIN}:$PATH"

    # Avoid cluster Python/Spack pollution
    unset PYTHONPATH
    unset PYTHONHOME

    # Use conda libmamba solver
    unset CONDA_NO_PLUGINS
    export CONDA_SOLVER=libmamba
else
    if command -v module >/dev/null 2>&1; then
        module load "${SVA_SNAKEMAKE_MODULE:-snakemake/9.6.3-python-3.11.5}"
    elif [[ -f /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
        module load "${SVA_SNAKEMAKE_MODULE:-snakemake/9.6.3-python-3.11.5}"
    else
        echo "ERROR: module command not available and SVA_CONDA_BIN was not provided." >&2
        exit 1
    fi
fi

echo "Using Snakemake:"
which snakemake
snakemake --version

echo "Using conda:"
which conda || true
conda --version || true
conda config --show solver || true

CMD=(
    bin/StructuralVarAnnotation
    --samples "${SVA_SAMPLES}"
    --config "${SVA_CONFIG}"
    --outdir "${SVA_OUTDIR}"
    --jobs "${SVA_JOBS}"
    --partition "${SVA_PARTITION}"
    --cluster-jobs
)

if [[ -n "${SVA_TARGETS:-}" ]]; then
    read -r -a TARGET_ARRAY <<< "${SVA_TARGETS}"
    for target in "${TARGET_ARRAY[@]}"; do
        CMD+=(--target "$target")
    done
fi

if [[ -n "${SVA_CONDA_BIN:-}" ]]; then
    CMD+=(--conda-bin "${SVA_CONDA_BIN}")
fi

if [[ "${SVA_DRYRUN:-0}" -eq 1 ]]; then
    CMD+=(--dry-run)
fi

echo "Running inside SLURM controller:"
printf ' %q' "${CMD[@]}"
echo
echo

"${CMD[@]}"