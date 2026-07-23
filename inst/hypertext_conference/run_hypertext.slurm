#!/bin/bash
#SBATCH --job-name=hypertext_hawkes
#SBATCH --output=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/logs/hypertext_%j.out
#SBATCH --error=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/logs/hypertext_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=128G

# Hypertext 2009 conference: hawkesNet publication fit.
# Submit from package root:
#   source inst/cluster_env.sh
#   sbatch inst/hypertext_conference/run_hypertext.slurm

set -euo pipefail
PACKAGE_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [ ! -f "$PACKAGE_ROOT/DESCRIPTION" ]; then
  echo "ERROR: DESCRIPTION file not found. Run from package root."
  exit 1
fi
cd "$PACKAGE_ROOT"
# shellcheck disable=SC1091
source inst/cluster_env.sh
mkdir -p "${HAWKESNET_OUTPUT_DIR}/logs"

module load R

echo "=========================================="
echo "Hypertext Conference Hawkes Study"
echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "Node:      $(hostname)"
echo "CPUs:      ${SLURM_CPUS_PER_TASK:-}"
echo "HAWKESNET_OUTPUT_DIR=${HAWKESNET_OUTPUT_DIR}"
echo "Started:   $(date)"
echo "=========================================="

export N_GOF=100
export GOF_CORES_OUTER=100
export K_FIXED=1
export RUN_FIT_A=TRUE
export RUN_GOF_A=TRUE
export RUN_FIT_B=TRUE
export RUN_GOF_B=TRUE
export RUN_FIT_C=TRUE
export RUN_GOF_C=TRUE
export RUN_FIT_D=TRUE
export RUN_GOF_D=TRUE
Rscript inst/hypertext_conference/hypertext_fits.R 2>&1

echo "=========================================="
echo "Finished:  $(date)"
echo "=========================================="
