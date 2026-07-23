#!/bin/bash
#SBATCH --job-name=openalex_hawkes
#SBATCH --output=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/logs/openalex_%j.out
#SBATCH --error=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/logs/openalex_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=128G

# Submit from package root:
#   source inst/cluster_env.sh
#   sbatch inst/openalex_study/run_openalex.slurm

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
echo "OpenAlex Hawkes Study"
echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "HAWKESNET_OUTPUT_DIR=${HAWKESNET_OUTPUT_DIR}"
echo "Started:   $(date)"
echo "=========================================="

Rscript inst/openalex_study/openalex_hawkes_study.R 2>&1

echo "=========================================="
echo "Finished:  $(date)"
echo "=========================================="
