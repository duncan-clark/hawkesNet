#!/bin/bash
#SBATCH --job-name=sim_study_BA
#SBATCH --output=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/logs/BA_sim_%j.out
#SBATCH --error=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/logs/BA_sim_%j.err
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=32G

# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   source inst/cluster_env.sh   # optional; sets HAWKESNET_OUTPUT_DIR
#   sbatch inst/simulation_study/run_BA.slurm
#
# Override logs at submit time if needed:
#   sbatch --output=$HAWKESNET_OUTPUT_DIR/logs/BA_sim_%j.out --error=$HAWKESNET_OUTPUT_DIR/logs/BA_sim_%j.err ...

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
# shellcheck disable=SC1091
source inst/cluster_env.sh
mkdir -p "${HAWKESNET_OUTPUT_DIR}/logs"

module load R

echo "Starting job at $(date)"
echo "HAWKESNET_OUTPUT_DIR=${HAWKESNET_OUTPUT_DIR}"
Rscript inst/simulation_study/simulation_study_BA.R 2>&1
echo "Finished job at $(date)"
