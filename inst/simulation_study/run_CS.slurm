#!/bin/bash
#SBATCH --job-name=sim_study_CS
#SBATCH --output=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/logs/CS_sim_%j.out
#SBATCH --error=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/logs/CS_sim_%j.err
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G

# Submit from the package root:
#   source inst/cluster_env.sh
#   sbatch inst/simulation_study/run_CS.slurm

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
# shellcheck disable=SC1091
source inst/cluster_env.sh
mkdir -p "${HAWKESNET_OUTPUT_DIR}/logs"

module load R

echo "Starting job at $(date)"
echo "HAWKESNET_OUTPUT_DIR=${HAWKESNET_OUTPUT_DIR}"
Rscript inst/simulation_study/simulation_study_CS.R 2>&1
echo "Finished job at $(date)"
