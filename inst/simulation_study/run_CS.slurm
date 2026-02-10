#!/bin/sh
#SBATCH --job-name=sim_study_CS
#SBATCH --output=cluster_output/CS_sim_%j.out
#SBATCH --error=cluster_output/CS_sim_%j.err
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128G

# 128 cpus: script uses 26 outer x 4 inner = 104 (no oversubscription). If squeue shows 256, set USE_256_WHEN_128=1 to use 51x5=255.
# Override: CORES_OVERRIDE=128|256, USE_256_WHEN_128=1 (128 request, 256 allocated), CORES_OUTER, CORES_INNER.
#
# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   mkdir -p cluster_output && sbatch inst/simulation_study/run_CS.slurm

module load R

echo "Starting job at $(date)"
Rscript inst/simulation_study/simulation_study_CS.R
echo "Finished job at $(date)"
