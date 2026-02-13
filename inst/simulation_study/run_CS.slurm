#!/bin/sh
#SBATCH --job-name=sim_study_CS
#SBATCH --output=cluster_output/CS_sim_%j.out
#SBATCH --error=cluster_output/CS_sim_%j.err
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G

# NeSI scheduling note:
# - Very large CPU/RAM requests often wait much longer in queue.
# - Defaulting to 64 CPUs / 64G tends to schedule faster than 128-core jobs; scale up via overrides if needed.
#
# Overrides:
#   CORES_OVERRIDE=<int>        (force N_CORES in the R script)
#   CORES_OUTER=<int>           (PSOCK outer workers for fitting)
#   CORES_INNER=<int>           (inner cores per fit for intensity cache)
#   SIM_WORKERS=<int>           (PSOCK workers for simulation phase)
#
# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   mkdir -p cluster_output && sbatch inst/simulation_study/run_CS.slurm

module load R

echo "Starting job at $(date)"
# Avoid nested BLAS/OpenMP threads stealing CPUs (and helps fork safety).
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript inst/simulation_study/simulation_study_CS.R
echo "Finished job at $(date)"
