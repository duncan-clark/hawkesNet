#!/bin/sh
#SBATCH --job-name=sim_study_CS
#SBATCH --output=cluster_output/CS_sim_%j.out
#SBATCH --error=cluster_output/CS_sim_%j.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=200G

# NeSI (Mahuika) notes:
# - Do not set --partition: let Slurm choose (default "large").
# - "milan" is opt-in only; requesting it can cause "configuration not available".
# - large: 72 CPUs/node, 105 GB/node, 3 days max → 128 CPUs = 2 nodes, 210 GB max; we request 200G.
# - For more memory (e.g. 256G) use: --partition=bigmem (if available), same --mem and --cpus-per-task.
# - Very large CPU/RAM requests often wait longer in queue.
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

# Use 120 workers for simulation phase (leaving buffer for R connection limit)
export SIM_WORKERS=120

Rscript inst/simulation_study/simulation_study_CS.R 2>&1
echo "Finished job at $(date)"
