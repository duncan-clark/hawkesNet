#!/bin/sh
#SBATCH --job-name=sim_study_CS
#SBATCH --output=cluster_output/CS_sim_%j.out
#SBATCH --error=cluster_output/CS_sim_%j.err
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=32G

# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   mkdir -p cluster_output && sbatch inst/simulation_study/run_CS.slurm

module load R

echo "Starting job at $(date)"
Rscript inst/simulation_study/simulation_study_CS.R
echo "Finished job at $(date)"
