#!/bin/sh
#SBATCH --job-name=openalex_hawkes
#SBATCH --output=openalex_%j.out
#SBATCH --error=openalex_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# OpenAlex Hawkes study: inhomogeneous fit, temporal fit, KS tests, GOF.
# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   sbatch inst/openalex_study/run_openalex.slurm
#
# Optional env (set before sbatch or in script):
#   OPENALEX_EMAIL=your@email
#   OPENALEX_STRING=Hawkes
#   OPENALEX_PAGES=50

module load R

echo "Starting job at $(date)"
Rscript inst/openalex_study/openalex_hawkes_study.R
echo "Finished job at $(date)"
