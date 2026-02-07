#!/bin/sh
#SBATCH --job-name=openalex_hawkes
#SBATCH --output=cluster_output/openalex_%j.out
#SBATCH --error=cluster_output/openalex_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# OpenAlex Hawkes study: inhomogeneous fit, temporal fit, KS tests, GOF.
# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   mkdir -p cluster_output && sbatch inst/openalex_study/run_openalex.slurm
#
# Optional env (set before sbatch or in script):
#   OPENALEX_EMAIL=your@email
#   OPENALEX_STRING=Hawkes
#   OPENALEX_PAGES=50

mkdir -p cluster_output

module load R

echo "=========================================="
echo "OpenAlex Hawkes Study"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $(hostname)"
echo "CPUs:      $SLURM_CPUS_PER_TASK"
echo "Memory:    $SLURM_MEM_PER_NODE MB"
echo "Started:   $(date)"
echo "=========================================="

Rscript inst/openalex_study/openalex_hawkes_study.R

echo "=========================================="
echo "Finished:  $(date)"
echo "=========================================="
