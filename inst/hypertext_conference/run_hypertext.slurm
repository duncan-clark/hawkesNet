#!/bin/sh
#SBATCH --job-name=hypertext_hawkes
#SBATCH --output=cluster_output/hypertext_%j.out
#SBATCH --error=cluster_output/hypertext_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=128G

# Hypertext 2009 conference: hawkesNet fit + GOF (inhomogeneous background, gwdegree + gwesp).
# 100 cores: fit uses all for intensity cache; GOF runs 100 sims in parallel.
#
# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   mkdir -p cluster_output && sbatch inst/hypertext_conference/run_hypertext.slurm
#
# Full run: MAX_EDGES=10000, no truncation, N_GOF=100, MAX_ITER=5000 (defaults when under SLURM).
# Optional env (set before sbatch or in script):
#   MAX_ITER=5000  N_GOF=100  GOF_CORES=100

mkdir -p cluster_output

module load R

echo "=========================================="
echo "Hypertext Conference Hawkes Study"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $(hostname)"
echo "CPUs:      $SLURM_CPUS_PER_TASK"
echo "Memory:    $SLURM_MEM_PER_NODE MB"
echo "Started:   $(date)"
echo "=========================================="

# Ensure we're in the package root directory
PACKAGE_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [ ! -f "$PACKAGE_ROOT/DESCRIPTION" ]; then
  echo "ERROR: DESCRIPTION file not found. Please run from package root directory."
  echo "Current directory: $(pwd)"
  exit 1
fi

cd "$PACKAGE_ROOT"
echo "Working directory: $(pwd)"
echo "Git status:"
git status --short || echo "Not a git repository or git not available"

echo "=========================================="
Rscript inst/hypertext_conference/hypertext_conference_modern_fit.R

echo "=========================================="
echo "Finished:  $(date)"
echo "=========================================="
