#!/bin/sh
#SBATCH --job-name=hypertext_hawkes
#SBATCH --output=cluster_output/hypertext_%j.out
#SBATCH --error=cluster_output/hypertext_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=32G

# Hypertext 2009 conference: hawkesNet fits (3 models, same formula).
# Formula A (tri+star), varying what is fixed vs free:
#   Fit A: K=0.5 fixed, edges=-5 fixed
#   Fit B: K=0.5 fixed, edges free
#   Fit C: K free, edges free
#
# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   mkdir -p cluster_output && sbatch inst/hypertext_conference/run_hypertext.slurm
#
# Full run defaults: N_GOF=25, MAX_ITER=5000, GOF for all 3 fits.
# Optional env (set before sbatch or in script):
#   MAX_ITER=5000  N_GOF=25  GOF_CORES_OUTER=25
#   RUN_FIT_A=TRUE  RUN_GOF_A=TRUE  (etc. for B, C)

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
export N_GOF=25
Rscript inst/hypertext_conference/hypertext_fits.R 2>&1

echo "=========================================="
echo "Finished:  $(date)"
echo "=========================================="
