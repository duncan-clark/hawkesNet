#!/bin/sh
#SBATCH --job-name=openalex_hawkes
#SBATCH --output=cluster_output/openalex_%j.out
#SBATCH --error=cluster_output/openalex_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
# Override cores: CORES_OVERRIDE=50 or 100.
#SBATCH --mem=128G

# OpenAlex Hawkes study: inhomogeneous fit, temporal fit, KS tests, GOF.
# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   mkdir -p cluster_output && sbatch inst/openalex_study/run_openalex.slurm
#
# Optional env (set before sbatch or in script):
#   OPENALEX_EMAIL=your@email
#   OPENALEX_STRING=Hawkes
#   OPENALEX_PAGES=100

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
echo ""
echo "NOTE: Gender prediction requires 'gender' and 'genderdata' packages."
echo "      If not installed, install BEFORE submitting this job (requires internet):"
echo "        Rscript -e \"install.packages(c('gender', 'genderdata'))\""
echo "=========================================="

# Ensure we're in the package root directory (should contain R/, inst/, DESCRIPTION)
# If running from package root, this will be the current directory
# If running from elsewhere, adjust the path
PACKAGE_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [ ! -f "$PACKAGE_ROOT/DESCRIPTION" ]; then
  echo "ERROR: DESCRIPTION file not found. Please run from package root directory."
  echo "Current directory: $(pwd)"
  echo "Looking for: $PACKAGE_ROOT/DESCRIPTION"
  exit 1
fi

cd "$PACKAGE_ROOT"
echo "Working directory: $(pwd)"
echo "Git status:"
git status --short || echo "Not a git repository or git not available"

# Pull latest code (optional - uncomment if you want to auto-pull)
# NOTE: Using devtools::load_all() loads from source code, so multiple jobs can run
# simultaneously without conflicts. However, you need to ensure the source code
# directory has the latest changes before submitting jobs.
# echo "Pulling latest code..."
# git pull origin paper_1_update || echo "Git pull failed or not needed"

echo "=========================================="
Rscript inst/openalex_study/openalex_hawkes_study.R

echo "=========================================="
echo "Finished:  $(date)"
echo "=========================================="
