#!/bin/bash
# One-shot: git pull, load R, install package (no save workspace), sbatch OpenAlex job.
# Run from package root (directory containing R/, inst/, DESCRIPTION), or from this script's directory.
set -e

# Find package root (directory containing DESCRIPTION and inst/openalex_study/run_openalex.slurm)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Script lives in inst/openalex_study/, so package root is two levels up
if [ -f "$SCRIPT_DIR/../../DESCRIPTION" ] && [ -f "$SCRIPT_DIR/run_openalex.slurm" ]; then
  PACKAGE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
elif [ -f "$(pwd)/DESCRIPTION" ] && [ -f "$(pwd)/inst/openalex_study/run_openalex.slurm" ]; then
  PACKAGE_ROOT="$(pwd)"
else
  echo "ERROR: Run from package root (directory containing R/, inst/, DESCRIPTION) or from inst/openalex_study/"
  exit 1
fi

cd "$PACKAGE_ROOT"
echo "Package root: $PACKAGE_ROOT"

echo "--- git pull ---"
git pull

echo "--- module load R ---"
module load R

echo "--- devtools::install() (no workspace save) ---"
R --no-save -e "devtools::install()"

echo "--- sbatch inst/openalex_study/run_openalex.slurm ---"
sbatch inst/openalex_study/run_openalex.slurm

echo "Done. Job submitted."
