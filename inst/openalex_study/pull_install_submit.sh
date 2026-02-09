#!/bin/bash
# One-shot: git pull, load R, install package (no save workspace), sbatch chosen job.
# Usage: pull_install_submit.sh <BA|CS|openalex>
#   BA      - simulation_study (Barabási–Albert)
#   CS      - simulation_study (Change Statistics)
#   openalex - OpenAlex Hawkes study
# Run from package root or from this script's directory.
set -e

STUDY="${1:-}"
STUDY_LOWER="$(echo "$STUDY" | tr '[:upper:]' '[:lower:]')"
case "$STUDY_LOWER" in
  ba)       SLURM_SCRIPT="inst/simulation_study/run_BA.slurm" ;;
  cs)       SLURM_SCRIPT="inst/simulation_study/run_CS.slurm" ;;
  openalex) SLURM_SCRIPT="inst/openalex_study/run_openalex.slurm" ;;
  *)
    echo "Usage: $0 <BA|CS|openalex>"
    echo "  BA       - simulation study (Barabási–Albert)"
    echo "  CS       - simulation study (Change Statistics)"
    echo "  openalex - OpenAlex Hawkes study"
    exit 1
    ;;
esac

# Find package root (directory containing DESCRIPTION and the slurm script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
echo "Study: $STUDY_LOWER -> $SLURM_SCRIPT"

echo "--- git pull ---"
git pull

echo "--- module load R ---"
module load R

echo "--- devtools::install() (no workspace save) ---"
R --no-save -e "devtools::install()"

if [ ! -f "$PACKAGE_ROOT/$SLURM_SCRIPT" ]; then
  echo "ERROR: Slurm script not found: $PACKAGE_ROOT/$SLURM_SCRIPT"
  exit 1
fi
echo "--- sbatch $SLURM_SCRIPT ---"
sbatch "$SLURM_SCRIPT"

echo "Done. Job submitted."
