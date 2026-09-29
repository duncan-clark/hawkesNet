#!/bin/bash
# One-shot: git pull, load R, install package (no save workspace), sbatch chosen jobs.
# Usage: pull_install_submit.sh <study> [study2 study3 ...]
#        pull_install_submit.sh all
#
# Studies:
#   BA        - simulation study (Barabási–Albert)
#   CS        - simulation study (Change Statistics)
#   openalex  - OpenAlex Hawkes study
#   hypertext - Hypertext conference Hawkes study
#   all       - submit all of the above
#
# Examples:
#   ./pull_install_submit.sh BA CS          # submit two jobs
#   ./pull_install_submit.sh all            # submit all four jobs
#   ./pull_install_submit.sh hypertext      # submit one job (as before)
#
# Run from package root or from inst/.
set -e

resolve_slurm_script() {
  local study_lower="$(echo "$1" | tr '[:upper:]' '[:lower:]')"
  case "$study_lower" in
    ba)        echo "inst/simulation_study/run_BA.slurm" ;;
    cs)        echo "inst/simulation_study/run_CS.slurm" ;;
    openalex)  echo "inst/openalex_study/run_openalex.slurm" ;;
    hypertext) echo "inst/hypertext_conference/run_hypertext.slurm" ;;
    *)         echo "" ;;
  esac
}

ALL_STUDIES="BA CS openalex hypertext"

if [ $# -eq 0 ]; then
  echo "Usage: $0 <study> [study2 study3 ...]"
  echo "       $0 all"
  echo ""
  echo "Studies: BA, CS, openalex, hypertext, all"
  exit 1
fi

# Expand "all" and validate each study name
STUDIES=()
for arg in "$@"; do
  arg_lower="$(echo "$arg" | tr '[:upper:]' '[:lower:]')"
  if [ "$arg_lower" = "all" ]; then
    STUDIES=($ALL_STUDIES)
    break
  fi
  script="$(resolve_slurm_script "$arg")"
  if [ -z "$script" ]; then
    echo "ERROR: Unknown study '$arg'. Valid: BA, CS, openalex, hypertext, all"
    exit 1
  fi
  STUDIES+=("$arg")
done

# Find package root (directory containing DESCRIPTION and inst/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "$SCRIPT_DIR/../DESCRIPTION" ] && [ -d "$SCRIPT_DIR/../R" ]; then
  PACKAGE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
elif [ -f "$(pwd)/DESCRIPTION" ] && [ -d "$(pwd)/inst" ]; then
  PACKAGE_ROOT="$(pwd)"
else
  echo "ERROR: Run from package root (directory containing R/, inst/, DESCRIPTION) or from inst/"
  exit 1
fi

cd "$PACKAGE_ROOT"
echo "Package root: $PACKAGE_ROOT"
echo "Jobs to submit: ${STUDIES[*]}"
echo ""

echo "--- git pull ---"
git pull

echo "--- module load R ---"
module load R

echo "--- devtools::install() (no workspace save) ---"
R --no-save -e "devtools::install()"

# Submit each job
SUBMITTED=0
FAILED=0
for study in "${STUDIES[@]}"; do
  SLURM_SCRIPT="$(resolve_slurm_script "$study")"
  if [ ! -f "$PACKAGE_ROOT/$SLURM_SCRIPT" ]; then
    echo "WARNING: Slurm script not found: $SLURM_SCRIPT (skipping $study)"
    FAILED=$((FAILED + 1))
    continue
  fi
  echo "--- sbatch $SLURM_SCRIPT ---"
  sbatch "$SLURM_SCRIPT"
  SUBMITTED=$((SUBMITTED + 1))
done

echo ""
echo "Done. $SUBMITTED job(s) submitted, $FAILED skipped."
