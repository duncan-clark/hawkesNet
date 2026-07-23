#!/usr/bin/env bash
# Source from SLURM scripts (from package root):
#   source inst/cluster_env.sh
#
# Override by exporting HAWKESNET_OUTPUT_DIR before sbatch if needed.

if [[ -z "${HAWKESNET_OUTPUT_DIR:-}" ]]; then
  # Package root = directory containing DESCRIPTION (SLURM_SUBMIT_DIR or cwd)
  _PKG_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
  if [[ ! -f "${_PKG_ROOT}/DESCRIPTION" ]]; then
    _PKG_ROOT="$(pwd)"
  fi
  # Durable NeSI/project default: sibling of the git clone
  export HAWKESNET_OUTPUT_DIR="$(cd "${_PKG_ROOT}/.." && pwd)/cluster_output"
fi

mkdir -p "${HAWKESNET_OUTPUT_DIR}"/{runs,logs,paper_figures,diagnostics}
echo "HAWKESNET_OUTPUT_DIR=${HAWKESNET_OUTPUT_DIR}"
