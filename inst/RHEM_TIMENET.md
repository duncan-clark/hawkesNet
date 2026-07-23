# RHEM / timeNet dependency (experimental branch)

This branch (`rhem-experimental`) adds repeated-hit / transaction-network (RHEM-style) experiments **on top of JASA Paper 1 `main`** (package `hawkesNet`, tag `v0.1.0-jasa`).

## Required companion package: **timeNet**

Valued temporal change statistics used by the Ethereum / RHEM mark PMF backend need the separate R package **timeNet**:

- Local checkout (this project): `../timeNet` (sibling of this package)
- GitHub: `https://github.com/duncan-clark/timeNet`

Install before running timeNet-backed studies:

```r
# from GitHub
devtools::install_github("duncan-clark/timeNet")

# or from a local sibling checkout
devtools::install("../timeNet")
```

`DESCRIPTION` lists `timeNet` under **Suggests** and `Remotes: github::duncan-clark/timeNet`.

Without `timeNet`, simple/ERNM backends may still run; formulas that route to the `timeNet` backend will error with a clear missing-package message.

## Studies on this branch

- `inst/amlsim_study/` — AMLSim SAR / transaction experiments
- `inst/ethereum_study/` — Ethereum stablecoin dyad experiments (timeNet search)

Result folders (`amlsim_study_results_*`, `ethereum_study_results_*`) are gitignored; keep them outside the repo or sync separately.
