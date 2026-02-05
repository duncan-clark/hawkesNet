# hawkesGrowthNet

[![R-CMD-check](https://github.com/duncan-clark/hawkesGrowthNet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/duncan-clark/hawkesGrowthNet/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

**hawkesGrowthNet** is an R package that provides tools for simulating and analyzing networks that grow with **Hawkes process** arrival times.

- Generate events under various Hawkesian network growth formulations
- Estimate model parameters for hawkesGrowthNet models
- Visualize network growth over time (see `make_network_growth_animation`; requires suggested packages `networkDynamic` and `animation`)

---

## Installation

You can install the development version of **hawkesGrowthNet** from GitHub:

```r
# install.packages("devtools")  # if needed
devtools::install_github("duncan-clark/hawkesGrowthNet")
```

---

## Running the simulation studies

After `git pull`, you can run the BA and CS simulation studies from an **interactive RStudio session** (e.g. in the cloud) or via **SLURM** on a cluster.

- **Interactive (RStudio):** From the package root, run `devtools::load_all()` then `source("inst/simulation_study/simulation_study_BA.R")` or `simulation_study_CS.R`.
- **SLURM:** From the package root, run `sbatch inst/simulation_study/run_BA.slurm` or `sbatch inst/simulation_study/run_CS.slurm`.

Full details (paths, options, outputs) are in **[inst/simulation_study/README.md](inst/simulation_study/README.md)**.
