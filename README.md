# hawkesGrowthNet

[![R-CMD-check](https://github.com/duncan-clark/hawkesGrowthNet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/duncan-clark/hawkesGrowthNet/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

**New users:** try the minimal examples first — see [**Try the examples first**](#new-users-try-the-examples-first) below.

**hawkesGrowthNet** is an R package that provides tools for simulating and analyzing networks that grow with **Hawkes process** arrival times.

- Generate events under various Hawkesian network growth formulations
- Estimate model parameters for hawkesGrowthNet models (homogeneous and inhomogeneous background rates)
- Goodness-of-fit (GOF) diagnostics with simulations matching the fitted model
- Visualize network growth over time (see `make_network_growth_animation`; requires suggested packages `networkDynamic` and `animation`)

**Note:** When using inhomogeneous background rates (KDE-based), GOF simulations automatically use `cond_intensity_inhom` with time-varying background to match the fitted model exactly.

---

## New users: try the examples first

Before running the full simulation studies, run the **minimal working examples** to simulate and fit one BA and one CS network (a few minutes total). From the package root in R or RStudio:

```r
# Install or load the package first (see Installation below)
# Then, with working directory = package root:
source("inst/examples/example_BA.R")   # BA (degree-weighted) model — ~1 min
source("inst/examples/example_CS.R")   # CS (change-statistic) model — ~2 min
```

Details: **[inst/examples/README.md](inst/examples/README.md)**. For full studies (many replicates, SLURM), see **Running the simulation studies** below.

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
