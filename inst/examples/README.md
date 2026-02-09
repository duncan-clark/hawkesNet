# Minimal working examples

**Try these first** after installing the package. Short, runnable scripts that simulate one network and fit the model so you can see the workflow in under a few minutes.

Run from the **package root** (directory containing `R/`, `inst/`, `DESCRIPTION`). In R or RStudio:

```r
# From package root (e.g. setwd("/path/to/hawkesNet") or open project)
source("inst/examples/example_BA.R")   # Barabási–Albert (degree-weighted) model
source("inst/examples/example_CS.R")  # Change-statistic (ERGM-style) model
```

| Script | Model | What it does | Approx. time |
|--------|--------|----------------|--------------|
| **example_BA.R** | BA (degree-weighted attachment) | One simulation + one fit; `m` = expected edges per event (Poisson); prints true vs fitted parameters | ~1 min |
| **example_CS.R** | CS (change statistics / ERGM-style) | One simulation + one fit (K fixed); prints true vs fitted parameters | ~2 min |

For full simulation studies (many replicates, consistency, explosive regime, SLURM), see [inst/simulation_study/README.md](../simulation_study/README.md).
