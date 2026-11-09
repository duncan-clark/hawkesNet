# hawkesNet

`hawkesNet` models the timing and content of network updates. An event can
contain several simultaneous edits, such as a node birth and changes to its
incident edges. The model combines a Hawkes intensity for event times with a
probability law for each complete update.

## Install

The package builds from source on macOS and Linux with R and a C++17 compiler.

```r
install.packages("remotes")
remotes::install_github("duncan-clark/hawkesNet")
```

## Simulate and fit a small network

The example starts with a single undirected edge. The update law proposes
births, additions and deletions. With M2 timing, an anchor's degree can change
the excitation from its past updates.

```r
library(hawkesNet)

model <- network_update_model(mark = "PA", timing = "M2",
                              summary = "anchor_threshold")
set.seed(12)
path <- simulate_network_updates(model, horizon = 5)
loglik_network_updates(path, model)[c("loglik", "timing", "marks")]

# Estimate the background rate on this short path.
fit <- fit_network_updates(path, model,
                           fixed = setdiff(names(model$parameters), "mu"))
fit$parameters
fit$convergence
```

Use `network_update_data()` to supply observed updates and their observation
horizon. Choose `mark = "CS-1"` for change-statistic proposals, and omit `fixed`
to estimate all active parameters. The guide covers the M0–M2 timing choices,
the feedback information check and a runnable installed example.
[Complete-update guide](inst/examples/NETWORK_UPDATES.md)

## Specify an update law

`poisson_proposals()` takes admissible edit groups and their Poisson rates.
An update contains the union of groups proposed at least once. Use
`sample_update()` to draw an update and `logprob_update()` to evaluate its
probability, including the contributions of overlapping groups.
[Grouped proposals](inst/examples/NETWORK_UPDATES.md#arbitrary-poisson-edit-set-proposals)

## Refit the CRAN application

Use `cran_update_data()` and `fit_cran_updates()` for the directed
package-dependency network. The workflow also provides release-type Hawkes
timing, conditional forward simulations and joint forecasts. It uses prepared
CRAN caches and a publication warehouse. Graph import requires the `sqlite3`
command-line utility. [CRAN workflow guide](inst/cran/README.md)

## Archived studies

The previous package revision and study instructions cover addition-only
growth models, older CS mark laws, RHEM and dyadic fits. Use that revision to
run those studies or replay saved fits. [Archive](archive/README.md)
