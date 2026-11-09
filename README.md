# hawkesNet

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

`hawkesNet` models when a network changes and what changes together. An event
is a complete update, which can contain several simultaneous edits. The
likelihood separates a normalized update probability from a ground intensity
for event timing. Simulation generates both the event times and updates.

## Installation

The source package currently supports macOS and Linux and requires R and a
C++17 compiler. Install the development
version, including its GitHub dependency `ernm`, with:

```r
install.packages("remotes")  # once
remotes::install_github("duncan-clark/hawkesNet")
library(hawkesNet)
```

For a local source checkout, run `R CMD INSTALL hawkesNet` from its parent
directory. Developers can instead use `pkgload::load_all("hawkesNet", debug = FALSE)`.
The CRAN graph import workflow also uses the `sqlite3` command-line utility.

## Start with a complete-update model

This small reference example generates births and edge additions/deletions on
an undirected graph. M2 lets later changes to an anchor's degree revise the
excitation from its earlier updates.

```r
library(hawkesNet)
model <- network_update_model(mark = "PA", timing = "M2",
                              summary = "anchor_threshold")
model$parameters  # explicit defaults on their natural scales
set.seed(12)
path <- simulate_network_updates(model, horizon = 5)
loglik_network_updates(path, model)[c("loglik", "timing", "marks")]

# A small interface demonstration: estimate only the background rate.
fit <- fit_network_updates(path, model,
                           fixed = setdiff(names(model$parameters), "mu"))
fit$parameters
fit$convergence  # zero indicates optimizer convergence
```

Omit `fixed` to estimate all active parameters. The default initial graph is a
single edge; supply a symmetric adjacency matrix through `initial` to choose
another graph. `path$times` and `path$updates` retain complete updates, including
their births, anchor and added/deleted targets. The likelihood includes survival
to `path$horizon`. A short run demonstrates the interface; it does not establish
parameter recovery. This dense reference implementation is for small networks.

The model has separate `gamma_B` and `gamma_Delta` coefficients. The kernel is
`K * exp(-beta * lag)`, with integrated mass `K / beta` under M0. Bounded
summaries give a non-explosive model. See the
[complete-update guide](inst/examples/NETWORK_UPDATES.md) for model choices,
custom paths, fitting and the feedback information check.

Runnable scripts are installed with the package:

```r
source(system.file("examples", "example_network_updates.R", package = "hawkesNet"))
```

## Choose the interface for your data

| Task | Entry points | Scope |
|---|---|---|
| Births and anchored edge turnover | `network_update_model`, `network_update_data`, `simulate_network_updates`, `fit_network_updates` | Undirected complete updates and separate initial/later feedback coefficients; small-network reference implementation. |
| Arbitrary grouped edit proposals | `poisson_proposals`, `sample_update`, `logprob_update` | Exact Poisson union laws for caller-defined admissible edit sets. |
| Addition-only growth example | `sim_hawkesNet`, `fit_hawkesNet`, `loglik_hawkesNet` | Addition histories stored in a `network` object; select the mark law explicitly. |
| Growth-history feedback | `hawkes_timing`, `hawkes_ground_intensity`, `hawkes_ground_compensator` | Bounded degree or triangle summaries and one feedback coefficient. |
| CRAN complete update law | `cran_update_data`, `fit_cran_updates`, `loglik_cran_updates` | Directed births, additions, deletions, shared personnel and two-component mixtures. |
| CRAN timing | `cran_timing_data`, `fit_cran_timing`, `loglik_cran_timing` | Calendar background, release types, initial topology and clipped later-topology feedback. |
| CRAN forward checks | `cran_graph_data`, `simulate_cran_updates`, `simulate_cran_joint` | Generate updates at fixed times, or generate both times and updates. |
| Event-time controls | `fit_temporal_hawkes`, `loglik_temporal_hawkes` | Unmarked temporal benchmarks. |

Use `help("cran_updates", package = "hawkesNet")` and the
[CRAN workflow guide](inst/cran/README.md) for input formats and the manuscript
specification. Its primary CS laws use degree, age, copying, migration and
shared-person effects. Transitive triples, open diamonds and shared-consumer
rectangles are held-out goodness-of-fit summaries. The complete archived
CRAN inputs are external to the package; installation alone does not download
or reconstruct the application dataset.

The [mark-law guide](inst/examples/MARK_PMF_CONSTRUCTIONS.md) distinguishes
single-edit proposals, paired proposals and the direct size-conditional law.
The [growth timing guide](inst/examples/TIMING_FEEDBACK.md) explains the
additions-only timing interface. Its single-coefficient M2 is a restricted
model. Use `network_update_model` or the CRAN timing interface for separate
initial- and later-topology coefficients.

## Reproducibility and development

An [R-CMD-check workflow template](inst/ci/R-CMD-check.yaml) is included for
Linux and macOS. To activate it, copy it to `.github/workflows/` using a GitHub
login with permission to manage workflows.

Keep the event definition, candidate sets, nonempty conditioning, time units,
and model options with each fitted object. A finite likelihood on observed
data does not by itself establish that a model can be simulated indefinitely.
Use bounded summaries or the manuscript's non-explosion conditions when
choosing feedback specifications.

The `PMF_mark_CS` modes `joint` and `size_conditional`, RHEM, and the dyadic
`joint_marked_hawkes` interfaces remain available for earlier studies. Their
model definitions differ from the primary CRAN specification. Scripts under
`inst/simulation_study` and the other historical study directories retain
those earlier designs. The `network_update_model` reference API follows the
current recovery study's birth/turnover mechanism; it does not bundle its
100-replicate experiment or archived numerical results.

From the package's parent directory:

```sh
R CMD build hawkesNet
R CMD check hawkesNet_0.2.0.tar.gz
```

For source-level checks:

```r
pkgload::load_all("hawkesNet", debug = FALSE)
testthat::test_dir("hawkesNet/tests/testthat")
```

Report reproducible issues at [GitHub Issues](https://github.com/duncan-clark/hawkesNet/issues).
