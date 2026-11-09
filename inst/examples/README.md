# Runnable example

```r
library(hawkesNet)
source(system.file("examples", "example_network_updates.R", package = "hawkesNet"))
```

The script simulates a network with births, edge additions and edge deletions.
It evaluates the likelihood, fits the background rate and shows a grouped
proposal law. It runs without external data and prints the optimizer's
convergence code.

The complete-update guide explains the model, how to supply observations and
how to combine Poisson proposals for overlapping edit groups.
[Complete-update guide](NETWORK_UPDATES.md)

The CRAN workflow covers directed turnover, paired migrations, shared personnel,
release-type timing and forward simulation. [CRAN guide](../cran/README.md)
