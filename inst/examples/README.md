# Small runnable examples

After installation, run an example from any working directory:

```r
library(hawkesNet)
source(system.file("examples", "example_network_updates.R", package = "hawkesNet"))
```

For a source checkout, load that checkout with `pkgload::load_all("hawkesNet", debug = FALSE)`
and source `hawkesNet/inst/examples/example_network_updates.R`.

| Script | What it demonstrates |
|---|---|
| `example_network_updates.R` | Births, anchored edge additions/deletions, separate feedback coefficients, and a small fit. Start here. |
| `example_BA.R` | Preferential-attachment growth, M0 event times and a one-parameter fit. |
| `example_CS.R` | Explicit CS-1 exponential weights, complete-update simulation and a one-parameter fit. |
| `example_mark_constructions.R` | A reusable function comparing single-proposal growth laws and the historical direct size-conditional variant. Sourcing defines the function; call it to run. |
| `example_timing_feedback.R` | Bounded degree feedback on an undirected addition history, same-history timing comparisons and rescaled waits. |

These small examples are interface demonstrations. Their restricted fits are
not parameter-recovery studies. They do not fetch data or write output files.
Fitting errors are allowed to surface, and optimizer convergence is reported.

The [complete-update guide](NETWORK_UPDATES.md) covers the reference turnover
model and generic grouped proposals. The [mark-law guide](MARK_PMF_CONSTRUCTIONS.md) explains candidate support and
normalization. The [timing guide](TIMING_FEEDBACK.md) gives the scope of the
single-coefficient growth interface. For directed turnover, paired migrations,
shared personnel, and the separate initial/later feedback coefficients used
in the application, see the [CRAN workflow](../cran/README.md).
