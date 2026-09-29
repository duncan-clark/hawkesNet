# Exact mark PMF constructions

Use the same explicit mode, candidate window, age convention and node-arrival cutoff in simulation, likelihood evaluation and fitting. Load this source checkout with `pkgload::load_all("hawkesNet")` from the repository root.

| Paper example | Code | Edge-set law | Interpretation of `m` |
|---|---|---|---|
| BA-1 | `PMF_mark_BA` | Fixed degree/age weights, Poisson attempts, collapsed independent inclusion | Mean attempts |
| CS-1 | `PMF_mark_CS(..., cs_mode="independent")` | Fixed single-edge logistic change-statistic/age weights, Poisson attempts, collapsed independent inclusion | Mean attempts |
| CS-2 | `PMF_mark_CS(..., cs_mode="size_conditional")` | Poisson count conditioned on `K <= D`, then whole-update structural/age weights within size `K` | Count-law Poisson parameter, before truncation |
| Earlier redevelopment CS | `PMF_mark_CS(..., cs_mode="joint")` | Collapsed-Poisson temporal reference tilted by whole-update statistics over all sizes | Reference attempt mean |

Both new CS laws propose Poisson node births and condition the complete node/edge update to be nonempty by default. Independence in CS-1 holds before that conditioning, conditional on the graph and proposed births. A birth-only event remains nonempty. `condition_nonempty=FALSE` is useful for inspecting the proposal law, but invisible empty updates are unsuitable for ordinary network-history simulation/estimation.

The earlier default remains unchanged for reproducibility: omitting `cs_mode` chooses `joint` when `m` is present, otherwise the pre-redevelopment legacy implementation. Explicit new modes require a parameter named exactly `m`; the temporal parameter `mu` is not substituted. The older no-`m` implementation is retained for historical use and is not validated as one of the rewritten examples. Two former simulation/fit tests now exercise CS-1 explicitly rather than relying on that older path.

## Reproducible short example

```r
pkgload::load_all("hawkesNet", quiet=TRUE)
source("hawkesNet/inst/examples/example_mark_constructions.R")
result <- run_mark_construction_example("CS-2")
result$fit$fit_table
```

The example estimates only `m` and is an interface smoke check. Select `"BA-1"` or `"CS-1"` to exercise the other construction.

For a structural CS fit, use an ERNM formula such as `"edges + triangles + star(c(2,3))"`, provide the matching `CS_params`, and fix `CS_params1=0` via `fixed_params="CS_params1"`. In CS-2 the edge coefficient cancels exactly conditional on size; the fitter refuses to optimize this flat direction. Node-only structural coefficients also cancel conditional on births. Other parameters can still be weakly identified on a given support. In CS-1, triangle effects need candidate edges whose frozen single-edge triangle changes vary.

## Calculation and support

CS-1 computes each edge's statistics against the same baseline, restoring the graph between calculations; its PMF requires no subset enumeration. CS-2 evaluates complete size-K sets using telescoping ERNM calculations, without any generative edge order. Its enumerator visits `choose(D,K)` states, not all `2^D` states. The implementation retains a `max_candidates=12` guard (configurable up to 20) for CS-2 and the earlier joint model. Errors do not silently shrink the scientific candidate set. This guard does not apply to CS-1.

`truncation` is a declared candidate *node* window, not an edge-count cap. It defaults to four nodes for nonlegacy CS modes. `mark_decay="node_entrance"` selects most recent vertex indices and uses the older-index endpoint's birth age; `"activity"` selects vertices by recent activity and uses the most recent activity of either endpoint. Edges already present, loops, deletions, repeated edges and changes to modeled attributes of existing nodes are not valid updates for the new CS implementations.

For CS-2, returned `edge_probs` condition on both the realized node and edge counts; `edge_count_probs` describes the separate count distribution conditional on births and applicable nonempty conditioning. For CS-1, `edge_probs` condition on births and applicable nonempty conditioning, while `reference_edge_probs` are before nonempty conditioning. These quantities should not be compared as if they condition on the same information.

Each mode returns reusable numeric `log_density_func(params)` and `density_func(params)` caches; candidate graphs and ERNM external pointers are not retained in the new caches. Parameter changes reevaluate the correct normalizer. Categorical birth attributes use exact simplex probabilities, including zero probability for forbidden levels.

## Evidence and remaining work

The September 22 validation directory contains exhaustive small-support tests and nine matched-law simulation/refit checks for BA-1/CS-1/CS-2. Every smoke attempt is retained with source hashes and seeds. These short runs establish implementation compatibility, not consistency. The September 21 recovery plots, including the roughly 5,000-event CS runs, concern the **earlier joint CS model**. They must not be relabeled as CS-1 or CS-2 results.

The older study-specific accelerated CS fitter implements only that earlier joint law. Use the production `fit_hawkesNet` with an explicit mode for the new models until equivalent accelerated evaluators are separately derived and validated.
