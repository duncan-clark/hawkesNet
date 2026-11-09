# Complete-update reference workflow

Run the installed example for a short simulation and fit.

```r
library(hawkesNet)
source(system.file("examples", "example_network_updates.R", package = "hawkesNet"))
```

`network_update_model()` specifies a small undirected birth/turnover model.
At each event it draws Poisson births, chooses an anchor, and makes independent
Poisson proposals for adding and deleting incident edges. The edge draw is
conditioned on at least one edit; the birth count keeps its unconditional
Poisson distribution. Edits are applied simultaneously. New vertices receive
successive integer labels. There are no vertex deletions.

## Model choices

```r
library(hawkesNet)
model <- network_update_model(mark = "CS-1", timing = "M2",
                              summary = "anchor_threshold")
model$parameters
```

`mark="PA"` uses target-degree weights. `mark="CS-1"` uses exponential weights
for signed changes in triangles, two-stars divided by 6, and three-stars
divided by 18. Addition and deletion coefficients are separate. All candidate
statistics are computed before the complete update. These simulation-study
CS terms differ from the motif-free CRAN application specification.

M0 assigns every update amplitude one. M1 uses
`exp(gamma_B * initial_summary)`. M2 additionally multiplies by
`exp(gamma_Delta * clipped_summary_change)`. Summaries refer to the same anchor
throughout its update's excitation history. Choices are `anchor_threshold`,
`degree_exposure` (anchor degree divided by `scale + degree`), and `log_degree`.
The last is unbounded; the reference simulator requires nonpositive `gamma_B`
and finite M2 clipping as a sufficient non-explosion condition. This restriction
is conservative. Evaluation of a finite path's likelihood is still possible
outside that simulation restriction.

Defaults use the recovery study's coefficients, but the default initial graph
is one edge, not the study's random graph on 100 vertices. `repeat_probability`
is a known design choice, not a fitted parameter. It defaults to 0.9. The
threshold, scale, clipping and repeat probability stay fixed during fitting.
The reference fitter uses strictly positive baseline, kernel and proposal
parameters; its search bounds are reported with the fit.

## Data, likelihood and fitting

```r
set.seed(12)
path <- simulate_network_updates(model, horizon = 5)
path$updates[[1]]
loglik_network_updates(path, model)
```

The path contains `times`, `updates`, `initial` and `horizon`. For external data,
use `network_update_data(times, updates, initial, horizon)`. Each update is a
list with `births`, `anchor`, `add` and `delete`; the last two are integer target
vectors. The initial graph must be a symmetric zero-diagonal 0/1 adjacency
matrix with at least two vertices. Times are strictly increasing in
`(0, horizon)`. Likelihood replay checks that all updates are admissible.

`loglik_network_updates` returns the total log likelihood, timing and update
contributions, event intensities and the integrated intensity. It includes
survival to the observation endpoint.

`fit_network_updates(path, model)` fits all active parameters. Use `fixed` to
hold named parameters constant and `starts` to supply a list of complete named
parameter vectors. The result includes the fitted `model`, `parameters`,
`loglik`, optimizer `convergence`, `boundary` parameters and numerical `bounds`.
With `hessian=TRUE` (the default), `standard_errors` and `intervals` contain
local observed-information uncertainty for free parameters when the fit is
interior, converged and has suitable information. Otherwise they are missing;
these approximations have no general guarantee for a growing network.
Inspect convergence and boundary flags and compare starting values. Local
convergence alone does not establish identification or a global maximum.

After fitting M1, `feedback_score(path, fit$model)` computes the later-feedback
score, efficient outer-product information, and predicted standard error near
zero. It assumes an interior M1 timing maximum. A rank-deficient nuisance
information matrix yields missing information and an explanation. A short
path often contains too little information for this check to be useful.

The implementation uses explicit adjacency matrices and quadratic timing
storage. It is a transparent small-network reference, not the accelerated
engine for the complete CRAN archive. Event/vertex computation guards raise
errors rather than returning truncated simulation paths.

## Arbitrary Poisson edit-set proposals

`poisson_proposals(groups, rates)` implements both manuscript constructions.
Each group has an independent Poisson count; the update is the union of groups
whose counts are positive. For a size-conditional channel, set a group's rate
to `omega * pi[size] * weight / sum(weights_of_that_size)`. Concatenate groups
and rates to combine channels. Define every group on the same pre-update graph.

```r
law <- poisson_proposals(
  list("delete_old", "add_new", c("delete_old", "add_new")),
  c(0.2, 0.05, 0.15)
)
sample_update(law)
exp(logprob_update(law, c("delete_old", "add_new")))
```

Edit identifiers are strings; the caller is responsible for their graph
meaning and admissibility. The probability sums over all compatible proposal
combinations. `nonempty=TRUE` conditions the whole union on being nonempty.
Large overlapping groups can make exact evaluation expensive; the configurable
`max_component` guard raises an error without changing model support.
