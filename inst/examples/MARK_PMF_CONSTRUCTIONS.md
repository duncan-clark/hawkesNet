# Update laws and their interfaces

One mark is a complete unordered update. Keep the same candidate sets, birth
law, age convention and nonempty conditioning in simulation and likelihood
evaluation.

## Poisson proposals in the manuscript

Construction 1 draws a Poisson number of single-edit proposals and retains
only distinct edits. If an edit has selection probability `p` and the proposal
mean is `omega`, its inclusion probability is `1 - exp(-omega * p)`. Edits are
independent before nonempty conditioning, given the proposed births and graph.
The proposal mean is not the expected number of observed edits.

Construction 2 draws Poisson proposals of edit sets, selecting each set by its
size and weight, then takes their union. Overlapping proposals induce dependence
between edits. `poisson_proposals()`, `sample_update()` and `logprob_update()`
provide a generic interface for caller-defined edit groups; see the
[complete-update guide](NETWORK_UPDATES.md). The CRAN `CS-pair` law implements
this construction for migration
pairs together with single additions and deletions. A recorded migration can
arise from either the pair channel or two single proposals; its likelihood sums
over both mechanisms. See the [CRAN guide](../cran/README.md) for fitting this law.

## Addition-only growth interfaces

These interfaces operate on `network` objects with vertex and edge addition
times. They do not represent deletions or revisions.

| Interface | Edge law | Meaning of `m` |
|---|---|---|
| `PMF_mark_BA` | Degree/age-weighted single proposals with Poisson collapsing | Proposal mean |
| `PMF_mark_CS(..., cs_mode="independent")` | Frozen exponential single-edge change-statistic weights, with Poisson collapsing | Proposal mean |
| `PMF_mark_CS(..., cs_mode="size_conditional")` | One direct draw of the complete size-`K` edge set; `K` is truncated Poisson | Poisson parameter before truncation |
| `PMF_mark_CS(..., cs_mode="joint")` | A collapsed-Poisson reference tilted by whole-update statistics | Reference proposal mean |

`size_conditional` is the direct whole-update variant described in the paper's
supplement. Its historical code/example name is CS-2. It does **not** apply the
Poisson-set union step of the manuscript's Construction 2. The `joint` mode is
also a separate earlier specification. Always supply `cs_mode` for new work:
the historical default selects `joint` when `m` is present and the older
legacy law otherwise.

The explicit CS modes use Poisson vertex births and condition the complete
birth/edge update to be nonempty by default. Birth-only updates are allowed.
`condition_nonempty=FALSE` is useful for inspecting a proposal law, but empty
updates cannot be recovered from an ordinary graph history.

## Runnable comparison

```r
library(hawkesNet)
source(system.file("examples", "example_mark_constructions.R",
                   package = "hawkesNet"))
result <- run_mark_construction_example("CS-1")
result$fit$fit_table
```

This fits only `m`. The function also accepts `"BA-1"` and the historical
`"CS-2"` direct variant. For a structural CS fit, supply an ERNM formula such as
`"edges + triangles + star(c(2,3))"` and matching `CS_params`. An edge-count
coefficient cancels from normalized exponential single weights and from
size-conditional weights, so fix it at zero. The fitter rejects this flat
parameter direction. Other coefficients require variation in their candidate
statistics to be identifiable.

## Support and computation

CS-1 evaluates every candidate edit against the same pre-edit graph. The direct
size-conditional mode instead evaluates each complete candidate set. That
calculation visits `choose(D,K)` sets for `D` candidate edges and realized size
`K`; the `joint` law enumerates all subsets. Both enumerators have a
`max_candidates=12` guard, configurable up to 20. The guard raises an error;
it does not remove candidates. CS-1 does not enumerate subsets.

`truncation` declares a candidate **node** window. It is part of the model,
not an edge-count cap or an approximation to the unrestricted likelihood.
It defaults to four nodes in explicit CS modes. `mark_decay="node_entrance"`
selects recent vertex indices and uses the older endpoint's birth age;
`mark_decay="activity"` selects vertices by recent activity and uses the most
recent activity of either endpoint. Edges already present, loops, repeated
edges, deletions and revisions of existing attributes are excluded.

The returned `log_density_func(params)` reevaluates the complete normalized
log probability. Prefer it to exponentiating tiny update probabilities.
In CS-1, `edge_probs` condition on births and applicable nonempty conditioning;
`reference_edge_probs` are before conditioning. In the direct size-conditional
mode, `edge_probs` also condition on the realized edge count, while
`edge_count_probs` gives that count's distribution. These are different
conditional quantities.
