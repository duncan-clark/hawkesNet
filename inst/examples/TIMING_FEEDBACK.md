# Network feedback in event timing

`timing = hawkes_timing()` selects M0 by default. Pass the same `timing` object to
simulation, likelihood, fitting, conditional-intensity evaluation and residual
calculations. Numerical feedback strength belongs in `params$feedback_gamma`
(default zero), so it can be fitted alongside the existing numerical parameters.

| Specification | Meaning |
| --- | --- |
| `hawkes_timing("M0")` | Each previous event contributes `K * exp(-beta_overall * age)`. |
| `hawkes_timing("M1", feature = "degree")` | Multiply each contribution by alpha evaluated immediately after its complete update; retain that value forever. |
| `hawkes_timing("M2", feature = "degree")` | Reevaluate the structural feature of every previous complete update in the graph immediately before the evaluation time. |

The mark is the **complete simultaneous network update**, including all new
vertices and edges. Features average over that update; edges are not treated as
separate Hawkes arrivals. Simulation samples one update from the supplied PMF
at each accepted event, updates the graph immediately and uses that graph to
generate subsequent times. The full PMF, such as `PMF_mark_BA`, is retained in
the marked likelihood.

## Two feasible alpha specifications

**Degree importance with a bounded feature.** Use
`hawkes_timing("M1", feature = "degree", scale = s)` or the M2 equivalent.
Let `dbar` be mean current degree across the distinct vertices touched by the
update: its new vertices together with endpoints of all its added edges.
The feature and default link are

```
x = dbar / (s + dbar)
alpha = exp(feedback_gamma * x)
```

`s > 0` is a fixed, specified degree scale: `dbar = s` gives `x = 1/2`.
For an update touching no vertices, the feature is zero. Alpha equals one at
zero feedback, is positive, and is bounded between `min(1, exp(gamma))` and
`max(1, exp(gamma))` for any fixed finite gamma. Positive gamma gives more
excitation to updates involving well-connected vertices. Under M2, later growth
can strengthen an earlier peripheral update. A difference of `delta` in the
feature gives an excitation-amplitude ratio of `exp(gamma * delta)`.

**Triangle participation with a bounded link.** Use
`hawkes_timing("M1", feature = "triangles")` or the M2 equivalent.
The feature is the fraction of the update's added edges that currently belong
to at least one triangle. A node-only update has feature zero. The default is

```
x = number of added edges currently in a triangle / number of added edges
alpha = 2 * plogis(feedback_gamma * x)
```

This alpha lies strictly between zero and two and equals one at zero gamma or
zero triangle participation. Positive gamma increases the influence of an
update whose edges are embedded in triangles; negative gamma decreases it.
Under M2, closing a triangle can change the influence of older updates even
when those updates initially created no triangle. The feature averages over
the entire added-edge collection, including triangles completed by multiple
edges in the same update.

These are modeling choices, not established empirical improvements. Degree
measures involvement in well-connected vertices; triangle participation
measures local closure. Fix the feature and scale before comparison, or choose
them using training data. Fit gamma only when the feature varies sufficiently;
nearly constant alpha is confounded with K, and a graph with no triangles
cannot identify triangle feedback. Compare M0, M1 and M2 on the same event
definition and update PMF, with held-out timing calibration, predictive scores
and uncertainty. `link = "exp"` or `link = "logistic"` can explicitly override
the default link while retaining the selected feature.

## APIs and interpretation

`sim_hawkesNet`, `loglik_hawkesNet`, `fit_hawkesNet`, `cond_intensity`,
`cond_intensity_inhom`, `compensators_hawkesNet` and `ks_test_pval_hawkesNet`
accept `timing`. `gof` defaults to `fit$timing` when available and retains the
fitted gamma; an explicit `timing` overrides it. Supply the same CS `cs_mode`
and `max_candidates` used in fitting when generating CS goodness-of-fit
replicates. GOF validates the fitted numerical parameters without clipping K or
changing valid zero excitation, edge-decay or vertex-birth rates.
Ground-only evaluation is also available:

```r
timing <- hawkes_timing("M2", feature = "degree", scale = 2)
params$feedback_gamma <- 0.5
hawkes_ground_intensity(t, params, mark_filtration = net, timing = timing)
hawkes_ground_compensator(params, c(start, end), net, timing = timing)
```

`hawkes_ground_intensity` accepts a vector of times and uses only events strictly
before each time. `hawkes_ground_compensator` returns the integrated ground
intensity over the requested interval, including the effect of any supplied
prehistory. The log-likelihood interface retains its existing requirement that
observed event times lie inside its observation window. `compensators_hawkesNet`
returns cumulative integrated intensity from the observation-window start to
each event; residual waiting times are `diff(c(0, compensator))`.

The code uses **K as exponential-kernel amplitude**, not branching ratio:
the M0 kernel has integral `K / beta_overall`. M1 adds a mark-specific factor;
M2 changes old contributions when the graph changes, so the M0 ratio alone is
not a branching-ratio interpretation for feedback models. Bounded alpha is
useful for controlling rates but does not, by itself, verify the paper's
full marked-history contraction assumption.

Homogeneous M1/M2 simulation uses sequential thinning with the exact current
graph state. Inhomogeneous simulation requires the existing
`inhom_bg$mu_fit$mu_fun` baseline function and a supplied finite global upper
bound, `inhom_bg$mu_fit$mu_bound` (or `inhom_bg$mu_bound`), valid over the entire
simulation window. A baseline evaluated at a grid is insufficient to certify
such a bound. Inhomogeneous fitting uses `mu_vec` and `integral_bg` as before;
both must describe the same baseline and observation window.

This implementation is a reference computation for **simple undirected
networks with recorded vertex and edge addition times**. Directed networks,
loops, parallel edges, edge deletions and edge revisions are outside the
implemented feedback feature support. A complete-update PMF can be richer
scientifically, but these feature functions cannot reconstruct unsupported
changes from the current stored growth filtration. Use an explicit extended
history representation and additional feature implementation for those models.

M2 recalculates old-update importance across observed graph states. Its
event-history storage/work can be quadratic, and graph feature computation can
cost more, especially for triangles. Cached likelihood closures avoid repeated
mark-PMF setup and must still reevaluate alpha when gamma changes. This path
is intended for small and moderate comparisons; it is not a scalable fit
implementation for very large networks.

See `example_timing_feedback.R` for a short BA simulation, same-history model
comparison, gamma-only fit and residual calculation.
