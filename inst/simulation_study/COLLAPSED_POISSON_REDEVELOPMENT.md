# Simulation study after the mark-law correction

> Version note (September 22, 2026): the CS law discussed below is the earlier `cs_mode="joint"` model. New CS-1 and CS-2 examples have explicit modes; see [the construction guide](../examples/MARK_PMF_CONSTRUCTIONS.md). These historical study instructions and outputs do not validate those new laws.

Historical BA/CS simulations and estimates were generated under the former mark
law. They do not validate the collapsed-Poisson specification. Preserve those
files and write every new study to a distinct output directory. For revised BA,
`m` is the expected number of attachment attempts, not the expected number of
distinct edges. For revised CS it is the attempt rate of the collapsed-Poisson
reference distribution, before whole-update interaction and conditioning; it is
not the final mean number of edges. The reference candidate weights remain fixed
within an update. The CS law tilts that reference using the exact unordered
whole-update statistic, so its final edge indicators are dependent.

## Bounded smoke and conditional recovery run

Use the source checkout rather than an older installed package:

```sh
Rscript inst/simulation_study/smoke_collapsed_poisson.R \
  --output-dir=/tmp/hawkesnet-collapsed-smoke --models=BA
```

The runner can be called from any working directory using its absolute path.
It requires an explicit empty/new output directory, loads the package using
`pkgload::load_all()`, uses one core, and records package source hashes, parameters,
seeds and session information. Defaults are two replicates at each of `T=2,4`.
Use `--models=BA,CS` to include the revised CS model. Its four-vertex candidate
window has at most six candidate edges; `max_candidates=12` enforces a hard
enumeration bound. This small-support experiment is not the original 100-node
study and must not be presented as validating that larger specification.

The default `--fit=mark_scale` estimates only `m`, with the remaining parameters
fixed exactly at their generating values. This checks a targeted part of the new
likelihood and is not a demonstration of joint parameter recovery. An exploratory
joint fit is available with `--fit=joint`; its CS edges coefficient is fixed at the
generating value. Tiny samples can have weak identification and boundary fits.
No standard errors or consistency claims are produced by this runner.

The checks include exhaustive BA mark normalization on a tiny network, cached vs
direct likelihood agreement, finite simulation/refitting, and agreement between
the optimizer objective and a fresh fitted likelihood. All attempted replicates,
warnings and optimization failures are saved. Event-count reconciliation detects
an unrecorded empty-update problem instead of silently fitting fewer events.
Parameters at optimization bounds are named explicitly, even when the optimizer
reports convergence; convergence alone is not a recovery or identifiability check.
An iteration-limit exit is recorded and makes the run unsuccessful; it is not
removed from the report. Generated estimates are descriptive smoke output only.

The runner writes `checks.csv`, `attempts.csv`, `estimates.csv`, a source/session
manifest, and a result RDS per attempted replicate. It does not read old study
results, modify paper figures, install packages, or launch the cluster scripts.

## Conditions for a new full study

1. Validate the exact CS law: the collapsed-Poisson reference is tilted by the
   statistic of the whole unordered update. Sequential toggles are only an exact
   way to calculate this statistic, never a latent event order. The tilt requires
   a normalizer over possible update subsets. Its current exact implementation
   enumerates those subsets and scales as `2^D` for `D` candidate edges. The old
   100-node window can have 4,950 candidate edges and is computationally impossible
   with this enumerator. A larger full CS study is blocked until a scientifically
   justified bounded candidate design or a validated scalable normalizer is
   available. Do not silently replace the interacting law with independent edges.
2. Decide how all events are represented, including zero-node/zero-edge updates,
   empty candidate sets, the first event and the seed network. Verify that fitting
   reconstructs exactly the event history simulated. State any conditioning rule.
3. Specify candidate support and truncation as part of the model, using identical
   settings in simulation and likelihood. Reassess identifiability of normalized
   CS weights (especially a common edges term) using likelihood profiles and
   multiple starts before treating every coefficient as an estimand.
4. Run increasing observation horizons with prespecified replicates and distinct
   deterministic seeds per `(model, horizon, replicate)`. If using parallel
   workers, initialize independent RNG streams or seed each task; seeding only the
   parent process does not fix worker randomness. Keep fixed coefficients equal
   to their intended values rather than perturbing them with the free starts.
5. Save every attempted replicate, error, optimizer status, boundary estimate and
   objective value. Report failure rates and unconditional estimation diagnostics
   alongside conditional-on-convergence summaries. Do not filter estimates by
   closeness to truth or arbitrary magnitude thresholds to make RMSE decrease.
6. Report bias, RMSE, uncertainty of those Monte Carlo summaries, interval
   coverage where standard errors are trustworthy, and run-time scaling. An
   observed RMSE curve need not be monotonically decreasing with finite replicates;
   parameter recovery over a grid is evidence, not a proof of consistency.
7. Use time-only, mark-only and simpler network-update baselines, and evaluate
   held-out time/mark prediction and graph calibration. Separate the effects of
   additional event exposure from changes in network size and candidate support.
8. Keep temporal parameterizations consistent: the network model uses the kernel
   `K * exp(-beta_overall * lag)`, with branching ratio `K / beta_overall`.
   The separate temporal-fitting code uses a normalized exponential kernel, so its
   `K` denotes the branching ratio. Transform parameters before comparing fits or
   reporting a common recovery target. Use the corrected temporal simulator for
   all new runs; a correct mark PMF alone cannot validate the old timing algorithm.

The existing long study scripts need these changes before reuse. In particular,
the current CS consistency block perturbs `CS_params1` before holding it fixed,
both studies drop failed runs through `NULL`, and their `keep` rules exclude
large estimates. The older `diagnose_RMSE_consistency.R` does not map the true `m`
value and sets monotone RMSE as a target. Those scripts and historical results are
left intact here; the bounded runner provides a clean initial validation path.
