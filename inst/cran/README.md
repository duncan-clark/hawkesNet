# CRAN application

The manuscript models one complete ten-minute network update at a time.
Package-major lines are vertices; edges are dependencies. Each update first
draws its releasing lines and births, then their simultaneous dependency edits.
An existing line must change at least one edge; a newborn line may have none.

## Primary specifications

| Manuscript name | `law` | Parameters, including releases and mixture |
|---|---|---:|
| PA | `"PA"` | 33 |
| CS | `"CS-1"` | 36 |
| CS-pair | `"CS-pair"` | 40 |
| TERGM | `"TERGM-logdegree"` | 36 |

Counts are free coefficients. Fitted vectors also retain fixed-zero slots
for compatibility; inspect `parameter_count` and `fixed_parameters` rather
than using the vector length.

The defaults are motif-free weights, a two-component mixture, exponential
CS weights and covariate proposal means. Supply dated personnel to CS,
CS-pair and TERGM. PA is the baseline without personnel terms.

Each mixture draws **one component for the whole update**, with shared target
preferences and component-specific release and edit rates. Each component is
conditioned separately before mixing. CS-pair proposes migration pairs that
delete an old-major dependency and add a higher major of the same family.
Its likelihood sums over single and paired proposals producing the observed
union. Shared personnel affect ordinary additions, deletions and migration
rates separately. Target degree, age, copying and release-class/out-degree
rate terms remain. No transitive-triple, diamond or shared-consumer term is
fitted; these motifs are forward goodness-of-fit targets.

```r
library(hawkesNet)
people <- readRDS("personnel.rds")
data <- cran_update_data("CRAN_CACHE", law = "CS-pair", people = people)
fit <- fit_cran_updates(data, window = "train", starts = 2)
fit$parameter_count
fit$predictive
saveRDS(fit, "CS-pair_train.rds")
```

`window="train"` estimates on 2020–2023; `"full"` estimates on 2020–2025.
Predictive contributions condition on the observed preceding history. Fitting
is unpenalized and reports boundary flags, convergence and a projected-score
check. Archived fit files, if available, supply warm starts; otherwise generic
starts are used. Existing fitted objects can be passed as `parameters`.

An installed command-line entry point fits all four laws in both windows:

```sh
Rscript /PATH/TO/INSTALLED/hawkesNet/cran/refit_primary.R CRAN_CACHE OUTPUT personnel.rds
```

Locate it with `system.file("cran", "refit_primary.R", package="hawkesNet")`.
This is a substantial full-data calculation; use a new output directory.

## Data contract

`cran_update_data()` imports the parameter-independent numerical cache used
by the paper. It does not download or reconstruct the CRAN archive. The cache
contains `marks/`, `actor_selection/`, `cs_popularity_20260924/`,
`cs2_paired_20260924/`, `marks.jsonl`, CSV snapshots in `summaries/` and a
read-only `../cran_history.sqlite` publication warehouse. Timing additionally
uses the family trigger/change cache and `current_fits/calendar_values.csv`.
Reading the warehouse requires the `sqlite3` executable. The numerical
NPY/NPZ files are read by R/C++; no Python runtime is required.

The full empirical data and reconstruction workflow are separate research
artifacts, not bundled with this package. The package's small tests and
examples run without them. Native memory-mapped candidate caches currently
support macOS and Linux.

Personnel data contain sorted integer `bin` values (UTC ten-minute intervals
since 1970), character `package` names, and a list column `people` of positive
integer person IDs. Each row replaces that family's personnel set. IDs should
represent the union of maintainers, authors and contributors after identity
resolution; matching names/emails and resolving missing or ambiguous identities
are preprocessing choices. No package code from metadata is executed.

Observed likelihoods use personnel through the current bin. Forward scenarios
freeze lists strictly before their origin; new majors inherit their family's
set and new families have empty sets. Edges may include more than one major of
a dependency family, so simulated support is broader than DESCRIPTION records.

## Timing

```r
timing_data <- cran_timing_data("CRAN_CACHE", release_types = TRUE)
# One value per day: 0 ordinary, 1 publication break, 2 reopening.
timing_data$calendar_phase <- publication_calendar_phase
clock_fit <- fit_cran_timing(timing_data, model = "M1b",
                             calendar_effects = TRUE)
```

The daily profile is fixed from 2019. `calendar_effects=TRUE` jointly estimates
the break multiplier and the separate reopening level. Supply phase labels
from the publication calendar; the package does not infer breaks from future
updates. For the paper's M2, use `model="M2", release_types=TRUE,
type_by_kernel=TRUE, calendar_effects=TRUE` (14 parameters). M1b has 12,
M1c 13 and M1a 8. `M0`, `one_exponential` and `Poisson` provide matched
unmarked controls (7, 5 and 3 with calendar effects). `batch` and `M1` preserve
the older size-only and initial-topology interfaces. See `?cran_timing` for
parameter names, bounds and data requirements.

## Simulation and checks

```r
initial <- cran_graph_data("CRAN_CACHE", "2025-01-01", people = people)
paths <- simulate_cran_updates(initial$state, fit, times = recorded_times,
                               checkpoints = max(recorded_times), seed = 1)
```

`simulate_cran_updates()` checks structure at supplied times.
`simulate_cran_joint()` generates times and updates together using a specified
daily calendar and complete pre-origin excitation history. It returns first
migration and motif-creation waits; use `directed_waits=TRUE` for the three
held-out directed motifs. A first creation compares the graphs before and
after the complete update, even when other motif instances are destroyed.
Simulation guards raise errors rather than recoding failed paths as censoring.

Timing fits return `timing`, `calendar` and `simulation_model` fields ready
for `simulate_cran_joint()`:

```r
forecast <- simulate_cran_joint(
  initial$state, history = pre_origin_history, fit = fit,
  timing = clock_fit$timing, calendar = clock_fit$calendar,
  origin = initial$origin, model = clock_fit$simulation_model,
  directed_waits = TRUE, seed = 1
)
```

A scalar graph summary is not sufficient to initialize M2; each historical
release retains its own summary and timestamp. Future calendar phases are
conditional scenario inputs. See `?simulate_cran_joint`.

## Archived exploratory specifications

The paper's structural checks call `gof.R` from the Git checkout. Earlier
experiments with motif preferences and paired additions use the scripts from
the [archived CRAN checkout](https://github.com/duncan-clark/hawkesNet/tree/1484437a3a491b875de3ad5d48ad4d9af2fbb393/inst/cran).

Saved diagnostic fits distinguish `legacy`, `reduced`, `pattern_preferences`,
`shared_consumer` and `pair_channels`. These explicit term selections retain
their single-component default.
For addition-pair experiments, `single`, `capped2` and `poisson` are different
proposal-count laws and require separate fits. The exact unrestricted union
likelihood may exceed its computation budget; it raises an error without
truncating candidates. Existing fits without a pair-law flag retain their
archived at-most-one-pair interpretation.
