# Package audit, November 2026

This audit used the current manuscript as a specification and changed package
code, documentation and tests. The manuscript and its reported analyses were
not edited or rerun.

## Supported workflows

| Manuscript method | Package entry point |
|---|---|
| Independent Poisson single-edit proposals | `poisson_proposals()` with singleton groups |
| Size-conditional Poisson edit-set proposals and their unordered union | `poisson_proposals()`, `sample_update()`, `logprob_update()` |
| PA and CS-1 birth/turnover models with M0, M1 or M2 timing | `network_update_model()`, `simulate_network_updates()`, `fit_network_updates()` |
| Efficient-information check for later feedback | `feedback_score()` at an interior fitted M1 model |
| CRAN release/edit mixtures, personnel covariates and paired migrations | `cran_update_data(..., term_specification="motif_free", mixture="two_component")` |
| Calendar timing, release types and initial/later topology | `fit_cran_timing()` with manuscript model aliases |
| Structural and joint CRAN forecasts | `simulate_cran_updates()`, `simulate_cran_joint()` |

The reference recovery model uses explicit undirected adjacency matrices and
anchored edge turnover, with the paper's Poisson births, candidate statistics
and two feedback coefficients. Its default initial graph is one edge. It is
suited to small examples and independent verification; the large CRAN workflow
has native implementations. General vertex deletion is not exposed by these
network simulators. The generic proposal law accepts caller-defined edits,
whose graph admissibility remains the caller's responsibility.

The addition-only `network` interfaces and experimental relational interfaces
remain available. Their documentation identifies their narrower scope. The
historical `cs_mode="size_conditional"` draws a final update of a chosen size;
it is distinct from the current Poisson edit-set proposal construction.

## Verification

- Independently enumerated all proposal selections for overlapping groups and
  compared their union probabilities; checked normalization, rare-event
  numerical stability and simulated frequencies.
- Compared CS change statistics with direct graph recounting on every simple
  undirected graph with four vertices. Checked pre-event M2 intensities and
  interval compensators against closed-form turnover examples.
- Checked nesting of the reference M0/M1/M2 likelihoods, replay of simulated
  paths, fitting and failure on inadmissible data or computation guards.
- Checked CRAN mixture normalization, likelihood/sampler agreement, gradients,
  motif detection and shared-person effects. The public motif-free constructors
  reproduce the saved primary CS, CS-pair and TERGM complete log likelihoods;
  the timing interface reproduces all eight saved full-window timing fits to
  floating-point precision. These external-data comparisons supplement the
  self-contained package tests.
- Corrected inhomogeneous thinning bounds, temporal KDE normalization and tied
  event handling, and an experimental dyadic score recursion and amount law.
  Regression tests cover these changes.
- Added checks for native array dimensions, indices and numerical cache files,
  regenerated help, and exercised installed examples and the package check.

The package currently builds on macOS and Linux. CRAN archive data and prepared
caches are external inputs; the package does not silently fetch them. Numerical
bounds, optimizer convergence and boundary flags remain visible. Wald and
score approximations do not establish general asymptotic validity on growing
networks, and point-estimate simulations do not include parameter uncertainty.
Raw archive reconstruction, personnel identity resolution and the application's
four-week block-score effect intervals remain in the separate research scripts.
The CRAN guide describes the prepared-data contract. Mixture fits report point
estimates and do not claim to supply those robust intervals automatically.

Run `R CMD build hawkesNet` and `R CMD check --no-manual hawkesNet_0.2.0.tar.gz`
from the checkout's parent directory. The GitHub Actions workflow performs the
package check on Linux and macOS once activated. The configuration is supplied
as `inst/ci/R-CMD-check.yaml`: the publishing login lacked GitHub's `workflow`
permission, so this audit did not activate it. See `README.md`, installed examples and the
CRAN workflow guide for smaller entry points.
