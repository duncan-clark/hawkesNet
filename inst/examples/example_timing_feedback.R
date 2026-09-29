# Run after installing the updated package, or first use
# pkgload::load_all("hawkesNet") from the repository root.
# This small demonstration checks the workflow; it is not a recovery study.
library(hawkesNet)

params <- list(mu = 2, K = .2, beta_overall = 1.3, beta_edges = .3,
               m = 2, feedback_gamma = .7)
window <- c(0, 8)
timing_models <- list(
  M0 = hawkes_timing("M0"),
  M1 = hawkes_timing("M1", feature = "degree", scale = 1),
  M2 = hawkes_timing("M2", feature = "degree", scale = 1)
)
# Alternative feasible alpha: 2 * plogis(gamma * triangle participation).
triangle_feedback <- hawkes_timing("M2", feature = "triangles")

simulations <- lapply(seq_along(timing_models), function(i) {
  set.seed(813)
  sim_hawkesNet(params, window, PMF_mark_BA, cond_intensity,
               timing = timing_models[[i]], verbose = FALSE)
})
names(simulations) <- names(timing_models)
print(data.frame(model = names(simulations),
                 events = vapply(simulations, function(x) x$events$n, integer(1))))

# Compare all three specifications on the SAME M2-generated marked history.
# Keep the BA PMF and its parameters identical. Gamma-only estimation makes
# this a small runnable example; practical work must assess joint uncertainty.
observed <- simulations$M2$net
initial <- params
initial$feedback_gamma <- .1
fixed <- setdiff(names(initial), "feedback_gamma")
fits <- lapply(timing_models[c("M1", "M2")], function(timing) {
  fit_hawkesNet(initial, window, observed, PMF_mark_BA, timing = timing,
    fixed_params = fixed, method = "L-BFGS-B", maxit = 80,
    cache_intensity = TRUE, get_hessian = FALSE, verbose = FALSE)
})
comparison <- data.frame(
  model = c("M0", "M1", "M2"),
  gamma = c(0, fits$M1$params$feedback_gamma, fits$M2$params$feedback_gamma),
  loglik = c(loglik_hawkesNet(params, window, observed, PMF_mark_BA,
                            timing = timing_models$M0)$loglik,
             fits$M1$fit$value, fits$M2$fit$value)
)
print(comparison)

# Also evaluate the second alpha choice on this same complete marked history.
print(data.frame(model = "M2", feature = "triangles", gamma = params$feedback_gamma,
  loglik = loglik_hawkesNet(params, window, observed, PMF_mark_BA,
                          timing = triangle_feedback)$loglik))

compensator <- compensators_hawkesNet(fits$M2$params, window, observed,
                                     timing = timing_models$M2)
rescaled_waits <- diff(c(0, compensator))
print(data.frame(time = get_times(observed)$times, rescaled_wait = rescaled_waits))
# A fitted-sample KS p-value is descriptive; use simulation/refitting or held-out
# calibration for an inferential goodness-of-fit assessment.
