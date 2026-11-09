# Small addition-only preferential-attachment example.
# Run with source(system.file("examples", "example_BA.R", package = "hawkesNet"))
library(hawkesNet)

set.seed(1)
params <- list(mu = 4, K = 0.25, beta_overall = 1.2,
               beta_edges = 0.5, m = 1.4)
window <- c(0, 5)
sim <- sim_hawkesNet(params, window, PMF_mark_BA, cond_intensity,
                     truncation = 12L, verbose = FALSE)

# Estimate only the proposal mean m. This small run demonstrates the interface;
# it is not evidence of recovery or a substitute for fitting all unknown terms.
initial <- params
initial$m <- 0.9
fit <- fit_hawkesNet(
  initial, window, sim$net, PMF_mark_BA,
  truncation = 12L, fixed_params = setdiff(names(params), "m"),
  method = "L-BFGS-B", maxit = 100, get_hessian = FALSE, verbose = FALSE
)
print(fit$fit_table)
print(c(events = sim$events$n, truth_m = params$m,
        estimate_m = fit$params$m, convergence = fit$fit$convergence))
