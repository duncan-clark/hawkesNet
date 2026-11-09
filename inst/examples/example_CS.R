# Small CS-1 example: frozen exponential weights and collapsed Poisson proposals.
# Run with source(system.file("examples", "example_CS.R", package = "hawkesNet"))
library(hawkesNet)

set.seed(1)
params <- list(mu = 4, K = 0.25, beta_overall = 1.2,
               beta_edges = 0.5, m = 1.4, node_lambda = 0.8,
               CS_params = c(0, 0.6))
window <- c(0, 5)
sim <- sim_hawkesNet(
  params, window, PMF_mark_CS, cond_intensity,
  formula_RHS = "edges + triangles", cs_mode = "independent",
  truncation = 4L, verbose = FALSE
)

# The edge-count coefficient cancels from normalized exponential weights.
# Here it and all parameters except m are fixed for a small interface example.
initial <- params
initial$m <- 0.9
fit <- fit_hawkesNet(
  initial, window, sim$net, PMF_mark_CS,
  formula_RHS = "edges + triangles", cs_mode = "independent",
  truncation = 4L, fixed_params = setdiff(names(params), "m"),
  method = "L-BFGS-B", maxit = 100, get_hessian = FALSE, verbose = FALSE
)
print(fit$fit_table)
print(c(events = sim$events$n, truth_m = params$m,
        estimate_m = fit$params$m, convergence = fit$fit$convergence))
