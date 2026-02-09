# Minimal example: Barabási–Albert (BA) mark model
# Run from package root: source("inst/examples/example_BA.R")
# Or in R: setwd("/path/to/hawkesGrowthNet"); source("inst/examples/example_BA.R")

library(hawkesGrowthNet)

# Short time window so this finishes in under a minute
TIME <- 50
params <- list(mu = 10, beta_overall = 1, K = 0.5, beta_edges = 1, m = 1.5)

cat("Simulating one BA network...\n")
set.seed(1)
sim <- sim_hawkesGrowthNet(
  params = params,
  time_window = c(0, TIME),
  PMF_mark = PMF_mark_BA,
  cond_intensity = cond_intensity,
  hashed_edges = TRUE,
  verbose = TRUE,
  mu_multiplier = 3,
  truncation = 100
)

cat("Number of events:", length(sim$events$t), "\n")
cat("Network size (vertices):", network::network.size(sim$net), "\n")

# Fit the model to the simulated network
cat("Fitting BA model...\n")
params_init <- list(mu = 0.1, beta_overall = 0.5, beta_edges = 0.5, K = 0.5, m = 0.8)
fit <- tryCatch(
  fit_hawkesGrowthNet(
    params_init = params_init,
    time_window = c(0, TIME),
    mark_filtration = sim$net,
    PMF_mark = PMF_mark_BA,
    maxit = 500,
    truncation = 100,
    verbose = FALSE
  ),
  error = function(e) {
    cat("Fit failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(fit) && length(fit$fit) > 0) {
  cat("\nTrue vs fitted parameters:\n")
  print(data.frame(
    param = names(fit$fit$par),
    true = unlist(params)[names(fit$fit$par)],
    fitted = fit$fit$par
  ))
} else {
  cat("(Fit did not converge or failed; try increasing maxit or checking the network.)\n")
}
cat("\nDone. Try example_CS.R for the change-statistic (CS) model.\n")

