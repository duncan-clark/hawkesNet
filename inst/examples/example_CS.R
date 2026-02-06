# Minimal example: Change-statistic (CS) mark model
# Run from package root: source("inst/examples/example_CS.R")
# Or in R: setwd("/path/to/hawkesGrowthNet"); source("inst/examples/example_CS.R")

library(hawkesGrowthNet)

# Short time window so this finishes in a few minutes
TIME <- 50
TRUNCATION <- 100
params <- list(
  mu = 10,
  beta_overall = 2,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-6.7, 2, 0.1, -0.1)
)
formula_RHS <- "edges + triangles() + star(c(2,3))"

cat("Simulating one CS network...\n")
set.seed(1)
sim <- sim_hawkesGrowthNet(
  params = params,
  time_window = c(0, TIME),
  PMF_mark = PMF_mark_CS,
  cond_intensity = cond_intensity,
  hashed_edges = TRUE,
  verbose = FALSE,
  mu_multiplier = 3,
  truncation = TRUNCATION,
  formula_RHS = formula_RHS
)

cat("Number of events:", length(sim$events$t), "\n")
cat("Network size (vertices):", network::network.size(sim$net), "\n")

# Fit the model (K fixed for identifiability)
cat("Fitting CS model...\n")
params_init <- list(
  mu = 10,
  beta_overall = 1,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-10, 0, 0, 0)
)
fit <- tryCatch(
  fit_hawkesGrowthNet(
    params_init = params_init,
    time_window = c(0, TIME),
    mark_filtration = sim$net,
    PMF_mark = PMF_mark_CS,
    formula_RHS = "edges + triangles + star(c(2,3))",
    maxit = 500,
    truncation = TRUNCATION,
    fixed_params = c("K"),
    verbose = FALSE
  ),
  error = function(e) {
    cat("Fit failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(fit) && length(fit$fit) > 0) {
  cat("\nTrue vs fitted parameters (K fixed):\n")
  print(data.frame(
    param = names(fit$fit$par),
    true = unlist(params)[names(fit$fit$par)],
    fitted = fit$fit$par
  ))
} else {
  cat("(Fit did not converge or failed; try increasing maxit or checking the network.)\n")
}
cat("\nDone. See inst/simulation_study/ for full studies.\n")
