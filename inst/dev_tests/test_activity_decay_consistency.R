## =============================================================================
## Activity-decay consistency test
## =============================================================================
## Tests that sim_hawkesNet -> fit_hawkesNet round-trips correctly when using
## mark_decay = "activity".  The CS simulation study uses "node_entrance" and
## works fine; this test isolates "activity" to check for inconsistencies.
##
## Target: < 2 minutes on a laptop (maxit kept low, small T).
## =============================================================================

library(hawkesNet)
library(network)
library(sna)

set.seed(42)
cat("=== Activity-decay consistency test ===\n\n")

# -----------------------------------------------------------------------------
# 1. Ground-truth parameters
# -----------------------------------------------------------------------------
params_true <- list(
  mu           = 10,
  beta_overall = 2,
  K            = 0.5,
  beta_edges   = 1,
  node_lambda  = 1,
  CS_params    = c(-6.7, 2, 0.1, -0.1)
)
FORMULA_RHS <- "edges + triangles + star(c(2,3))"
TRUNCATION  <- 100
T_SIM       <- 5        # short window for speed
MARK_DECAY  <- "activity"

cat("Ground truth:\n")
str(params_true)
cat("mark_decay:", MARK_DECAY, "\n")
cat("T:", T_SIM, "| truncation:", TRUNCATION, "\n\n")

# -----------------------------------------------------------------------------
# 2. Simulate from ground truth with mark_decay = "activity"
# -----------------------------------------------------------------------------
cat("--- Simulating from ground truth ---\n")
t0 <- proc.time()
sim <- sim_hawkesNet(
  params         = params_true,
  time_window    = c(0, T_SIM),
  PMF_mark       = PMF_mark_CS,
  cond_intensity = cond_intensity,
  formula_RHS    = FORMULA_RHS,
  truncation     = TRUNCATION,
  mark_decay     = MARK_DECAY,
  growth_only    = FALSE,
  max_node_time  = Inf,
  hashed_edges   = TRUE,
  verbose        = 2,
  mu_multiplier  = 5
)
t_sim <- (proc.time() - t0)[3]
cat(sprintf("Simulation: %d events, %d nodes, %d edges (%.1f s)\n",
            length(sim$events$t), network.size(sim$net),
            network.edgecount(sim$net), t_sim))

stopifnot(length(sim$events$t) > 20)
stopifnot(network.edgecount(sim$net) > 5)

n_events <- length(sim$events$t)
n_nodes  <- network.size(sim$net)
n_edges  <- network.edgecount(sim$net)
cat(sprintf("Edges per event: %.2f\n", n_edges / n_events))
cat(sprintf("Nodes per event: %.2f\n\n", n_nodes / n_events))

# -----------------------------------------------------------------------------
# 3. Evaluate log-likelihood at true params
# -----------------------------------------------------------------------------
cat("--- Log-likelihood at true params ---\n")
ll_true <- loglik_hawkesNet(
  params           = params_true,
  time_window      = c(0, T_SIM),
  mark_filtration  = sim$net,
  PMF_mark         = PMF_mark_CS,
  formula_RHS      = FORMULA_RHS,
  truncation       = TRUNCATION,
  mark_decay       = MARK_DECAY,
  growth_only      = FALSE,
  max_node_time    = max(get_times(sim$net)$node_times),
  cores            = 1
)
cat("  loglik(true):", ll_true$loglik, "\n")
stopifnot(is.finite(ll_true$loglik))

# Compare with node_entrance on same data
ll_ne <- loglik_hawkesNet(
  params           = params_true,
  time_window      = c(0, T_SIM),
  mark_filtration  = sim$net,
  PMF_mark         = PMF_mark_CS,
  formula_RHS      = FORMULA_RHS,
  truncation       = TRUNCATION,
  mark_decay       = "node_entrance",
  growth_only      = FALSE,
  max_node_time    = max(get_times(sim$net)$node_times),
  cores            = 1
)
cat("  loglik(true, node_entrance):", ll_ne$loglik, "\n")
cat("  Difference (activity - node_entrance):", ll_true$loglik - ll_ne$loglik, "\n\n")

# -----------------------------------------------------------------------------
# 4. Fit with mark_decay = "activity" (short run)
# -----------------------------------------------------------------------------
cat("--- Fitting with mark_decay = 'activity' (maxit=200) ---\n")
params_init <- list(
  mu           = 10,
  beta_overall = 1,
  K            = 0.5,
  beta_edges   = 0.5,
  node_lambda  = 0.5,
  CS_params    = c(-7, 1, 0, 0)
)
p_scale <- c(
  beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.1,
  CS_params1 = 1, CS_params2 = 0.1, CS_params3 = 0.1, CS_params4 = 0.1
)

t0 <- proc.time()
fit_activity <- fit_hawkesNet(
  params_init      = params_init,
  time_window      = c(0, T_SIM),
  mark_filtration  = sim$net,
  PMF_mark         = PMF_mark_CS,
  formula_RHS      = FORMULA_RHS,
  truncation       = TRUNCATION,
  mark_decay       = MARK_DECAY,
  growth_only      = FALSE,
  max_node_time    = max(get_times(sim$net)$node_times),
  method           = "Nelder-Mead",
  maxit            = 200,
  verbose          = TRUE,
  fixed_params     = c("K", "mu"),
  parscale         = p_scale,
  cores            = 1,
  cache_intensity  = TRUE,
  combine_intensity = FALSE
)
t_fit <- (proc.time() - t0)[3]
cat(sprintf("\nFit completed in %.1f s\n", t_fit))
cat("Convergence:", fit_activity$fit$convergence, "\n")
cat("loglik(fitted):", -fit_activity$fit$value, "\n")
cat("loglik(true):  ", ll_true$loglik, "\n\n")

# Check fitted params are in the ballpark
pfit <- fit_activity$params
cat("--- Parameter comparison ---\n")
cat(sprintf("  %-15s  %8s  %8s\n", "Parameter", "True", "Fitted"))
cat(sprintf("  %-15s  %8.3f  %8.3f\n", "beta_overall", params_true$beta_overall, pfit$beta_overall))
cat(sprintf("  %-15s  %8.3f  %8.3f\n", "beta_edges",   params_true$beta_edges,   pfit$beta_edges))
cat(sprintf("  %-15s  %8.3f  %8.3f\n", "node_lambda",  params_true$node_lambda,  pfit$node_lambda))
for (i in seq_along(params_true$CS_params)) {
  cat(sprintf("  %-15s  %8.3f  %8.3f\n", paste0("CS_params[", i, "]"),
              params_true$CS_params[i], pfit$CS_params[i]))
}

# -----------------------------------------------------------------------------
# 5. Simulate from FITTED params and compare
# -----------------------------------------------------------------------------
cat("\n--- Simulating from FITTED params ---\n")
N_RESIM <- 3
for (r in seq_len(N_RESIM)) {
  sim_fitted <- tryCatch({
    sim_hawkesNet(
      params         = pfit,
      time_window    = c(0, T_SIM),
      PMF_mark       = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS    = FORMULA_RHS,
      truncation     = TRUNCATION,
      mark_decay     = MARK_DECAY,
      growth_only    = FALSE,
      max_node_time  = Inf,
      hashed_edges   = TRUE,
      verbose        = FALSE,
      mu_multiplier  = 5
    )
  }, error = function(e) { cat("  Sim", r, "FAILED:", e$message, "\n"); NULL })

  if (!is.null(sim_fitted) && !is.null(sim_fitted$net)) {
    n_ev <- length(sim_fitted$events$t)
    n_nd <- network.size(sim_fitted$net)
    n_ed <- network.edgecount(sim_fitted$net)
    cat(sprintf("  Resim %d: %d events, %d nodes, %d edges (%.2f edges/event)\n",
                r, n_ev, n_nd, n_ed, n_ed / max(n_ev, 1)))
  }
}
cat(sprintf("\n  Original sim: %d events, %d nodes, %d edges (%.2f edges/event)\n",
            n_events, n_nodes, n_edges, n_edges / n_events))

# -----------------------------------------------------------------------------
# 6. Also simulate from TRUE params again for comparison
# -----------------------------------------------------------------------------
cat("\n--- Re-simulating from TRUE params ---\n")
for (r in seq_len(N_RESIM)) {
  sim_true2 <- tryCatch({
    sim_hawkesNet(
      params         = params_true,
      time_window    = c(0, T_SIM),
      PMF_mark       = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS    = FORMULA_RHS,
      truncation     = TRUNCATION,
      mark_decay     = MARK_DECAY,
      growth_only    = FALSE,
      max_node_time  = Inf,
      hashed_edges   = TRUE,
      verbose        = FALSE,
      mu_multiplier  = 5
    )
  }, error = function(e) { cat("  Sim", r, "FAILED:", e$message, "\n"); NULL })

  if (!is.null(sim_true2) && !is.null(sim_true2$net)) {
    n_ev <- length(sim_true2$events$t)
    n_nd <- network.size(sim_true2$net)
    n_ed <- network.edgecount(sim_true2$net)
    cat(sprintf("  True-resim %d: %d events, %d nodes, %d edges (%.2f edges/event)\n",
                r, n_ev, n_nd, n_ed, n_ed / max(n_ev, 1)))
  }
}

# -----------------------------------------------------------------------------
# 7. Compare with mark_decay = "node_entrance" fit on same data
# -----------------------------------------------------------------------------
cat("\n--- Fitting with mark_decay = 'node_entrance' for comparison (maxit=200) ---\n")
t0 <- proc.time()
fit_ne <- fit_hawkesNet(
  params_init      = params_init,
  time_window      = c(0, T_SIM),
  mark_filtration  = sim$net,
  PMF_mark         = PMF_mark_CS,
  formula_RHS      = FORMULA_RHS,
  truncation       = TRUNCATION,
  mark_decay       = "node_entrance",
  growth_only      = FALSE,
  max_node_time    = max(get_times(sim$net)$node_times),
  method           = "Nelder-Mead",
  maxit            = 200,
  verbose          = FALSE,
  fixed_params     = c("K", "mu"),
  parscale         = p_scale,
  cores            = 1,
  cache_intensity  = TRUE,
  combine_intensity = FALSE
)
t_fit_ne <- (proc.time() - t0)[3]
pfit_ne <- fit_ne$params
cat(sprintf("Fit (node_entrance) completed in %.1f s\n", t_fit_ne))
cat("loglik(fitted, node_entrance):", -fit_ne$fit$value, "\n\n")

cat("--- Parameter comparison: activity vs node_entrance fits ---\n")
cat(sprintf("  %-15s  %8s  %10s  %12s\n", "Parameter", "True", "Fit(activ)", "Fit(node_ent)"))
cat(sprintf("  %-15s  %8.3f  %10.3f  %12.3f\n", "beta_overall", params_true$beta_overall, pfit$beta_overall, pfit_ne$beta_overall))
cat(sprintf("  %-15s  %8.3f  %10.3f  %12.3f\n", "beta_edges",   params_true$beta_edges,   pfit$beta_edges,   pfit_ne$beta_edges))
cat(sprintf("  %-15s  %8.3f  %10.3f  %12.3f\n", "node_lambda",  params_true$node_lambda,  pfit$node_lambda,  pfit_ne$node_lambda))
for (i in seq_along(params_true$CS_params)) {
  cat(sprintf("  %-15s  %8.3f  %10.3f  %12.3f\n", paste0("CS_params[", i, "]"),
              params_true$CS_params[i], pfit$CS_params[i], pfit_ne$CS_params[i]))
}

cat("\n=== Activity-decay consistency test complete ===\n")
