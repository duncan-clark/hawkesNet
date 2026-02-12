# =============================================================================
# Manual GOF diagnostic script for OpenAlex results
# =============================================================================
# Run interactively from package root:
#   source("inst/openalex_study/manual_gof.R")
#
# This script loads saved results, reconstructs the fitted parameters exactly
# as gof() would, runs ONE verbose simulation, and lets you poke around.
# =============================================================================

library(hawkesNet)
library(network)
library(sna)
library(ernm)

# --- 1. Load saved results ---
PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) Sys.getenv("SLURM_SUBMIT_DIR") else getwd()
rds_path <- file.path(PKG_ROOT, "cluster_output", "results_openalex_full.RDS")
if (!file.exists(rds_path)) {
  rds_path <- "cluster_output/results_openalex_full.RDS"
}
stopifnot(file.exists(rds_path))
cat("Loading results from:", rds_path, "\n")
dat <- readRDS(rds_path)

# Extract key objects
net_obs          <- dat$net_raw
inhom_bg         <- dat$inhom_bg
TRUNCATION       <- 300L
time_window      <- c(0, 1)

# --- 2. Pick which model to inspect ---
# Change this to inspect a different model
MODEL <- "structural"  # Options: "structural", "nodematch", "nodemix"

if (MODEL == "structural") {
  fit         <- dat$fit_inhom_structural
  params_init <- dat$params_init_structural
  formula_RHS <- dat$FORMULA_RHS_STRUCTURAL
} else if (MODEL == "nodematch") {
  fit         <- dat$fit_inhom_nodematch
  params_init <- dat$params_init_nodematch
  formula_RHS <- dat$FORMULA_RHS_NODEMATCH
} else if (MODEL == "nodemix") {
  fit         <- dat$fit_inhom_nodemix
  params_init <- dat$params_init_nodemix
  formula_RHS <- dat$FORMULA_RHS_NODEMIX
} else {
  stop("Unknown MODEL: ", MODEL)
}

cat("\n=== Manual GOF for:", MODEL, "model ===\n")
cat("Formula:", formula_RHS, "\n")

if (is.null(fit)) stop("No fit available for model: ", MODEL)

# --- 3. Print the fit table ---
cat("\nFit table:\n")
if (!is.null(fit$fit_table)) {
  print(fit$fit_table, row.names = FALSE)
} else {
  cat("  (no fit_table; raw par:)\n")
  print(fit$fit$par)
}
cat("\nConvergence:", fit$fit$convergence, "\n")
cat("Iterations:", fit$fit$counts[1], "\n")

# --- 4. Reconstruct parameters exactly as gof() does ---
cat("\n--- Parameter reconstruction ---\n")
fixed <- fit$fixed_params
cat("Fixed params:", if (is.null(fixed)) "NONE (BUG?)" else paste(fixed, collapse = ", "), "\n")

# Build skeleton (strip fixed params + vertex_categorical_levels)
skel <- params_init
if (!is.null(fixed)) {
  for (p in fixed) skel[[p]] <- NULL
}
skel$vertex_categorical_levels <- NULL

cat("\nfit$par length:", length(fit$fit$par), "\n")
cat("skeleton length:", length(unlist(skel)), "\n")
cat("fit$par names:", paste(names(fit$fit$par), collapse = ", "), "\n")
cat("skeleton names:", paste(names(unlist(skel)), collapse = ", "), "\n")

# Relist
pfit_vals <- tryCatch({
  relist(fit$fit$par, skeleton = skel)
}, error = function(e) {
  cat("\n*** RELIST FAILED:", e$message, "***\n")
  cat("This means GOF would fall back to params_init!\n")
  NULL
})

# Merge into full params
pfit <- params_init
if (!is.null(pfit_vals)) {
  for (n in names(pfit_vals)) pfit[[n]] <- pfit_vals[[n]]
  cat("\nRelist succeeded. Merged fitted values into params_init.\n")
} else {
  cat("\nWARNING: Using params_init as fallback!\n")
}

# Restore metadata
pfit$vertex_categorical_levels <- params_init$vertex_categorical_levels

# Handle inhomogeneous background
use_inhom <- !is.null(inhom_bg) && !is.null(inhom_bg$mu_fit) && !is.null(inhom_bg$mu_fit$mu_fun)
if (use_inhom) {
  Tval <- time_window[2] - time_window[1]
  pfit$mu <- inhom_bg$integral_bg / Tval
  cat("Using inhomogeneous background (avg mu =", round(pfit$mu, 4), ")\n")
} 

# Clamp to valid ranges
pfit$K <- min(max(pfit$K, 0.001), 0.999)
pfit$mu <- max(pfit$mu, 0.001)
pfit$node_lambda <- max(pfit$node_lambda, 0.1)
pfit$beta_overall <- max(pfit$beta_overall, 0.001)
pfit$beta_edges <- max(pfit$beta_edges, 0.001)

# Repair vertex_categorical if present
if (!is.null(pfit$vertex_categorical)) {
  pfit <- repair_vertex_categorical_params(pfit, eps = 1e-6)
}

# --- 5. Print final simulation parameters ---
cat("\n=== Parameters for GOF simulation ===\n")
scalar_params <- c("mu", "beta_overall", "K", "beta_edges", "node_lambda")
for (p in scalar_params) {
  if (!is.null(pfit[[p]])) {
    src <- if (!is.null(fixed) && p %in% fixed) "(FIXED)" else "(fitted)"
    cat(sprintf("  %-15s = %12.6f  %s\n", p, pfit[[p]], src))
  }
}
if (!is.null(pfit$CS_params)) {
  cat("  CS_params      =", paste(round(pfit$CS_params, 6), collapse = ", "), " (fitted)\n")
}
if (!is.null(pfit$vertex_categorical)) {
  for (attr_name in names(pfit$vertex_categorical)) {
    cat("  vertex_categorical$", attr_name, " =", 
        paste(round(pfit$vertex_categorical[[attr_name]], 4), collapse = ", "), "\n")
  }
}

# --- 6. Compare with observed network ---
cat("\n=== Observed network ===\n")
cat("  Nodes:", network.size(net_obs), "\n")
cat("  Edges:", network.edgecount(net_obs), "\n")
cat("  Mean degree:", round(mean(degree(net_obs, gmode = "graph")), 2), "\n")
cat("  Events:", length(get_times(net_obs)$times), "\n")

# --- 7. Run ONE verbose simulation ---
cat("\n=== Running 1 verbose simulation ===\n")
cat("  (This may take a few minutes for large networks)\n\n")
t_sim <- proc.time()

sim_result <- tryCatch({
  sim_hawkesNet(
    params = pfit,
    time_window = time_window,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = formula_RHS,
    truncation = TRUNCATION,
    mark_decay = "activity",
    max_node_time = 1,
    hashed_edges = TRUE,
    verbose = TRUE,
    mu_multiplier = 5,
    stop_on_full_network = FALSE,
    inhom_bg = inhom_bg
  )
}, error = function(e) {
  cat("\n*** SIMULATION FAILED:", e$message, "***\n")
  NULL
})

elapsed_sim <- (proc.time() - t_sim)[3]
cat("\nSimulation took:", round(elapsed_sim, 1), "s\n")

# --- 8. Compare simulated vs observed ---
if (!is.null(sim_result) && !is.null(sim_result$net)) {
  net_sim <- sim_result$net
  cat("\n=== Simulated network ===\n")
  cat("  Nodes:", network.size(net_sim), "\n")
  cat("  Edges:", network.edgecount(net_sim), "\n")
  cat("  Mean degree:", round(mean(degree(net_sim, gmode = "graph")), 2), "\n")
  cat("  Events:", length(sim_result$events$t), "\n")
  
  cat("\n=== Comparison ===\n")
  cat(sprintf("  %-20s %10s %10s\n", "", "Observed", "Simulated"))
  cat(sprintf("  %-20s %10d %10d\n", "Nodes", network.size(net_obs), network.size(net_sim)))
  cat(sprintf("  %-20s %10d %10d\n", "Edges", network.edgecount(net_obs), network.edgecount(net_sim)))
  cat(sprintf("  %-20s %10.2f %10.2f\n", "Mean degree", 
              mean(degree(net_obs, gmode = "graph")), mean(degree(net_sim, gmode = "graph"))))
  cat(sprintf("  %-20s %10d %10d\n", "Events", 
              length(get_times(net_obs)$times), length(sim_result$events$t)))
  
  # Ratio
  node_ratio <- network.size(net_sim) / network.size(net_obs)
  edge_ratio <- network.edgecount(net_sim) / network.edgecount(net_obs)
  cat(sprintf("\n  Node ratio (sim/obs): %.2f\n", node_ratio))
  cat(sprintf("  Edge ratio (sim/obs): %.2f\n", edge_ratio))
  
  if (node_ratio < 0.5 || node_ratio > 2.0) {
    cat("\n  *** WARNING: Simulated network size is very different from observed! ***\n")
    cat("  This suggests the parameters may not be correctly reconstructed.\n")
  } else {
    cat("\n  Network sizes are in reasonable agreement.\n")
  }
} else {
  cat("\n  No simulated network to compare.\n")
}

cat("\n=== Manual GOF script complete ===\n")
cat("Objects available for inspection: pfit, net_obs, sim_result, dat\n")
