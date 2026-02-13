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
GROWTH_ONLY      <- TRUE
time_window      <- c(0, 1)

# --- 2. Pick which model to inspect ---
# Change this to inspect a different model
MODEL <- "structural"  # Options: "structural", "nodematch", "nodemix"
SEED_EVENTS <- 20      # Set to 0 for cold-start, or e.g. 20 for conditional simulation

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

# --- 4. Reconstruct parameters (name-based, matches gof.R logic) ---
cat("\n--- Parameter reconstruction (name-based) ---\n")
fixed <- fit$fixed_params
cat("Fixed params:", if (is.null(fixed)) "NONE (possible old RDS?)" else paste(fixed, collapse = ", "), "\n")

par_vec <- fit$fit$par
par_names <- names(par_vec)
cat("\nfit$par (", length(par_vec), "values):\n")
for (i in seq_along(par_vec)) {
  cat(sprintf("  [%d] %s = %.8f\n", i, par_names[i], par_vec[i]))
}

# Start from params_init and overwrite with fitted values by name
pfit <- params_init
pfit$vertex_categorical_levels <- params_init$vertex_categorical_levels

# Map scalar parameters
scalar_names <- c("mu", "beta_overall", "K", "beta_edges", "node_lambda", "m")
for (nm in scalar_names) {
  if (nm %in% par_names) {
    cat(sprintf("  Mapping %s: %.8f (from fit$par) -> overwriting params_init %.8f\n",
                nm, par_vec[nm], if (!is.null(pfit[[nm]])) pfit[[nm]] else NA))
    pfit[[nm]] <- par_vec[nm]
  } else {
    cat(sprintf("  Keeping %s: %.8f (from params_init, %s)\n",
                nm, if (!is.null(pfit[[nm]])) pfit[[nm]] else NA,
                if (!is.null(fixed) && nm %in% fixed) "FIXED" else "not in fit$par"))
  }
}

# Map CS_params: look for CS_params1, CS_params2, ... in par_vec
cs_idx <- grep("^CS_params[0-9]+$", par_names)
if (length(cs_idx) > 0) {
  cs_nums <- as.integer(sub("^CS_params", "", par_names[cs_idx]))
  n_cs <- length(pfit$CS_params)
  cat(sprintf("\nMapping %d CS_params entries into %d slots:\n", length(cs_idx), n_cs))
  for (j in seq_along(cs_idx)) {
    k <- cs_nums[j]
    if (k >= 1L && k <= n_cs) {
      cat(sprintf("  CS_params[%d] = %.8f (was %.8f)\n",
                  k, par_vec[cs_idx[j]], pfit$CS_params[k]))
      pfit$CS_params[k] <- par_vec[cs_idx[j]]
    } else {
      cat(sprintf("  WARNING: CS_params%d out of range (n_cs=%d)\n", k, n_cs))
    }
  }
}

# Map vertex_categorical
vc_idx <- grep("^vertex_categorical\\.", par_names)
if (length(vc_idx) > 0 && !is.null(pfit$vertex_categorical)) {
  cat(sprintf("\nMapping %d vertex_categorical entries:\n", length(vc_idx)))
  for (j in vc_idx) {
    parts <- strsplit(par_names[j], "\\.")[[1]]
    if (length(parts) >= 3) {
      attr_name <- parts[2]
      level_name <- paste(parts[3:length(parts)], collapse = ".")
      cat(sprintf("  %s$%s = %.8f\n", attr_name, level_name, par_vec[j]))
      if (!is.null(pfit$vertex_categorical[[attr_name]])) {
        pfit$vertex_categorical[[attr_name]][level_name] <- par_vec[j]
      }
    }
  }
}

# Handle inhomogeneous background
use_inhom <- !is.null(inhom_bg) && !is.null(inhom_bg$mu_fit) && !is.null(inhom_bg$mu_fit$mu_fun)
if (use_inhom) {
  Tval <- time_window[2] - time_window[1]
  pfit$mu <- inhom_bg$integral_bg / Tval
  cat("\nUsing inhomogeneous background (avg mu =", round(pfit$mu, 4), ")\n")
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

# --- 5b. Diagnostic: edge intercept analysis ---
if (!is.null(pfit$CS_params) && length(pfit$CS_params) >= 1) {
  cs1 <- pfit$CS_params[1]
  p_baseline <- plogis(cs1)
  cat("\n=== Edges intercept diagnostic ===\n")
  cat(sprintf("  CS_params[1] (edges intercept) = %.4f\n", cs1))
  cat(sprintf("  Baseline edge probability (no structure) = plogis(%.4f) = %.2e\n", cs1, p_baseline))
  
  # What the observed edge density implies
  n_v <- network.size(net_obs)
  n_e <- network.edgecount(net_obs)
  obs_dens <- if (n_v > 1) 2 * n_e / (n_v * (n_v - 1)) else NA
  if (!is.na(obs_dens)) {
    obs_logit <- qlogis(obs_dens)
    cat(sprintf("  Observed edge density = %.6f  (logit = %.2f)\n", obs_dens, obs_logit))
    cat(sprintf("  Gap: edges intercept (%.2f) vs logit(obs density) (%.2f) = %.1f\n",
                cs1, obs_logit, cs1 - obs_logit))
    if (cs1 < obs_logit - 4) {
      cat("  *** WARNING: Edges intercept is much more negative than logit(density).\n")
      cat("  *** The model relies heavily on triangles/gwdegree to explain edges.\n")
      cat("  *** Cold-start simulation will produce very few edges.\n")
      cat("  *** Consider re-fitting with cs_intercept_floor = -8 and\n")
      cat("  *** initializing CS_params[1] at logit(edge_density).\n")
    }
  }
}

# --- 6. Compare with observed network ---
cat("\n=== Observed network ===\n")
cat("  Nodes:", network.size(net_obs), "\n")
cat("  Edges:", network.edgecount(net_obs), "\n")
cat("  Mean degree:", round(mean(degree(net_obs, gmode = "graph")), 2), "\n")
cat("  Events:", length(get_times(net_obs)$times), "\n")

# --- 7. Run TWO verbose simulations for comparison ---
cat("\n=== Running Comparison Simulations ===\n")
if (SEED_EVENTS > 0) {
  cat(sprintf("  (Conditional on first %d events)\n", SEED_EVENTS))
}
cat("  (This may take a few minutes for large networks)\n\n")

# Extract seed if requested
seed_net <- NULL
seed_times <- NULL
if (SEED_EVENTS > 0) {
  all_times <- get_times(net_obs)$times
    if (length(all_times) >= SEED_EVENTS) {
      t_seed <- all_times[SEED_EVENTS]
      seed_net <- filtration_to_net(net_obs, t_seed, equals = TRUE)
      seed_times <- all_times[1:SEED_EVENTS]
      cat(sprintf("Seeding with first %d events (up to t=%.4f)\n", SEED_EVENTS, t_seed))
      cat(sprintf("Seed network size: %d nodes, %d edges\n", 
                  network.size(seed_net), network.edgecount(seed_net)))
    }
}

# --- Simulation 1: Growth-Only (New Logic) ---
cat("\n>>> SIMULATION 1: Growth-Only (Citation Logic) <<<\n")
t_sim1 <- proc.time()
sim_growth <- tryCatch({
  sim_hawkesNet(
    params = pfit,
    time_window = time_window,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = formula_RHS,
    truncation = TRUNCATION,
    mark_decay = "activity",
    growth_only = TRUE,
    max_node_time = 1,
    hashed_edges = TRUE,
    verbose = TRUE,
    mu_multiplier = 5,
    stop_on_full_network = FALSE,
    inhom_bg = inhom_bg,
    seed_net = seed_net,
    seed_times = seed_times
  )
}, error = function(e) {
  cat("\n*** GROWTH-ONLY SIMULATION FAILED:", e$message, "***\n")
  NULL
})
elapsed_sim1 <- (proc.time() - t_sim1)[3]

# --- Simulation 2: All-Edges (Old Logic) ---
cat("\n>>> SIMULATION 2: All-Edges (Continuous Collaboration Logic) <<<\n")
t_sim2 <- proc.time()
sim_all <- tryCatch({
  sim_hawkesNet(
    params = pfit,
    time_window = time_window,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = formula_RHS,
    truncation = TRUNCATION,
    mark_decay = "activity",
    growth_only = FALSE,
    max_node_time = 1,
    hashed_edges = TRUE,
    verbose = TRUE,
    mu_multiplier = 5,
    stop_on_full_network = FALSE,
    inhom_bg = inhom_bg,
    seed_net = seed_net,
    seed_times = seed_times
  )
}, error = function(e) {
  cat("\n*** ALL-EDGES SIMULATION FAILED:", e$message, "***\n")
  NULL
})
elapsed_sim2 <- (proc.time() - t_sim2)[3]

# --- 8. Compare simulated vs observed ---
cat("\n=== FINAL COMPARISON ===\n")
cat(sprintf("  %-20s %10s %15s %15s\n", "", "Observed", "Growth-Only", "All-Edges"))

n_obs <- network.size(net_obs)
e_obs <- network.edgecount(net_obs)
d_obs <- mean(degree(net_obs, gmode = "graph"))
ev_obs <- length(get_times(net_obs)$times)

n_growth <- if(!is.null(sim_growth)) network.size(sim_growth$net) else NA
e_growth <- if(!is.null(sim_growth)) network.edgecount(sim_growth$net) else NA
d_growth <- if(!is.null(sim_growth)) mean(degree(sim_growth$net, gmode = "graph")) else NA
ev_growth <- if(!is.null(sim_growth)) length(sim_growth$events$t) else NA

n_all <- if(!is.null(sim_all)) network.size(sim_all$net) else NA
e_all <- if(!is.null(sim_all)) network.edgecount(sim_all$net) else NA
d_all <- if(!is.null(sim_all)) mean(degree(sim_all$net, gmode = "graph")) else NA
ev_all <- if(!is.null(sim_all)) length(sim_all$events$t) else NA

cat(sprintf("  %-20s %10d %15d %15d\n", "Nodes", n_obs, n_growth, n_all))
cat(sprintf("  %-20s %10d %15d %15d\n", "Edges", e_obs, e_growth, e_all))
cat(sprintf("  %-20s %10.2f %15.2f %15.2f\n", "Mean degree", d_obs, d_growth, d_all))
cat(sprintf("  %-20s %10d %15d %15d\n", "Events", ev_obs, ev_growth, ev_all))
cat(sprintf("  %-20s %10s %15.1fs %15.1fs\n", "Sim Time", "", elapsed_sim1, elapsed_sim2))

cat("\n=== Interpretation ===\n")
cat("Note: Both simulations use the SAME fitted parameters (which were likely\n")
cat("fitted under the old 'All-Edges' logic). You should see that 'All-Edges'\n")
cat("produces much thinner networks because the probability is diluted over\n")
cat("thousands of candidates, whereas 'Growth-Only' focuses probability on\n")
cat("new entrants.\n")

cat("\n=== Manual GOF script complete ===\n")
cat("Objects available for inspection: pfit, net_obs, sim_growth, sim_all, dat\n")
