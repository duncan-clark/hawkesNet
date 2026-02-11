# =============================================================================
# OpenAlex Interactive Script: Structural-only fit for interactive RStudio cluster
# =============================================================================
# This script runs only the structural fit (no gender terms) for quick testing
# Run interactively in RStudio on the cluster
# =============================================================================

# Load hawkesNet package
library(hawkesNet)

# Load required libraries
library(dplyr)
library(network)
library(sna)
library(ernm)
library(parallel)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
}

# Paths: run from package root (directory containing inst/)
PKG_ROOT <- getwd()
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) {
  parent_dir <- dirname(PKG_ROOT)
  if (file.exists(file.path(parent_dir, "DESCRIPTION"))) {
    PKG_ROOT <- parent_dir
  } else {
    if (file.exists(file.path(PKG_ROOT, "inst", "openalex_study", "openalex_hawkes_study.R"))) {
      # We're already in package root
    } else {
      warning("Could not find package root. Using current directory: ", PKG_ROOT)
    }
  }
}
cat("Package root:", PKG_ROOT, "\n")
cat("Working directory:", getwd(), "\n")
source(file.path(PKG_ROOT, "inst", "openalex_study", "get_network_openalex.R"))

# =============================================================================
# Config (override via env or edit)
# =============================================================================
EMAIL <- Sys.getenv("OPENALEX_EMAIL", "duncan-clark@outlook.com")
SEARCH_STRING <- Sys.getenv("OPENALEX_STRING", "Hawkes Process")
PAGES <- as.integer(Sys.getenv("OPENALEX_PAGES", 100))
PER_PAGE <- 100L
MIN_DATE <- "1971-04-01"
MAX_DATE <- "2020-01-01"
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 16))
MAX_ITER <- 5000
TRUNCATION <- 100L
GOF_TIME_WINDOW <- c(0, 1)  # Full time period for GOF simulations
N_GOF <- 25L   # number of simulated networks for goodness-of-fit
TOPIC <- "Point processes and geometric inequalities"

# =============================================================================
# 1. Fetch data and prepare network
# =============================================================================
t_total <- proc.time()
cat("=== OpenAlex Interactive Study (Structural-only) ===\n")
cat("  Search:", SEARCH_STRING, "| Pages:", PAGES, "| Cores:", N_CORES, "\n")
cat("  Date range:", MIN_DATE, "to", MAX_DATE, "\n")
cat("  Package root:", PKG_ROOT, "\n")
cat("  Working directory:", getwd(), "\n")

cat("\n--- Step 1: Fetch data and prepare network ---\n")
t_step <- proc.time()
out <- get_network(email = EMAIL, pages = PAGES, per_page = PER_PAGE,
                   string = SEARCH_STRING, min_date = MIN_DATE, max_date = MAX_DATE, topics_include = c(TOPIC))
net_raw <- out$net
edges <- out$edges
network::set.vertex.attribute(net_raw, "time", net_raw %v% "time_scaled")
network::set.edge.attribute(net_raw, "time", net_raw %e% "time_scaled")
net_raw <- hawkesNet::normalize_times_01(net_raw, attr = "time", keep_na = TRUE)
n_events <- length(hawkesNet::get_times(net_raw)$times)
n_nodes <- network::network.size(net_raw)
cat("  Network:", n_events, "events,", n_nodes, "nodes\n")
cat("  Step 1 took:", round((proc.time() - t_step)[3], 1), "s\n\n")

# ========================
# Inhom BG
# =======================

cat("--- Step 2: Structural-only fit (no gender) ---\n")
t_step <- proc.time()
time_window_01 <- c(0, 1)
cat("  Preparing inhomogeneous background (KDE)...\n")
t_kde <- proc.time()
inhom_bg <- tryCatch(
  prepare_inhomogeneous_background(net_raw, time_attr = "time", bw = NULL, grid_n = 2048),
  error = function(e) { cat("  ERROR: prepare_inhomogeneous_background failed:", e$message, "\n"); NULL }
)
cat("  KDE background:", round((proc.time() - t_kde)[3], 1), "s\n")

FORMULA_RHS_STRUCTURAL <- "edges + triangles + gwdegree(0.5)"

cat("\n  Formula:", FORMULA_RHS_STRUCTURAL, "\n")

# Get expected parameters for structural formula
exp_cs_structural <- expected_params_PMF_mark_CS(net_raw, FORMULA_RHS_STRUCTURAL)
n_cs_structural <- if (!is.na(exp_cs_structural$CS_params_length)) exp_cs_structural$CS_params_length else 4L



# Create parscale for structural model
p_scale_structural <- c(
  beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
  setNames(rep(0.1, n_cs_structural), paste0("CS_params", seq_len(n_cs_structural)))
)

cat("  Method: Nelder-Mead (max", MAX_ITER, "iterations)\n")

# Initialize parameters for structural-only model (no vertex_categorical)
params_init_structural <- list(
  mu = inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1]),
  beta_overall = 1,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-10, rep(0, n_cs_structural - 1))
)

# =============================================================================
# 2. Structural-only fit
# =============================================================================

if (!is.null(inhom_bg)) {

  t_fit_structural <- proc.time()
  
  fit_inhom_structural <- tryCatch(
    fit_hawkesNet_inhom(
      params_init = params_init_structural,
      time_window = time_window_01,
      mark_filtration = net_raw,
      PMF_mark = PMF_mark_CS,
      mu_vec = inhom_bg$mu_vec,
      integral_bg = inhom_bg$integral_bg,
      formula_RHS = FORMULA_RHS_STRUCTURAL,
      truncation = TRUNCATION,
      mark_decay = "activity",
      max_node_time = 1,
      method = "Nelder-Mead",
      maxit = MAX_ITER,
      trace = 1,
      reltol = 1e-8,
      verbose = TRUE,
      fixed_params = c("K", "mu"),
      parscale = p_scale_structural,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      cores = N_CORES
    ),
    error = function(e) { cat("  ERROR: Fit failed:", e$message, "\n"); NULL }
  )
  
  elapsed_fit_structural <- (proc.time() - t_fit_structural)[3]
  if (!is.null(fit_inhom_structural)) {
    cat("  Fit completed:", round(elapsed_fit_structural, 1), "s (", round(elapsed_fit_structural / 60, 1), "min)\n")
    cat("  Convergence:", fit_inhom_structural$fit$convergence, "\n")
    cat("  Iterations:", fit_inhom_structural$fit$counts[1], "\n")
    # Print fit table
    if (!is.null(fit_inhom_structural$fit_table)) {
      cat("\n  Structural-only fit results:\n")
      print(fit_inhom_structural$fit_table, max = NULL)
    }
  } else {
    cat("  Fit FAILED after", round(elapsed_fit_structural, 1), "s\n")
  }
  cat("  Step 2 total:", round((proc.time() - t_step)[3], 1), "s\n")
} else {
  cat("  No inhomogeneous background; skipping structural fit\n")
}

cat("\n=== Interactive study complete ===\n")
cat("  Total wall time:", round((proc.time() - t_total)[3], 1), "s (", round((proc.time() - t_total)[3] / 60, 1), "min)\n")


# =============================================================================
# 2b. nodeMatch fit (initialized from structural fit)
# =============================================================================

# --- between fit 1 and fit 2 ---
fit_inhom_structural$intens_funcs <- NULL  # strip heavy closures first
fit_inhom_structural <- NULL
rm(out)  # large OpenAlex API response
gc()     # CRITICAL: reclaim memory before pbmclapply forks 50 processes

fit_inhom_nodematch <- NULL
FORMULA_RHS_NODEMATCH <- "edges + triangles + gwdegree(0.5) + nodeMatch('gender')"
if (!is.null(inhom_bg)) {
  cat("\n--- Step 2b: nodeMatch fit (initialized from structural) ---\n")
  cat("  Formula:", FORMULA_RHS_NODEMATCH, "\n")
  t_step_nodematch <- proc.time()
  
  # Get expected parameters for nodeMatch formula (ensure length >= structural + 1 for nodeMatch term when NA)
  exp_cs_nodematch <- expected_params_PMF_mark_CS(net_raw, FORMULA_RHS_NODEMATCH)
  n_cs_nodematch <- if (!is.na(exp_cs_nodematch$CS_params_length)) {
    exp_cs_nodematch$CS_params_length
  } else {
    max(5L, (if (exists("n_cs_structural", inherits = FALSE)) n_cs_structural else 5L) + 1L)
  }
  
  # Initialize nodeMatch fit from structural fit results
  if (!is.null(fit_inhom_structural) && fit_inhom_structural$fit$convergence == 0) {
    # Extract structural parameters from structural fit
    # CRITICAL: Create skeleton that matches EXACTLY what was used during optimization
    skel_structural <- list(
      mu = params_init_structural$mu,
      beta_overall = params_init_structural$beta_overall,
      K = params_init_structural$K,
      beta_edges = params_init_structural$beta_edges,
      node_lambda = params_init_structural$node_lambda,
      CS_params = params_init_structural$CS_params
    )
    pfit_structural <- tryCatch({
      relist(fit_inhom_structural$fit$par, skeleton = skel_structural)
    }, error = function(e) {
      cat("  Warning: Failed to extract structural parameters:", e$message, "\n")
      cat("  Using independent initialization\n")
      NULL
    })
    
    if (!is.null(pfit_structural)) {
      # Validate extracted parameters - check each field individually
      mu_valid <- is.finite(pfit_structural$mu) && pfit_structural$mu > 0
      beta_overall_valid <- is.finite(pfit_structural$beta_overall)
      K_valid <- is.finite(pfit_structural$K)
      beta_edges_valid <- is.finite(pfit_structural$beta_edges)
      node_lambda_valid <- is.finite(pfit_structural$node_lambda) && pfit_structural$node_lambda > 0
      cs_valid <- all(is.finite(pfit_structural$CS_params)) && length(pfit_structural$CS_params) >= n_cs_structural
      
      if (mu_valid && beta_overall_valid && K_valid && beta_edges_valid && 
          node_lambda_valid && cs_valid) {
        # Initialize nodeMatch: use structural CS_params (edges, degree(0), triangles, stars), add nodeMatch term
        cs_structural <- pfit_structural$CS_params[seq_len(min(n_cs_structural, length(pfit_structural$CS_params)))]
        cs_padding <- rep(0, max(0L, n_cs_nodematch - length(cs_structural)))
        cs_init_nodematch <- c(cs_structural, cs_padding)[seq_len(n_cs_nodematch)]
        
        # Ensure mu is calculated correctly
        mu_val <- if (is.finite(pfit_structural$mu) && pfit_structural$mu > 0) {
          pfit_structural$mu
        } else {
          inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1])
        }
        
        params_init_nodematch <- list(
          mu = mu_val,
          beta_overall = pfit_structural$beta_overall,
          K = pfit_structural$K,
          beta_edges = pfit_structural$beta_edges,
          node_lambda = pfit_structural$node_lambda,
          CS_params = cs_init_nodematch,
          vertex_categorical = list(gender = c(female = 0.1, male = 0.5)),
          vertex_categorical_levels = list(gender = c("female", "male", "unknown"))
        )
        cat("  Initialized from structural fit\n")
      } else {
        pfit_structural <- NULL  # Force fallback
      }
    }
    
    if (is.null(pfit_structural)) {
      # Fallback: independent initialization
      params_init_nodematch <- list(
        mu = inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1]),
        beta_overall = 1,
        K = 0.5,
        beta_edges = 1,
        node_lambda = 1,
        CS_params = c(-10, rep(0, n_cs_nodematch - 1)),
        vertex_categorical = list(gender = c(female = 0.1, male = 0.5)),
        vertex_categorical_levels = list(gender = c("female", "male", "unknown"))
      )
      cat("  Structural fit parameters invalid; using independent initialization\n")
    }
  } else {
    # Fallback: independent initialization if structural fit failed
    params_init_nodematch <- list(
      mu = inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1]),
      beta_overall = 1,
      K = 0.5,
      beta_edges = 1,
      node_lambda = 1,
      CS_params = c(-10, rep(0, n_cs_nodematch - 1)),
      vertex_categorical = list(gender = c(female = 0.1, male = 0.5)),
      vertex_categorical_levels = list(gender = c("female", "male", "unknown"))
    )
    cat("  Structural fit not available; using independent initialization\n")
  }
  
  # Create parscale for nodeMatch
  p_scale_nodematch <- c(
    beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
    setNames(rep(0.1, n_cs_nodematch), paste0("CS_params", seq_len(n_cs_nodematch))),
    vertex_categorical.gender.female = 0.1, vertex_categorical.gender.male = 0.1
  )
  
  cat("  Method: Nelder-Mead (max", MAX_ITER, "iterations)\n")
  t_fit_nodematch <- proc.time()
  
  fit_inhom_nodematch <- tryCatch(
    fit_hawkesNet_inhom(
      params_init = params_init_nodematch,
      time_window = time_window_01,
      mark_filtration = net_raw,
      PMF_mark = PMF_mark_CS,
      mu_vec = inhom_bg$mu_vec,
      integral_bg = inhom_bg$integral_bg,
      formula_RHS = FORMULA_RHS_NODEMATCH,
      truncation = TRUNCATION,
      mark_decay = "activity",
      max_node_time = 1,
      method = "Nelder-Mead",
      maxit = MAX_ITER,
      trace = 1,
      reltol = 1e-8,
      verbose = FALSE,
      fixed_params = c("K", "mu"),
      parscale = p_scale_nodematch,
      cache_intensity = TRUE,
      cores = N_CORES,
      combine_intensity = TRUE
    ),
    error = function(e) { cat("  ERROR: Fit failed:", e$message, "\n"); NULL }
  )
  
  elapsed_fit_nodematch <- (proc.time() - t_fit_nodematch)[3]
  if (!is.null(fit_inhom_nodematch)) {
    cat("  Fit completed:", round(elapsed_fit_nodematch, 1), "s (", round(elapsed_fit_nodematch / 60, 1), "min)\n")
    cat("  Convergence:", fit_inhom_nodematch$fit$convergence, "\n")
    cat("  Iterations:", fit_inhom_nodematch$fit$counts[1], "\n")
    # Print fit table
    if (!is.null(fit_inhom_nodematch$fit_table)) {
      cat("\n  nodeMatch fit results:\n")
      print(fit_inhom_nodematch$fit_table, max = NULL)
    }
  } else {
    cat("  Fit FAILED after", round(elapsed_fit_nodematch, 1), "s\n")
  }
  cat("  Step 2b total:", round((proc.time() - t_step_nodematch)[3], 1), "s\n")
} else {
  cat("  No inhomogeneous background; skipping nodeMatch fit\n")
}

cat("\n  Step 2 total:", round((proc.time() - t_step)[3], 1), "s\n\n")


# =========================================
# GOF : 
# ========================================

cat("--- Step 4: Goodness-of-fit ---\n")
t_step <- proc.time()
GOF_results_structural <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                               geodist_obs = NULL, geodist_sim = NULL,
                               wait_obs = NULL, wait_sim = NULL)
GOF_results_nodematch <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                              geodist_obs = NULL, geodist_sim = NULL,
                              wait_obs = NULL, wait_sim = NULL)
GOF_results <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                    geodist_obs = NULL, geodist_sim = NULL,
                    wait_obs = NULL, wait_sim = NULL)

# GOF for structural-only model (first)
if (RUN_GOF && !is.null(fit_inhom_structural)) {
  cat("  GOF for structural-only model...\n")
  
  # Reconstruct params_init for structural model
  skel_structural_gof <- params_init_structural
  params_init_structural_gof <- relist(fit_inhom_structural$fit$par, skeleton = skel_structural_gof)
  params_init_structural_gof$K <- params_init_structural$K
  params_init_structural_gof$mu <- params_init_structural$mu
  
  # For GOF simulations, use cond_intensity_inhom to match the fitted inhomogeneous model
  # The gof() function will automatically use cond_intensity_inhom when inhom_bg is provided
  GOF_results_structural <- gof(
    fit = fit_inhom_structural,
    net_obs = net_raw,
    params_init = params_init_structural_gof,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,  # Will be overridden to cond_intensity_inhom by gof() when inhom_bg is provided
    formula_RHS = FORMULA_RHS_STRUCTURAL,
    time_window = GOF_TIME_WINDOW,
    truncation = TRUNCATION,
    mark_decay = "activity",
    max_node_time = 1,
    inhom_bg = inhom_bg,  # This enables inhomogeneous simulations matching the fitted model
    n_sim = N_GOF,
    cores = N_CORES,
    max_deg = 15,
    k_esp = 15,
    degree = 0,
    esp = 0,
    mu_multiplier = 5,
    verbose = TRUE
  )
} else {
  if (!RUN_GOF) cat("  RUN_GOF = FALSE; skipping structural GOF\n")
  if (is.null(fit_inhom_structural)) cat("  No structural fit available; skipping structural GOF\n")
}