# =============================================================================
# OpenAlex Hawkes study: inhomogeneous fit, temporal fit, KS tests, GOF
# =============================================================================
# Run from package root: Rscript inst/openalex_study/openalex_hawkes_study.R
# Or submit via SLURM: sbatch inst/openalex_study/run_openalex.slurm
#
# Requires: hawkesNet package (includes inhomogeneous fit and KDE background).
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
# Under SLURM, use submit dir so path stays valid if getwd() breaks later.
PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) Sys.getenv("SLURM_SUBMIT_DIR") else getwd()
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) {
  # Try parent directory
  parent_dir <- dirname(PKG_ROOT)
  if (file.exists(file.path(parent_dir, "DESCRIPTION"))) {
    PKG_ROOT <- parent_dir
  } else {
    # Try looking for inst/ directory
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
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 50))
# If request is 128 but squeue shows 256 (one node), use 256. Override: CORES_OVERRIDE=128 or 256.
if (nzchar(Sys.getenv("CORES_OVERRIDE"))) {
  N_CORES <- as.numeric(Sys.getenv("CORES_OVERRIDE"))
} else if (N_CORES == 128L) {
  N_CORES <- 256L
  cat("Request was 128; using 256 (typical when squeue shows 256)\n")
}
MAX_ITER <- 5000
TRUNCATION <- 200L
GOF_TIME_WINDOW <- c(0, 1)  # Full time period for GOF simulations
N_GOF <- 25L   # number of simulated networks for goodness-of-fit
PAPER_OUTPUT <- TRUE
RUN_GOF <- TRUE
TOPIC <- "Point processes and geometric inequalities"

# Reproducibility
set.seed(42L)

# Helper: run expression and return on_error (default NULL) if it fails
safe_run <- function(expr, label = "", on_error = NULL) {
  tryCatch(expr, error = function(e) {
    if (nzchar(label)) cat("  ERROR:", label, ":", e$message, "\n")
    on_error
  })
}

# Helper: default parameter initialization (structural or nodeMatch)
make_default_params <- function(n_cs, mu_init, include_gender = FALSE) {
  p <- list(
    mu = mu_init,
    beta_overall = 1,
    K = 0.5,
    beta_edges = 1,
    node_lambda = 1,
    CS_params = c(-10, rep(0, n_cs - 1))
  )
  if (include_gender) {
    p$vertex_categorical <- list(gender = c(female = 0.1, male = 0.5))
    p$vertex_categorical_levels <- list(gender = c("female", "male", "unknown"))
  }
  p
}

# Speed tips for Nelder-Mead: pass cores = N_CORES to fit_hawkesNet_inhom so each
# loglik evaluation parallelizes over cached intensity closures; lower TRUNCATION
# or fewer formula terms reduce work per evaluation; good parscale reduces iterations.

# =============================================================================
# GOF functions are now in R/gof.R (exported from package)
# =============================================================================

# =============================================================================
# 1. Fetch data and prepare network
# =============================================================================
t_total <- proc.time()
cat("=== OpenAlex Hawkes Study ===\n")
cat("  Search:", SEARCH_STRING, "| Pages:", PAGES, "| Cores:", N_CORES, "\n")
cat("  Date range:", MIN_DATE, "to", MAX_DATE, "\n")
cat("  Package root:", PKG_ROOT, "\n")
cat("  Working directory:", getwd(), "\n")
cat("  Cluster output dir:", file.path(PKG_ROOT, "cluster_output"), "\n")

# Test write permissions early
test_dir <- file.path(PKG_ROOT, "cluster_output")
if (!dir.exists(test_dir)) {
  test_create <- tryCatch({
    dir.create(test_dir, showWarnings = TRUE, recursive = TRUE)
    cat("  Created cluster_output directory\n")
  }, error = function(e) {
    cat("  WARNING: Cannot create cluster_output directory:", e$message, "\n")
  })
} else {
  cat("  cluster_output directory exists\n")
}

# Test write permissions
test_file <- file.path(test_dir, ".test_write")
test_write <- tryCatch({
  writeLines("test", test_file)
  unlink(test_file)
  cat("  Write permissions: OK\n")
}, error = function(e) {
  cat("  WARNING: Cannot write to cluster_output directory:", e$message, "\n")
})

cat("\n")

cat("--- Step 1: Fetch data and prepare network ---\n")
t_step <- proc.time()
out <- get_network(email = EMAIL, pages = PAGES, per_page = PER_PAGE,
                   string = SEARCH_STRING, min_date = MIN_DATE, max_date = MAX_DATE,topics_include = c(TOPIC))
net_raw <- out$net
edges <- out$edges
network::set.vertex.attribute(net_raw, "time", net_raw %v% "time_scaled")
network::set.edge.attribute(net_raw, "time", net_raw %e% "time_scaled")
net_raw <- hawkesNet::normalize_times_01(net_raw, attr = "time", keep_na = TRUE)
n_events <- length(hawkesNet::get_times(net_raw)$times)
n_nodes <- network::network.size(net_raw)
cat("  Network:", n_events, "events,", n_nodes, "nodes\n")
# Check gender distribution
if ("gender" %in% network::list.vertex.attributes(net_raw)) {
  gender_vals <- net_raw %v% "gender"
  gender_counts <- table(gender_vals, useNA = "ifany")
  cat("  Gender distribution:", paste(names(gender_counts), "=", gender_counts, collapse = ", "), "\n")
  if (all(gender_vals == "unknown", na.rm = TRUE)) {
    warning("⚠ WARNING: All genders are 'unknown'. This may affect nodeMatch and homophily statistics.")
    warning("⚠ Install 'gender' and 'genderdata' packages for gender prediction: install.packages(c('gender', 'genderdata'))")
  }
} else {
  warning("⚠ No 'gender' attribute found on network vertices")
}
cat("  Step 1 took:", round((proc.time() - t_step)[3], 1), "s\n\n")

# --- Cleanup: remove Step 1 temporaries before fitting ---
rm(out, test_dir, test_file)
if (exists("test_write"))  rm(test_write)
if (exists("test_create")) rm(test_create)
if (exists("gender_vals"))  rm(gender_vals)
if (exists("gender_counts")) rm(gender_counts)
gc()

# =============================================================================
# 2. Inhomogeneous (KDE) + CS fits: structural -> nodeMatch
# =============================================================================
cat("--- Step 2: Inhomogeneous (KDE) + CS fits ---\n")
t_step <- proc.time()
time_window_01 <- c(0, 1)
cat("  Preparing inhomogeneous background (KDE)...\n")
t_kde <- proc.time()
inhom_bg <- safe_run(
  prepare_inhomogeneous_background(net_raw, time_attr = "time", bw = NULL, grid_n = 2048),
  "prepare_inhomogeneous_background"
)
cat("  KDE background:", round((proc.time() - t_kde)[3], 1), "s\n")

# =============================================================================
# 2a. Structural-only fit (first, independent)
# =============================================================================
fit_inhom_structural <- NULL
# gwdegree(0.5) captures the full degree distribution with one parameter
# (geometrically weighted), replacing star(c(2,3,4,5)) which needed 4 params.
FORMULA_RHS_STRUCTURAL <- "edges + triangles + gwdegree(0.5)"
if (!is.null(inhom_bg)) {
  cat("\n--- Step 2a: Structural-only fit (no gender) ---\n")
  cat("  Formula:", FORMULA_RHS_STRUCTURAL, "\n")
  t_step_structural <- proc.time()
  
  # Get expected parameters for structural formula
  exp_cs_structural <- expected_params_PMF_mark_CS(net_raw, FORMULA_RHS_STRUCTURAL)
  n_cs_structural <- if (!is.na(exp_cs_structural$CS_params_length)) exp_cs_structural$CS_params_length else 3L
  
  # Initialize parameters for structural-only model (no vertex_categorical)
  mu_init <- inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1])
  params_init_structural <- make_default_params(n_cs_structural, mu_init)
  
  # Create parscale for structural model
  p_scale_structural <- c(
    beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
    setNames(rep(0.1, n_cs_structural), paste0("CS_params", seq_len(n_cs_structural)))
  )
  
  cat("  Method: Nelder-Mead (max", MAX_ITER, "iterations)\n")
  t_fit_structural <- proc.time()
  
  fit_inhom_structural <- safe_run(
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
      verbose = FALSE,
      fixed_params = c("K", "mu"),
      parscale = p_scale_structural,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      cores = N_CORES
    ),
    "Structural fit"
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
  cat("  Step 2a total:", round((proc.time() - t_step_structural)[3], 1), "s\n")
} else {
  cat("  No inhomogeneous background; skipping structural fit\n")
}

# =============================================================================
# 2b. nodeMatch fit (initialized from structural fit)
# =============================================================================
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
    max(4L, (if (exists("n_cs_structural", inherits = FALSE)) n_cs_structural else 3L) + 1L)
  }
  
  # Initialize nodeMatch fit from structural fit results
  if (!is.null(fit_inhom_structural) && fit_inhom_structural$fit$convergence == 0) {
    # CRITICAL: skeleton must exclude fixed_params (mu, K) — fit_hawkesNet_inhom strips them
    # before optim, so fit$par only has the non-fixed parameters.
    skel_structural <- list(
      beta_overall = params_init_structural$beta_overall,
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
      # Validate extracted parameters
      beta_overall_valid <- is.finite(pfit_structural$beta_overall)
      beta_edges_valid <- is.finite(pfit_structural$beta_edges)
      node_lambda_valid <- is.finite(pfit_structural$node_lambda) && pfit_structural$node_lambda > 0
      cs_valid <- all(is.finite(pfit_structural$CS_params)) && length(pfit_structural$CS_params) >= n_cs_structural
      
      if (beta_overall_valid && beta_edges_valid && node_lambda_valid && cs_valid) {
        # Initialize nodeMatch: use structural CS_params, add nodeMatch term (0 init)
        cs_structural <- pfit_structural$CS_params[seq_len(min(n_cs_structural, length(pfit_structural$CS_params)))]
        cs_padding <- rep(0, max(0L, n_cs_nodematch - length(cs_structural)))
        cs_init_nodematch <- c(cs_structural, cs_padding)[seq_len(n_cs_nodematch)]
        
        params_init_nodematch <- list(
          mu = params_init_structural$mu,         # fixed param: use init value
          beta_overall = pfit_structural$beta_overall,
          K = params_init_structural$K,            # fixed param: use init value
          beta_edges = pfit_structural$beta_edges,
          node_lambda = pfit_structural$node_lambda,
          CS_params = cs_init_nodematch,
          vertex_categorical = list(gender = c(female = 0.1, male = 0.5)),
          vertex_categorical_levels = list(gender = c("female", "male", "unknown"))
        )
        cat("  Initialized from structural fit (beta_overall=", round(pfit_structural$beta_overall, 4),
            ", beta_edges=", round(pfit_structural$beta_edges, 4),
            ", node_lambda=", round(pfit_structural$node_lambda, 4), ")\n", sep = "")
      } else {
        pfit_structural <- NULL  # Force fallback
      }
    }
    
    if (is.null(pfit_structural)) {
      params_init_nodematch <- make_default_params(n_cs_nodematch, mu_init, include_gender = TRUE)
      cat("  Structural fit parameters invalid; using independent initialization\n")
    }
  } else {
    params_init_nodematch <- make_default_params(n_cs_nodematch, mu_init, include_gender = TRUE)
    cat("  Structural fit not available; using independent initialization\n")
  }
  
  # Create parscale for nodeMatch
  p_scale_nodematch <- c(
    beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
    setNames(rep(0.1, n_cs_nodematch), paste0("CS_params", seq_len(n_cs_nodematch))),
    vertex_categorical.gender.female = 0.1, vertex_categorical.gender.male = 0.1
  )
  
  # Save structural fit to disk and drop from memory (rehydrate for GOF and save_list)
  structural_fit_cache <- file.path(PKG_ROOT, "cluster_output", "structural_fit_cache.RDS")
  if (!is.null(fit_inhom_structural)) {
    dir.create(file.path(PKG_ROOT, "cluster_output"), showWarnings = FALSE, recursive = TRUE)
    fit_to_cache <- fit_inhom_structural
    fit_to_cache$intens_funcs <- NULL
    tryCatch(
      saveRDS(fit_to_cache, structural_fit_cache),
      error = function(e) cat("  Warning: could not cache structural fit:", conditionMessage(e), "\n")
    )
    rm(fit_to_cache)
    fit_inhom_structural <- NULL
  }
  # --- Cleanup: remove structural-fit temporaries before nodeMatch fit ---
  rm(exp_cs_structural, t_step_structural, t_fit_structural, elapsed_fit_structural)
  gc()
  cat("  Structural fit cached to disk; environment cleaned before nodeMatch fit\n")
  
  # --- Cleanup: remove nodeMatch initialization temporaries ---
  if (exists("skel_structural", inherits = FALSE)) rm(skel_structural)
  if (exists("pfit_structural", inherits = FALSE)) rm(pfit_structural)
  if (exists("cs_structural", inherits = FALSE))   rm(cs_structural)
  if (exists("cs_padding", inherits = FALSE))       rm(cs_padding)
  if (exists("cs_init_nodematch", inherits = FALSE)) rm(cs_init_nodematch)
  if (exists("beta_overall_valid", inherits = FALSE)) rm(beta_overall_valid, beta_edges_valid, node_lambda_valid, cs_valid)

  cat("  Method: Nelder-Mead (max", MAX_ITER, "iterations)\n")
  t_fit_nodematch <- proc.time()
  
  fit_inhom_nodematch <- safe_run(
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
      combine_intensity = TRUE,
      cores = N_CORES
    ),
    "nodeMatch fit"
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

# --- Cleanup: remove Step 2 temporaries and strip heavy closures ---
# Strip intens_funcs from nodeMatch fit (not needed for GOF — only fit$par is used)
if (!is.null(fit_inhom_nodematch) && !is.null(fit_inhom_nodematch$intens_funcs)) {
  fit_inhom_nodematch$intens_funcs <- NULL
}
rm(t_step, t_kde)
if (exists("t_step_nodematch", inherits = FALSE)) rm(t_step_nodematch)
if (exists("t_fit_nodematch", inherits = FALSE)) rm(t_fit_nodematch)
if (exists("elapsed_fit_nodematch", inherits = FALSE)) rm(elapsed_fit_nodematch)
if (exists("exp_cs_nodematch", inherits = FALSE)) rm(exp_cs_nodematch)
gc()

# =============================================================================
# 3. Temporal Hawkes fit and KS test
# =============================================================================
cat("--- Step 3: Temporal Hawkes fit + KS test ---\n")
t_step <- proc.time()
t_events <- sort(unique(c(hawkesNet::get_times(net_raw)$node_times,
                         hawkesNet::get_times(net_raw)$edge_times)))
t_events <- t_events[!is.na(t_events)]
windowT <- c(min(t_events), max(t_events))
realiz <- data.frame(t = t_events)
fit_temporal <- NULL
init_gamma <- max(length(t_events) * 0.5, 10)
params_init_exp <- list(gamma = init_gamma, beta = 10, K = 0.2)
cat("  Fitting temporal Hawkes (exp kernel)...\n")
t_fit <- proc.time()
fit_temporal <- safe_run(
  hawkesNet::fit_temporal_hawkes(
    params_init = params_init_exp,
    realiz = realiz,
    windowT = windowT,
    method = "Nelder-Mead",
    maxit = 500,
    kernel = "exp",
    trace = 0
  ),
  "fit_temporal_hawkes"
)
cat("  Temporal fit:", round((proc.time() - t_fit)[3], 1), "s\n")
if (!is.null(fit_temporal)) {
  cat("  Temporal par:", paste(names(fit_temporal$par), "=", round(fit_temporal$par, 4), collapse = ", "), "\n")
}
ks_temporal_pval <- NA_real_
if (!is.null(fit_temporal) && exists("ks_test_pval_temporal")) {
  cat("  Computing KS test...\n")
  ks_temporal_pval <- safe_run(
    hawkesNet::ks_test_pval_temporal(
      realiz = realiz,
      windowT = windowT,
      hawkes_par = fit_temporal$par,
      kernel = "exp",
      use_kde = TRUE
    ),
    on_error = NA_real_
  )
  cat("  Temporal KS p-value:", ks_temporal_pval, "\n")
}
cat("  Step 3 total:", round((proc.time() - t_step)[3], 1), "s\n\n")

# --- Cleanup: remove Step 3 temporaries before GOF ---
rm(t_step, t_fit, init_gamma, params_init_exp)
if (exists("t_events", inherits = FALSE)) rm(t_events)
gc()

# =============================================================================
# 4. Goodness-of-fit: simulate from fitted model, compare degree/ESP/geodesic/waiting times
# =============================================================================
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

# Rehydrate structural fit from cache if it was dropped for memory
structural_fit_cache <- file.path(PKG_ROOT, "cluster_output", "structural_fit_cache.RDS")
if (is.null(fit_inhom_structural) && file.exists(structural_fit_cache)) {
  fit_inhom_structural <- tryCatch(readRDS(structural_fit_cache), error = function(e) NULL)
  if (!is.null(fit_inhom_structural)) cat("  Rehydrated structural fit from cache for GOF/save\n")
}

# GOF for structural-only model (first)
if (RUN_GOF && !is.null(fit_inhom_structural)) {
  cat("  GOF for structural-only model...\n")
  
  # Reconstruct params_init for structural model (skeleton must exclude fixed params mu, K)
  skel_structural_gof <- params_init_structural
  skel_structural_gof$mu <- NULL
  skel_structural_gof$K <- NULL
  params_init_structural_gof <- tryCatch(
    relist(fit_inhom_structural$fit$par, skeleton = skel_structural_gof),
    error = function(e) { cat("  Warning: GOF relist failed:", e$message, "\n"); skel_structural_gof }
  )
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

# GOF for nodeMatch model (second)
if (RUN_GOF && !is.null(fit_inhom_nodematch)) {
  cat("\n  GOF for nodeMatch model...\n")
  
  # Reconstruct params_init for nodeMatch (skeleton must exclude fixed params mu, K and vertex_categorical_levels)
  skel_nodematch_gof <- params_init_nodematch
  skel_nodematch_gof$vertex_categorical_levels <- NULL
  skel_nodematch_gof$mu <- NULL
  skel_nodematch_gof$K <- NULL
  params_init_nodematch_gof <- tryCatch(
    relist(fit_inhom_nodematch$fit$par, skeleton = skel_nodematch_gof),
    error = function(e) { cat("  Warning: GOF nodeMatch relist failed:", e$message, "\n"); skel_nodematch_gof }
  )
  params_init_nodematch_gof$vertex_categorical_levels <- params_init_nodematch$vertex_categorical_levels
  params_init_nodematch_gof$K <- params_init_nodematch$K
  params_init_nodematch_gof$mu <- params_init_nodematch$mu
  # Restore names and repair parameters before GOF
  params_init_nodematch_gof <- hawkesNet:::reconstruct_vertex_categorical_names(
    params_init_nodematch_gof, params_init_nodematch$vertex_categorical_levels)
  params_init_nodematch_gof <- hawkesNet:::repair_vertex_categorical_params(params_init_nodematch_gof, eps = 1e-6)
  
  # For GOF simulations, use cond_intensity_inhom to match the fitted inhomogeneous model
  # The gof() function will automatically use cond_intensity_inhom when inhom_bg is provided
  GOF_results_nodematch <- gof(
    fit = fit_inhom_nodematch,
    net_obs = net_raw,
    params_init = params_init_nodematch_gof,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,  # Will be overridden to cond_intensity_inhom by gof() when inhom_bg is provided
    formula_RHS = FORMULA_RHS_NODEMATCH,
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
  if (!RUN_GOF) cat("  RUN_GOF = FALSE; skipping nodeMatch GOF\n")
  if (is.null(fit_inhom_nodematch)) cat("  No nodeMatch fit available; skipping nodeMatch GOF\n")
}

# Use nodeMatch GOF as primary for display
if (RUN_GOF && !is.null(fit_inhom_nodematch)) {
  GOF_results <- GOF_results_nodematch
}

cat("  Step 4 total:", round((proc.time() - t_step)[3], 1), "s\n\n")

# =============================================================================
# 5. Save full state for rehydration
# =============================================================================
cat("--- Step 5: Save full state ---\n")
# Rehydrate structural fit from cache if needed (e.g. still null after Step 4)
if (is.null(fit_inhom_structural) && file.exists(structural_fit_cache)) {
  fit_inhom_structural <- tryCatch(readRDS(structural_fit_cache), error = function(e) NULL)
}
# Strip intensity caches from fits so RDS stays small (low cost to re-run intensity if needed)
fit_structural_for_save <- fit_inhom_structural
if (!is.null(fit_structural_for_save)) fit_structural_for_save$intens_funcs <- NULL
fit_nodematch_for_save <- fit_inhom_nodematch
if (!is.null(fit_nodematch_for_save)) fit_nodematch_for_save$intens_funcs <- NULL
save_list <- list(
  net_raw = net_raw,
  edges = edges,
  inhom_bg = inhom_bg,
  fit_inhom_nodematch = fit_nodematch_for_save,
  fit_inhom_structural = fit_structural_for_save,
  params_init_nodematch = params_init_nodematch,
  params_init_structural = params_init_structural,
  FORMULA_RHS_NODEMATCH = FORMULA_RHS_NODEMATCH,
  FORMULA_RHS_STRUCTURAL = FORMULA_RHS_STRUCTURAL,
  fit_temporal = fit_temporal,
  ks_temporal_pval = ks_temporal_pval,
  realiz = realiz,
  windowT = windowT,
  GOF_results = GOF_results,
  GOF = GOF_results,  # alias so dat$GOF$plots works
  GOF_results_nodematch = GOF_results_nodematch,
  GOF_results_structural = GOF_results_structural,
  N_GOF = N_GOF,
  SEARCH_STRING = SEARCH_STRING,
  time_window_01 = time_window_01
)
# Ensure output dir exists (getwd() can become invalid on some clusters)
cluster_output_dir <- file.path(PKG_ROOT, "cluster_output")
rds_path_primary <- file.path(cluster_output_dir, "results_openalex_full.RDS")
OPENALEX_RDS_PATH <- NULL  # set after save so Step 6 can find file

cat("  Saving to:", rds_path_primary, "\n")
cat("  Package root:", PKG_ROOT, "\n")
cat("  Cluster output dir:", cluster_output_dir, "\n")

save_ok <- tryCatch({
  dir.create(cluster_output_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(save_list, rds_path_primary)
  TRUE
}, error = function(e) {
  cat("  Save to ", rds_path_primary, " failed: ", conditionMessage(e), "\n")
  FALSE
})

if (save_ok) {
  OPENALEX_RDS_PATH <- rds_path_primary
  cat("  Successfully saved to:", rds_path_primary, "\n")
  if (file.exists(rds_path_primary)) {
    file_info <- file.info(rds_path_primary)
    cat("  File size:", round(file_info$size / 1024^2, 2), "MB\n")
  }
} else {
  fallback <- "results_openalex_full.RDS"
  tryCatch({
    saveRDS(save_list, fallback)
    OPENALEX_RDS_PATH <- fallback
    cat("  Saved to fallback:", normalizePath(fallback, mustWork = FALSE), "\n")
  }, error = function(e) {
    cat("  Fallback save also failed:", conditionMessage(e), "\n")
  })
}
if (is.null(OPENALEX_RDS_PATH)) cat("  WARNING: Failed to save results file!\n")

cat("\n")

# =============================================================================
# 6. PAPER_OUTPUT: rehydrate and produce figures/tables
# =============================================================================
if (PAPER_OUTPUT) {
  cat("--- Step 6: Paper output (figures & tables) ---\n")
  t_step <- proc.time()
  rds_file <- if (!is.null(OPENALEX_RDS_PATH)) OPENALEX_RDS_PATH else file.path(PKG_ROOT, "cluster_output", "results_openalex_full.RDS")
  if (!file.exists(rds_file)) rds_file <- "results_openalex_full.RDS"
  dat <- readRDS(rds_file)
  # Ensure GOF alias exists for older RDS that may not have it, and plots is never NULL
  if (is.null(dat$GOF) && !is.null(dat$GOF_results)) dat$GOF <- dat$GOF_results
  if (!is.null(dat$GOF_results) && is.null(dat$GOF_results$plots)) dat$GOF_results$plots <- list()
  if (!is.null(dat$GOF) && is.null(dat$GOF$plots)) dat$GOF$plots <- list()
  list2env(dat, envir = .GlobalEnv)
  cat("  Rehydrated; producing figures and tables.\n")

  # Print nodeMatch fit results
  if (!is.null(dat$fit_inhom_nodematch)) {
    cat("\n--- nodeMatch Fit Results ---\n")
    if (!is.null(dat$fit_inhom_nodematch$fit_table)) {
      cat("  Inhomogeneous fit (nodeMatch): parameter estimates and standard errors\n")
      print(dat$fit_inhom_nodematch$fit_table, max = NULL)
    } else {
      cat("  Inhomogeneous fit (nodeMatch): raw parameters\n")
      print(fit_inhom_nodematch$fit$par)
    }
    if (exists("params_init_nodematch") && !is.null(params_init_nodematch$vertex_categorical) &&
        !is.null(dat$fit_inhom_nodematch) && !is.null(dat$fit_inhom_nodematch$fit$par)) {
      # Reconstruct nodeMatch parameters for gender proportions using its own initialization
      # CRITICAL: Create skeleton that matches EXACTLY what was used during optimization
      skel_nodematch <- list(
        mu = params_init_nodematch$mu,
        beta_overall = params_init_nodematch$beta_overall,
        K = params_init_nodematch$K,
        beta_edges = params_init_nodematch$beta_edges,
        node_lambda = params_init_nodematch$node_lambda,
        CS_params = params_init_nodematch$CS_params,
        vertex_categorical = params_init_nodematch$vertex_categorical
      )
      
      pfit_nodematch <- tryCatch({
        relist(dat$fit_inhom_nodematch$fit$par, skeleton = skel_nodematch)
      }, error = function(e) NULL)
      
      if (!is.null(pfit_nodematch)) {
        pfit_nodematch$vertex_categorical_levels <- params_init_nodematch$vertex_categorical_levels
        # Restore names and repair parameters before expanding
        pfit_nodematch <- hawkesNet:::reconstruct_vertex_categorical_names(
          pfit_nodematch, params_init_nodematch$vertex_categorical_levels)
        pfit_nodematch <- hawkesNet:::repair_vertex_categorical_params(pfit_nodematch, eps = 1e-6)
        if (!is.null(pfit_nodematch$vertex_categorical$gender)) {
          levs <- params_init_nodematch$vertex_categorical_levels$gender
          pgender_nodematch <- expand_vertex_categorical_probs(pfit_nodematch$vertex_categorical$gender, levs)
          if (!is.null(pgender_nodematch)) {
            cat("  Fitted gender proportions (n-1 expanded):\n"); print(pgender_nodematch)
          } else {
            cat("  Fitted gender proportions: could not expand (invalid parameters)\n")
          }
        }
      }
    }
  }
  
  # Print structural-only fit results
  if (!is.null(dat$fit_inhom_structural)) {
    cat("\n--- Structural-Only Fit Results ---\n")
    if (!is.null(dat$fit_inhom_structural$fit_table)) {
      cat("  Inhomogeneous fit (structural-only): parameter estimates and standard errors\n")
      print(dat$fit_inhom_structural$fit_table, max = NULL)
    } else {
      cat("  Inhomogeneous fit (structural-only): raw parameters\n")
      print(fit_inhom_structural$fit$par)
    }
  }
  
  if (!is.null(dat$fit_temporal)) {
    cat("  Temporal Hawkes par:\n"); print(fit_temporal$par)
    cat("  Temporal KS p-value:", ks_temporal_pval, "\n")
  }

  # GOF plots: use plots from gof() function if available, otherwise generate here
  if (!is.null(dat$GOF_results)) {
    # First, try to use plots from gof() function (if available)
    if (!is.null(dat$GOF_results$plots) && length(dat$GOF_results$plots) > 0) {
      cat("  Displaying GOF plots from gof() function...\n")
      if (!is.null(dat$GOF_results$plots$degree_plot)) print(dat$GOF_results$plots$degree_plot)
      if (!is.null(dat$GOF_results$plots$esp_plot)) print(dat$GOF_results$plots$esp_plot)
      if (!is.null(dat$GOF_results$plots$geodist_plot)) print(dat$GOF_results$plots$geodist_plot)
      if (!is.null(dat$GOF_results$plots$waiting_times_plot)) print(dat$GOF_results$plots$waiting_times_plot)
    } else if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(dat$GOF_results$degree_obs)) {
      # Fallback: generate plots here (legacy code)
      cat("  Generating GOF plots (legacy method)...\n")
      gof <- dat$GOF_results
    gof <- dat$GOF_results
    max_deg <- length(gof$degree_obs) - 1
    # Degree: observed vs simulated boxplots
    deg_df <- rbind(
      data.frame(degree = 0:max_deg, count = gof$degree_obs, type = "Observed"),
      data.frame(degree = rep(0:max_deg, each = nrow(gof$degree_sim)),
                 count = as.vector(gof$degree_sim),
                 type = "Simulated")
    )
    p_deg <- ggplot(deg_df, aes(x = factor(degree), y = count, fill = type)) +
      geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
      labs(title = "GOF: Degree distribution", x = "Degree", y = "Count") +
      theme_minimal() + theme(legend.position = "bottom")
    print(p_deg)
    # ESP
    esp_obs <- gof$esp_obs
    esp_sim <- gof$esp_sim
    if (!is.null(esp_sim) && nrow(esp_sim) > 0) {
      esp_df <- rbind(
        data.frame(esp = 0:(length(esp_obs)-1), value = esp_obs, type = "Observed"),
        data.frame(esp = rep(0:(ncol(esp_sim)-1), each = nrow(esp_sim)),
                   value = as.vector(esp_sim), type = "Simulated")
      )
      p_esp <- ggplot(esp_df, aes(x = factor(esp), y = value, fill = type)) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
        labs(title = "GOF: ESP distribution", x = "ESP", y = "Count") +
        theme_minimal() + theme(legend.position = "bottom")
      print(p_esp)
    }
    # Geodesic: boxplot of pair counts at each distance + ECDF
    g_obs <- gof$geodist_obs
    g_sim <- gof$geodist_sim
    if (length(g_obs) > 0 && length(g_sim) > 0) {
      # Boxplot: tabulate counts at each integer distance
      max_geod <- min(max(c(g_obs, unlist(g_sim)), na.rm = TRUE), 20)
      geod_levels <- seq_len(max_geod)
      obs_tab <- table(factor(g_obs, levels = geod_levels))
      sim_tabs <- lapply(g_sim, function(g) {
        as.vector(table(factor(g, levels = geod_levels)))
      })
      sim_mat <- do.call(rbind, sim_tabs)
      geod_box_df <- rbind(
        data.frame(distance = geod_levels, count = as.vector(obs_tab), type = "Observed"),
        data.frame(distance = rep(geod_levels, each = nrow(sim_mat)),
                   count = as.vector(sim_mat), type = "Simulated")
      )
      p_geod_box <- ggplot(geod_box_df, aes(x = factor(distance), y = count, fill = type)) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
        labs(title = "GOF: Geodesic distance distribution",
             x = "Geodesic distance", y = "Number of pairs") +
        theme_minimal() + theme(legend.position = "bottom")
      print(p_geod_box)

      # ECDF version
      max_d <- max(c(g_obs, unlist(g_sim)), na.rm = TRUE)
      x_seq <- seq(0, min(max_d, 20), length.out = 200)
      ecdf_obs <- sapply(x_seq, function(x) mean(g_obs <= x, na.rm = TRUE))
      ecdf_sim <- sapply(x_seq, function(x) mean(unlist(g_sim) <= x, na.rm = TRUE))
      geod_df <- rbind(
        data.frame(dist = x_seq, ecdf = ecdf_obs, type = "Observed"),
        data.frame(dist = x_seq, ecdf = ecdf_sim, type = "Simulated")
      )
      p_geod <- ggplot(geod_df, aes(x = dist, y = ecdf, color = type)) +
        geom_line(linewidth = 1) +
        labs(title = "GOF: Geodesic distance (ECDF)", x = "Distance", y = "ECDF") +
        theme_minimal() + theme(legend.position = "bottom")
      print(p_geod)
    }
    # Waiting times between formations: triangle, 2-star, 3-star
    w_obs <- gof$wait_obs
    w_sim <- gof$wait_sim
    if (!is.null(w_obs) && is.list(w_obs) && !is.null(w_sim) && length(w_sim) > 0) {
      obs_vec <- c(w_obs$triangle, w_obs$star2, w_obs$star3)
      obs_metric <- rep(c("Triangle", "2-star", "3-star"),
                       c(length(w_obs$triangle), length(w_obs$star2), length(w_obs$star3)))
      sim_vec <- unlist(lapply(w_sim, function(x) c(x$triangle, x$star2, x$star3)))
      sim_metric <- unlist(lapply(w_sim, function(x) rep(c("Triangle", "2-star", "3-star"),
                         c(length(x$triangle), length(x$star2), length(x$star3)))))
      wait_df <- rbind(
        data.frame(metric = obs_metric, value = obs_vec, type = "Observed"),
        data.frame(metric = sim_metric, value = sim_vec, type = "Simulated")
      )
      wait_df <- wait_df[!is.na(wait_df$value), ]
      if (nrow(wait_df) > 0) {
        p_wait <- ggplot(wait_df, aes(x = metric, y = value, fill = type)) +
          geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
          labs(title = "GOF: Waiting time between formations (triangle / 2-star / 3-star)",
               x = "", y = "Waiting time") +
          theme_minimal() + theme(legend.position = "bottom")
        print(p_wait)
      }
    }
    } else {
      cat("  ggplot2 not available; skipping GOF plots\n")
    }
  } else {
    cat("  No GOF results available for plotting\n")
  }
  cat("  Step 6 total:", round((proc.time() - t_step)[3], 1), "s\n\n")
}

# =============================================================================
# Total elapsed time
# =============================================================================
cat("=== OpenAlex study complete ===\n")
cat("  Total wall time:", round((proc.time() - t_total)[3], 1), "s (",
    round((proc.time() - t_total)[3] / 60, 1), "min)\n")
