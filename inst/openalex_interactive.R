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
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 50))
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

# =============================================================================
# 2. Structural-only fit
# =============================================================================
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

FORMULA_RHS_STRUCTURAL <- "edges + triangles + star(c(2,3))"
if (!is.null(inhom_bg)) {
  cat("\n  Formula:", FORMULA_RHS_STRUCTURAL, "\n")
  
  # Get expected parameters for structural formula
  exp_cs_structural <- expected_params_PMF_mark_CS(net_raw, FORMULA_RHS_STRUCTURAL)
  n_cs_structural <- if (!is.na(exp_cs_structural$CS_params_length)) exp_cs_structural$CS_params_length else 4L
  
  # Initialize parameters for structural-only model (no vertex_categorical)
  params_init_structural <- list(
    mu = inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1]),
    beta_overall = 1,
    K = 0.5,
    beta_edges = 1,
    node_lambda = 1,
    CS_params = c(-10, rep(0, n_cs_structural - 1))
  )
  
  # Create parscale for structural model
  p_scale_structural <- c(
    beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
    setNames(rep(0.1, n_cs_structural), paste0("CS_params", seq_len(n_cs_structural)))
  )
  
  cat("  Method: Nelder-Mead (max", MAX_ITER, "iterations)\n")
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
      verbose = FALSE,
      fixed_params = c("K", "mu"),
      parscale = p_scale_structural,
      cache_intensity = TRUE,
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
