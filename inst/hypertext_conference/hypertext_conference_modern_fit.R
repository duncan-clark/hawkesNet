## =============================================================================
## Hypertext 2009 conference: modern hawkesNet fit + GOF
## =============================================================================
## Run from package root:
##   Rscript inst/hypertext_conference/hypertext_conference_modern_fit.R
##
## Quick local mode (keeps runtime under ~2 minutes):
##   LOCAL_QUICK=TRUE Rscript inst/hypertext_conference/hypertext_conference_modern_fit.R
##
## Cluster mode knobs:
##   SLURM_CPUS_PER_TASK=32 MAX_ITER=5000 N_GOF=50 OPENALEX_GOF_CORES=32 Rscript ...
##
## Notes:
## - This models *edge formation* (simple graph): we collapse repeated contacts to first contact per dyad.
## - We do NOT use growth_only; transitivity is included via gwesp.
## - Inhomogeneous background is estimated via KDE over the full time window (multi-day supported).
## =============================================================================

library(hawkesNet)
library(dplyr)
library(network)
library(ernm)
library(sna)

PKG_ROOT <- getwd()

# Use the repo (most up-to-date) implementation even if the installed package is older.
# This keeps GOF + truncation behavior consistent with your current hawkesNet code.
USE_REPO_CODE <- isTRUE(as.logical(Sys.getenv("USE_REPO_CODE", "TRUE")))
if (USE_REPO_CODE) {
  source(file.path(PKG_ROOT, "R", "utils.R"))
  source(file.path(PKG_ROOT, "R", "kde_background.R"))
  source(file.path(PKG_ROOT, "R", "temporal_hawkes.R"))
  source(file.path(PKG_ROOT, "R", "gof.R"))
  source(file.path(PKG_ROOT, "R", "mark_PMF.R"))
  source(file.path(PKG_ROOT, "R", "hawkesNet.R"))
}

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
default_local_quick <- if (nzchar(Sys.getenv("SLURM_JOB_ID")) || nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) "FALSE" else "TRUE"
LOCAL_QUICK <- isTRUE(as.logical(Sys.getenv("LOCAL_QUICK", default_local_quick)))

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", if (LOCAL_QUICK) 1L else 7L))
N_CORES <- max(1L, N_CORES)

MAX_ITER <- as.integer(Sys.getenv("MAX_ITER", if (LOCAL_QUICK) 200L else 5000L))
TRUNCATION <- as.integer(Sys.getenv("TRUNCATION", NA_integer_))  # if NA, choose below

# GOF
RUN_GOF <- isTRUE(as.logical(Sys.getenv("RUN_GOF", if (LOCAL_QUICK) "FALSE" else "TRUE")))
N_GOF <- as.integer(Sys.getenv("N_GOF", if (LOCAL_QUICK) 2L else 25L))
N_GOF <- max(1L, N_GOF)
N_CORES_GOF <- as.integer(Sys.getenv("GOF_CORES", N_CORES))
N_CORES_GOF <- max(1L, min(N_CORES_GOF, N_GOF))
SEED_EVENTS_GOF <- as.integer(Sys.getenv("SEED_EVENTS_GOF", 20L))

# Data shaping
USE_FIRST_CONTACT_ONLY <- isTRUE(as.logical(Sys.getenv("USE_FIRST_CONTACT_ONLY", "TRUE")))
MAX_EDGES <- as.integer(Sys.getenv("MAX_EDGES", if (LOCAL_QUICK) 400L else 0L)) # 0 = no cap

# Background KDE
GRID_N <- as.integer(Sys.getenv("KDE_GRID_N", if (LOCAL_QUICK) 1024L else 4096L))
BW <- suppressWarnings(as.numeric(Sys.getenv("KDE_BW", NA_real_))) # NA => default

# Model specification (transitivity + degree)
FORMULA_RHS <- Sys.getenv("FORMULA_RHS", "edges + gwdegree(0.1) + gwesp(0.1)")
MARK_DECAY <- Sys.getenv("MARK_DECAY", "activity")
GROWTH_ONLY <- FALSE

cat("=== Hypertext conference (modern) ===\n")
cat("  LOCAL_QUICK:", LOCAL_QUICK, "\n")
cat("  cores (fit):", N_CORES, "| cores (gof):", N_CORES_GOF, "\n")
cat("  MAX_ITER:", MAX_ITER, "| RUN_GOF:", RUN_GOF, "| N_GOF:", N_GOF, "\n")
cat("  FORMULA_RHS:", FORMULA_RHS, "\n")
cat("  USE_FIRST_CONTACT_ONLY:", USE_FIRST_CONTACT_ONLY, "\n\n")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
make_hypertext_net <- function(df, use_first_contact_only = TRUE, max_edges = 0L) {
  stopifnot(all(c("from", "to", "time") %in% names(df)))

  # Undirected: canonical ordering (from < to) and drop duplicates
  swap <- df$from > df$to
  df[swap, c("from", "to")] <- df[swap, c("to", "from")]
  df <- df %>% distinct()

  # Time: shift to start at 0 and convert to hours (more stable than seconds)
  df$time <- as.numeric(df$time)
  df <- df[is.finite(df$time), ]
  df$time <- df$time - min(df$time)
  df$time <- df$time / 3600

  if (use_first_contact_only) {
    df <- df %>%
      dplyr::group_by(.data$from, .data$to) %>%
      dplyr::summarise(time = min(.data$time), .groups = "drop")
  }

  df <- df %>% arrange(time)

  if (max_edges > 0L && nrow(df) > max_edges) {
    df <- df[seq_len(max_edges), , drop = FALSE]
  }

  # Map vertex ids to 1..n
  nodes <- sort(unique(c(df$from, df$to)))
  id_map <- setNames(seq_along(nodes), nodes)
  df$tail <- unname(id_map[as.character(df$from)])
  df$head <- unname(id_map[as.character(df$to)])

  # Node "entry" times: first appearance in any edge
  node_time <- rep(Inf, length(nodes))
  idx_tail <- split(df$time, df$tail)
  idx_head <- split(df$time, df$head)
  for (i in seq_along(node_time)) {
    tt <- idx_tail[[as.character(i)]]
    th <- idx_head[[as.character(i)]]
    node_time[i] <- min(c(tt, th), na.rm = TRUE)
  }
  node_time[!is.finite(node_time)] <- max(df$time)

  # Build a simple undirected network
  el <- as.matrix(df[, c("tail", "head")])
  net <- network::network(el, matrix.type = "edgelist", directed = FALSE)
  network::set.edge.attribute(net, "time", df$time)
  network::set.vertex.attribute(net, "time", node_time)
  network::delete.vertex.attribute(net, "vertex.names")

  list(net = net, edges = df, nodes = nodes, id_map = id_map)
}

diagnose_truncation <- function(net, trunc_grid = c(25L, 50L, 100L, 200L, 400L), mark_decay = "activity") {
  nv <- network::network.size(net)
  tmax <- max(get_times(net)$times)
  net_final <- filtration_to_net(net, tmax, equals = TRUE)
  old_nodes <- nv
  new_nodes <- nv

  # Version-tolerant call: some installs may not export get_truncated_candidates or accept growth_only.
  gtc <- if (exists("get_truncated_candidates", mode = "function")) get_truncated_candidates else hawkesNet:::get_truncated_candidates
  gtc_formals <- names(formals(gtc))

  out <- data.frame(truncation = trunc_grid, n_candidate_edges = NA_integer_)
  for (k in seq_along(trunc_grid)) {
    tr <- as.integer(trunc_grid[k])
    args <- list(net = net_final, new_nodes = new_nodes, old_nodes = old_nodes,
                 truncation = tr, mark_decay = mark_decay)
    if ("growth_only" %in% gtc_formals) args$growth_only <- FALSE
    cands <- do.call(gtc, args)
    out$n_candidate_edges[k] <- length(cands$heads)
  }
  out
}

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
raw <- read.table(system.file("extdata", "ht09_contact_list.dat", package = "hawkesNet"))
df <- data.frame(
  time = raw$V1 / 20,   # original timestamps are in 1/20 seconds
  from = raw$V2,
  to   = raw$V3
)

# Small jitter to break ties (helps optimization / KDE)
set.seed(1)
df$time <- df$time + rnorm(nrow(df), 0, 0.01)

obj <- make_hypertext_net(df, use_first_contact_only = USE_FIRST_CONTACT_ONLY, max_edges = MAX_EDGES)
net <- obj$net

times <- get_times(net)$times
time_window <- c(min(times), max(times))
cat("  Network:", length(times), "event times | nodes:", network.size(net), "| edges:", network.edgecount(net), "\n")
cat("  Time window (hours):", sprintf("[%.3f, %.3f]", time_window[1], time_window[2]), "\n\n")

# -----------------------------------------------------------------------------
# Choose truncation (diagnostic + default)
# -----------------------------------------------------------------------------
if (is.na(TRUNCATION)) {
  diag <- diagnose_truncation(net, mark_decay = MARK_DECAY)
  cat("--- Truncation diagnostic (candidate edge counts) ---\n")
  print(diag, row.names = FALSE)
  # Default heuristic: keep candidate set in the low tens of thousands
  # (depends on node count; safe choice for local quick is smaller).
  TRUNCATION <- if (LOCAL_QUICK) 50L else min(200L, network.size(net))
  cat("  Using TRUNCATION =", TRUNCATION, "\n\n")
} else {
  cat("  Using TRUNCATION (from env) =", TRUNCATION, "\n\n")
}

# -----------------------------------------------------------------------------
# Estimate inhomogeneous background (KDE) over the full multi-day window
# -----------------------------------------------------------------------------
cat("--- Estimating inhomogeneous background (KDE) ---\n")
t_bg <- proc.time()
inhom_bg <- prepare_inhomogeneous_background(
  net,
  time_attr = "time",
  bw = if (is.finite(BW)) BW else NULL,
  grid_n = GRID_N
)
cat("  KDE done in", round((proc.time() - t_bg)[3], 2), "s\n\n")

# -----------------------------------------------------------------------------
# Fit CS model (non-growth, includes transitivity)
# -----------------------------------------------------------------------------
cat("--- Fitting hawkesNet (CS, non-growth, transitivity) ---\n")
exp_cs <- expected_params_PMF_mark_CS(net, FORMULA_RHS)
n_cs <- if (!is.na(exp_cs$CS_params_length)) exp_cs$CS_params_length else 3L
cat("  CS_params length:", n_cs, "\n")

mu_init <- inhom_bg$integral_bg / (time_window[2] - time_window[1])
params_init <- list(
  mu = mu_init,
  beta_overall = 0.3,
  K = 0.5,
  beta_edges = 0.3,
  node_lambda = 0.1,
  CS_params = c(-10, rep(0, n_cs - 1))
)

p_scale <- c(
  beta_overall = 0.1,
  K = 0.1,
  beta_edges = 0.1,
  node_lambda = 0.5,
  setNames(rep(0.1, n_cs), paste0("CS_params", seq_len(n_cs)))
)

fit <- fit_hawkesNet(
  params_init = params_init,
  time_window = time_window,
  mark_filtration = net,
  PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec,
  integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_RHS,
  truncation = TRUNCATION,
  mark_decay = MARK_DECAY,
  growth_only = GROWTH_ONLY,
  max_node_time = max(get_times(net)$node_times),
  method = "Nelder-Mead",
  maxit = MAX_ITER,
  reltol = 1e-8,
  trace = 0,
  verbose = TRUE,
  fixed_params = c("mu"),  # mu absorbed by inhomogeneous background
  parscale = p_scale,
  cache_intensity = TRUE,
  combine_intensity = TRUE,
  cores = N_CORES
)

if (!is.null(fit$fit_table)) {
  cat("\n--- Fit results ---\n")
  print(fit$fit_table, max = NULL)
}

# -----------------------------------------------------------------------------
# GOF (fast checks only)
# -----------------------------------------------------------------------------
gof_res <- NULL
if (RUN_GOF) {
  cat("\n--- GOF ---\n")
  gof_res <- gof(
    fit = fit,
    net_obs = net,
    params_init = params_init,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = FORMULA_RHS,
    time_window = time_window,
    truncation = TRUNCATION,
    mark_decay = MARK_DECAY,
    growth_only = GROWTH_ONLY,
    max_node_time = max(get_times(net)$node_times),
    inhom_bg = inhom_bg,
    n_sim = N_GOF,
    cores = N_CORES_GOF,
    max_deg = 15,
    k_esp = 15,
    degree = 0,
    esp = 0,
    mu_multiplier = 5,
    seed_events = SEED_EVENTS_GOF,
    verbose = TRUE
  )
  if (!is.null(gof_res$plots) && requireNamespace("ggplot2", quietly = TRUE)) {
    if (!is.null(gof_res$plots$degree_plot)) print(gof_res$plots$degree_plot)
    if (!is.null(gof_res$plots$esp_plot)) print(gof_res$plots$esp_plot)
    if (!is.null(gof_res$plots$geodist_plot)) print(gof_res$plots$geodist_plot)
    if (!is.null(gof_res$plots$waiting_times_plot)) print(gof_res$plots$waiting_times_plot)
  }
}

# -----------------------------------------------------------------------------
# Save results for poking
# -----------------------------------------------------------------------------
out_path <- file.path(PKG_ROOT, "inst", "hypertext_conference", "hypertext_conference_modern_results.rds")
saveRDS(
  list(
    fit = fit,
    gof = gof_res,
    net = net,
    edges = obj$edges,
    inhom_bg = inhom_bg,
    params_init = params_init,
    formula_rhs = FORMULA_RHS,
    truncation = TRUNCATION,
    time_window = time_window
  ),
  out_path
)
cat("\nSaved:", out_path, "\n")
cat("Done.\n")

