## =============================================================================
## Hypertext 2009 conference: modern hawkesNet fit + GOF
## =============================================================================
## Run from package root:
##   Rscript inst/hypertext_conference/hypertext_conference_modern_fit.R
##
## Quick local mode (keeps runtime under ~2 minutes):
##   LOCAL_QUICK=TRUE Rscript inst/hypertext_conference/hypertext_conference_modern_fit.R
##
## Cluster mode (SLURM):
##   sbatch inst/hypertext_conference/run_hypertext.slurm
##
## Full run (default under SLURM): MAX_EDGES=10000, no truncation, N_GOF=100, MAX_ITER=5000, 100 cores.
## Cluster mode knobs:
##   SLURM_CPUS_PER_TASK=32 MAX_ITER=5000 N_GOF=50 GOF_CORES=32 Rscript ...
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

# Paths: run from package root. Under SLURM, use submit dir so path stays valid.
PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) Sys.getenv("SLURM_SUBMIT_DIR") else getwd()
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) {
  PKG_ROOT <- getwd()
}

# Use the repo (most up-to-date) implementation even if the installed package is older.
# This is helpful during development, but note: on macOS/Windows we use PSOCK workers for parallelism,
# and PSOCK workers do NOT automatically see functions you sourced into the master session.
# For robust parallel runs, prefer installing the package and leave USE_REPO_CODE=FALSE.
USE_REPO_CODE <- isTRUE(as.logical(Sys.getenv("USE_REPO_CODE", "FALSE")))
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

# Default to 7 cores even locally (override with SLURM_CPUS_PER_TASK / env var if desired).
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 7L))
N_CORES <- max(1L, N_CORES)

# If we're sourcing repo code in a PSOCK-only context (macOS/Windows/interactive), force serial
# to avoid "could not find function ..." errors inside workers.
if (USE_REPO_CODE) {
  os <- Sys.info()[["sysname"]]
  psock_only <- (os %in% c("Darwin", "Windows")) || interactive() || isTRUE(getOption("hawkesNet.force_psock", FALSE))
  if (psock_only && N_CORES > 1L) {
    message("NOTE: USE_REPO_CODE=TRUE under PSOCK; forcing N_CORES=1 for worker consistency. ",
            "Install hawkesNet and set USE_REPO_CODE=FALSE to use multiple cores.")
    N_CORES <- 1L
  }
}

MAX_ITER <- as.integer(Sys.getenv("MAX_ITER", if (LOCAL_QUICK) 200L else 5000L))
trunc_env <- Sys.getenv("TRUNCATION", "")
TRUNCATION <- if (nzchar(trunc_env)) suppressWarnings(as.integer(trunc_env)) else NA_integer_
if (length(TRUNCATION) != 1L || !is.finite(TRUNCATION)) TRUNCATION <- NA_integer_

# GOF: proper run uses 100 sims (matches 100 cores); quick run skips or uses 2
RUN_GOF <- isTRUE(as.logical(Sys.getenv("RUN_GOF", if (LOCAL_QUICK) "FALSE" else "TRUE")))
N_GOF <- as.integer(Sys.getenv("N_GOF", if (LOCAL_QUICK) 2L else 100L))
# Optional quick simulation test from fitted values (before GOF). Skip when under SLURM (saves time).
RUN_SIM_TEST <- isTRUE(as.logical(Sys.getenv("RUN_SIM_TEST", if (LOCAL_QUICK) "TRUE" else "FALSE")))
N_SIM_TEST <- as.integer(Sys.getenv("N_SIM_TEST", 3L))
N_GOF <- max(1L, N_GOF)
N_CORES_GOF <- as.integer(Sys.getenv("GOF_CORES", N_CORES))
N_CORES_GOF <- max(1L, min(N_CORES_GOF, N_GOF))
# PSOCK outer workers for GOF (avoids fork deadlocks; faster on 100 cores)
N_GOF_OUTER <- as.integer(Sys.getenv("GOF_CORES_OUTER", if (N_CORES >= 32L) min(50L, N_CORES) else 0L))
SEED_EVENTS_GOF <- as.integer(Sys.getenv("SEED_EVENTS_GOF", 20L))

# Data shaping
USE_FIRST_CONTACT_ONLY <- isTRUE(as.logical(Sys.getenv("USE_FIRST_CONTACT_ONLY", "TRUE")))
MAX_EDGES <- as.integer(Sys.getenv("MAX_EDGES", if (LOCAL_QUICK) 400L else 10000L))

# Background KDE
GRID_N <- as.integer(Sys.getenv("KDE_GRID_N", if (LOCAL_QUICK) 1024L else 4096L))
BW <- suppressWarnings(as.numeric(Sys.getenv("KDE_BW", NA_real_))) # NA => default

# Model specification (transitivity + degree).
# Include "edges" as baseline for interpretability of other CS stats. Fix K for stability.
# degree(0) penalizes degree-0 nodes (reduces excess isolates in simulations).
FORMULA_RHS <- Sys.getenv("FORMULA_RHS", "edges + degree(0) + gwdegree(0.1) + gwesp(0.1)")
FORMULA_RHS_STAR_ESP <- "edges + degree(0) + star(c(2,3,4)) + esp(2:4)"  # Alternative: star+esp instead of gwdegree+gwesp
FORMULA_RHS_DECAY001 <- "edges + degree(0) + gwdegree(0.01) + gwesp(0.01)"  # Smaller decay = stronger penalty on fat tails
MARK_DECAY <- Sys.getenv("MARK_DECAY", "activity")
GROWTH_ONLY <- FALSE

cat("=== Hypertext conference (modern) ===\n")
cat("  LOCAL_QUICK:", LOCAL_QUICK, "\n")
cat("  cores (fit):", N_CORES, "| cores (gof):", N_CORES_GOF, "| gof_outer:", N_GOF_OUTER, "\n")
cat("  MAX_ITER:", MAX_ITER, "| RUN_GOF:", RUN_GOF, "| N_GOF:", N_GOF, "| RUN_SIM_TEST:", RUN_SIM_TEST, "\n")
cat("  FORMULA_RHS:", FORMULA_RHS, "\n")
cat("  USE_FIRST_CONTACT_ONLY:", USE_FIRST_CONTACT_ONLY, "| MAX_EDGES:", if (MAX_EDGES > 0) MAX_EDGES else "no cap", "\n\n")

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

  # Map vertex ids to 1..n in *entry-time* order so node_entrance truncation is meaningful.
  # Entry time = first appearance in any edge.
  nodes_raw <- sort(unique(c(df$from, df$to)))
  entry_time_by_node <- vapply(nodes_raw, function(v) {
    min(df$time[df$from == v | df$to == v], na.rm = TRUE)
  }, numeric(1))
  ord_nodes <- order(entry_time_by_node, nodes_raw)
  nodes <- nodes_raw[ord_nodes]
  entry_time_by_node <- entry_time_by_node[ord_nodes]

  id_map <- setNames(seq_along(nodes), nodes)
  df$tail <- unname(id_map[as.character(df$from)])
  df$head <- unname(id_map[as.character(df$to)])

  # Node "entry" times in ID order
  node_time <- as.numeric(entry_time_by_node)

  # Build a simple undirected network
  el <- as.matrix(df[, c("tail", "head")])
  net <- network::network(el, matrix.type = "edgelist", directed = FALSE)
  network::set.edge.attribute(net, "time", df$time)
  network::set.vertex.attribute(net, "time", node_time)
  network::delete.vertex.attribute(net, "vertex.names")

  list(net = net, edges = df, nodes = nodes, id_map = id_map)
}

diagnose_truncation <- function(edges_df, node_entry_time,
                                trunc_grid = c(25L, 50L, 100L, 200L, 400L),
                                mark_decay = c("node_entrance", "activity"),
                                n_nodes = NULL) {
  if (!is.null(n_nodes)) trunc_grid <- sort(unique(c(trunc_grid, n_nodes)))
  mark_decay <- match.arg(mark_decay)
  stopifnot(all(c("tail", "head", "time") %in% names(edges_df)))
  stopifnot(is.numeric(node_entry_time) && length(node_entry_time) >= max(edges_df$tail, edges_df$head))

  edges_df <- edges_df[order(edges_df$time), , drop = FALSE]
  n_all <- length(node_entry_time)

  # Edges that are *impossible regardless of truncation*: endpoint not yet present at its first-contact time.
  # (Also catches inconsistent node-time assignment.)
  impossible_any <- edges_df$time < pmax(node_entry_time[edges_df$tail], node_entry_time[edges_df$head])
  # In a simple graph, re-adding an existing edge is impossible under any truncation.
  u <- pmin(edges_df$tail, edges_df$head)
  v <- pmax(edges_df$tail, edges_df$head)
  edge_key <- u + n_all * v
  impossible_any <- impossible_any | duplicated(edge_key)
  n_impossible_any <- sum(impossible_any, na.rm = TRUE)

  # For "activity" truncation: maintain last-activity times as we replay events.
  last_activity <- node_entry_time
  n_edges <- nrow(edges_df)
  n_impossible <- integer(length(trunc_grid))
  n_candidate_edges <- numeric(length(trunc_grid))

  # Precompute, per edge, the number of present nodes at that time (entry time <= t)
  # Since node IDs are entry-ordered, we can advance a pointer.
  entry_sorted <- node_entry_time
  if (!isTRUE(all(diff(entry_sorted) >= -1e-12))) {
    # Safety: if entry times are not monotone, fall back to ranking by time
    entry_sorted <- sort(entry_sorted)
  }
  present_count_at_edge <- integer(n_edges)
  ptr <- 0L
  for (i in seq_len(n_edges)) {
    t_i <- edges_df$time[i]
    while (ptr < n_all && entry_sorted[ptr + 1L] <= t_i) ptr <- ptr + 1L
    present_count_at_edge[i] <- ptr
  }

  # Helper: compute active set (IDs) for current time given truncation
  get_active_nodes <- function(n_present, truncation) {
    if (n_present <= 1L) return(integer(0))
    if (mark_decay == "node_entrance") {
      window_start <- max(1L, n_present - truncation + 1L)
      return(window_start:n_present)
    }
    # activity: top-k by last_activity among present nodes; tie-break by id
    ids <- seq_len(n_present)
    ord <- order(last_activity[ids], ids, decreasing = TRUE)
    sort(ord[seq_len(min(truncation, n_present))])
  }

  # Main pass for each truncation value.
  # Optimized enough for this dataset size; avoids constructing networks or ERNM objects.
  for (k in seq_along(trunc_grid)) {
    tr <- as.integer(trunc_grid[k])
    if (!is.finite(tr) || tr < 2L) tr <- 2L

    # Reset activity tracker for each truncation run.
    last_activity <- node_entry_time
    imp <- 0L
    cand_total <- 0

    for (i in seq_len(n_edges)) {
      if (impossible_any[i]) {
        imp <- imp + 1L
        next
      }
      n_present <- present_count_at_edge[i]
      active <- get_active_nodes(n_present, tr)

      # Candidate edges count at this time (ignores existing-edge filtering; OK for rough scaling).
      cand_total <- cand_total + (length(active) * (length(active) - 1L) / 2)

      tail_i <- edges_df$tail[i]
      head_i <- edges_df$head[i]
      if (!(tail_i %in% active && head_i %in% active)) {
        imp <- imp + 1L
      }

      # Update activity after processing this edge (so "pre-edge" activity defines candidate set).
      if (mark_decay == "activity") {
        last_activity[tail_i] <- edges_df$time[i]
        last_activity[head_i] <- edges_df$time[i]
      }
    }

    n_impossible[k] <- imp
    n_candidate_edges[k] <- cand_total
  }

  data.frame(
    mark_decay = mark_decay,
    truncation = as.integer(trunc_grid),
    n_impossible_edges = n_impossible,
    prop_impossible = round(n_impossible / n_edges, 4),
    pct_impossible_edges = round(100 * n_impossible / n_edges, 1),
    n_candidate_edges_total = n_candidate_edges,
    n_impossible_any_truncation = n_impossible_any,
    pct_impossible_any_truncation = round(100 * n_impossible_any / n_edges, 1),
    stringsAsFactors = FALSE
  )
}

#' Verify that the chosen truncation makes all observed edges possible.
#' Stops with error if any edge would be impossible under the given truncation.
verify_truncation <- function(edges_df, node_entry_time, truncation, mark_decay) {
  diag <- diagnose_truncation(edges_df, node_entry_time,
                              trunc_grid = as.integer(truncation),
                              mark_decay = mark_decay)
  n_imp <- diag$n_impossible_edges[1]
  n_any <- diag$n_impossible_any_truncation[1]
  if (n_any > 0) {
    stop("Truncation verification failed: ", n_any, " edges are impossible under ANY truncation ",
         "(endpoint not yet present or duplicate edge). Check node entry times.")
  }
  if (n_imp > 0) {
    stop("Truncation verification failed: TRUNCATION=", truncation, " makes ", n_imp,
         " observed edges impossible. Increase TRUNCATION or use network.size(net) for no truncation.")
  }
  invisible(TRUE)
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
# Choose truncation
# -----------------------------------------------------------------------------
# Hypertext: truncation = number of nodes (all nodes in candidate set at all times).
n_nodes <- network.size(net)
if (is.na(TRUNCATION)) {
  TRUNCATION <- n_nodes
  cat("  Using TRUNCATION =", TRUNCATION, "(= network size; all edges possible)\n\n")
} else {
  cat("  Using TRUNCATION (from env) =", TRUNCATION, "\n\n")
}

# Verify: all edges must be possible under the chosen truncation
verify_truncation(obj$edges, node_entry_time = net %v% "time",
                  truncation = TRUNCATION, mark_decay = MARK_DECAY)

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
  CS_params = c(-8, -5, rep(0, n_cs - 2))  # edges, degree(0), then gwdegree/gwesp
)

p_scale <- c(
  beta_overall = 0.1,
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
  fixed_params = c("mu", "K"),  # mu absorbed by inhomogeneous background; K fixed for stability
  parscale = p_scale,
  cache_intensity = TRUE,
  combine_intensity = TRUE,
  cores = N_CORES
)

if (!is.null(fit$fit_table)) {
  cat("\n--- Fit results (gwdegree + gwesp) ---\n")
  print(fit$fit_table, max = NULL)
}

# -----------------------------------------------------------------------------
# Fit alternative: star(c(2,3,4)) + esp(2:4) instead of gwdegree + gwesp
# -----------------------------------------------------------------------------
cat("\n--- Fitting hawkesNet (CS, star+esp variant) ---\n")
exp_cs_alt <- expected_params_PMF_mark_CS(net, FORMULA_RHS_STAR_ESP)
n_cs_alt <- if (!is.na(exp_cs_alt$CS_params_length)) exp_cs_alt$CS_params_length else 7L
cat("  CS_params length:", n_cs_alt, "\n")

params_init_star_esp <- list(
  mu = mu_init,
  beta_overall = 0.3,
  K = 0.5,
  beta_edges = 0.3,
  node_lambda = 0.1,
  CS_params = c(-8, -5, rep(0, n_cs_alt - 2))  # edges, degree(0), then star/esp
)

p_scale_star_esp <- c(
  beta_overall = 0.1,
  beta_edges = 0.1,
  node_lambda = 0.5,
  setNames(rep(0.1, n_cs_alt), paste0("CS_params", seq_len(n_cs_alt)))
)

fit_star_esp <- fit_hawkesNet(
  params_init = params_init_star_esp,
  time_window = time_window,
  mark_filtration = net,
  PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec,
  integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_RHS_STAR_ESP,
  truncation = TRUNCATION,
  mark_decay = MARK_DECAY,
  growth_only = GROWTH_ONLY,
  max_node_time = max(get_times(net)$node_times),
  method = "Nelder-Mead",
  maxit = MAX_ITER,
  reltol = 1e-8,
  trace = 0,
  verbose = TRUE,
  fixed_params = c("mu", "K"),
  parscale = p_scale_star_esp,
  cache_intensity = TRUE,
  combine_intensity = TRUE,
  cores = N_CORES
)

if (!is.null(fit_star_esp$fit_table)) {
  cat("\n--- Fit results (star+esp) ---\n")
  print(fit_star_esp$fit_table, max = NULL)
}

# -----------------------------------------------------------------------------
# Fit alternative: decay=0.01 (stronger penalty on fat tails)
# -----------------------------------------------------------------------------
cat("\n--- Fitting hawkesNet (CS, decay=0.01) ---\n")
exp_cs_decay <- expected_params_PMF_mark_CS(net, FORMULA_RHS_DECAY001)
n_cs_decay <- if (!is.na(exp_cs_decay$CS_params_length)) exp_cs_decay$CS_params_length else 3L
cat("  CS_params length:", n_cs_decay, "\n")

params_init_decay001 <- list(
  mu = mu_init,
  beta_overall = 0.3,
  K = 0.5,
  beta_edges = 0.3,
  node_lambda = 0.1,
  CS_params = c(-8, -5, rep(0, n_cs_decay - 2))  # edges, degree(0), then gwdegree/gwesp
)

p_scale_decay001 <- c(
  beta_overall = 0.1,
  beta_edges = 0.1,
  node_lambda = 0.5,
  setNames(rep(0.1, n_cs_decay), paste0("CS_params", seq_len(n_cs_decay)))
)

fit_decay001 <- fit_hawkesNet(
  params_init = params_init_decay001,
  time_window = time_window,
  mark_filtration = net,
  PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec,
  integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_RHS_DECAY001,
  truncation = TRUNCATION,
  mark_decay = MARK_DECAY,
  growth_only = GROWTH_ONLY,
  max_node_time = max(get_times(net)$node_times),
  method = "Nelder-Mead",
  maxit = MAX_ITER,
  reltol = 1e-8,
  trace = 0,
  verbose = TRUE,
  fixed_params = c("mu", "K"),
  parscale = p_scale_decay001,
  cache_intensity = TRUE,
  combine_intensity = TRUE,
  cores = N_CORES
)

if (!is.null(fit_decay001$fit_table)) {
  cat("\n--- Fit results (decay=0.01) ---\n")
  print(fit_decay001$fit_table, max = NULL)
}

# -----------------------------------------------------------------------------
# Optional: quick simulation test from fitted values (before GOF)
# -----------------------------------------------------------------------------
sim_test_res <- NULL
if (RUN_SIM_TEST && !is.null(fit$params)) {
  cat("\n--- Simulation test (fitted model sanity check) ---\n")
  pfit <- fit$params
  seed_net_test <- NULL
  seed_times_test <- NULL
  if (SEED_EVENTS_GOF > 0) {
    all_times <- get_times(net)$times
    if (length(all_times) >= SEED_EVENTS_GOF) {
      t_seed <- all_times[SEED_EVENTS_GOF]
      seed_net_test <- filtration_to_net(net, t_seed, equals = TRUE)
      seed_times_test <- all_times[1:SEED_EVENTS_GOF]
      cat("  Seeding with first", SEED_EVENTS_GOF, "events (up to t =", round(t_seed, 4), ")\n")
    }
  }
  sims_test <- vector("list", N_SIM_TEST)
  for (i in seq_len(N_SIM_TEST)) {
    sims_test[[i]] <- tryCatch({
      sim_hawkesNet(
        params = pfit,
        time_window = time_window,
        PMF_mark = PMF_mark_CS,
        cond_intensity = cond_intensity,
        formula_RHS = FORMULA_RHS,
        truncation = TRUNCATION,
        mark_decay = MARK_DECAY,
        growth_only = GROWTH_ONLY,
        max_node_time = max(get_times(net)$node_times),
        hashed_edges = TRUE,
        verbose = FALSE,
        mu_multiplier = 5,
        stop_on_full_network = FALSE,
        inhom_bg = inhom_bg,
        seed_net = seed_net_test,
        seed_times = seed_times_test
      )
    }, error = function(e) {
      list(net = NULL, events = list(t = numeric(0)), error = e$message)
    })
  }
  n_obs <- network::network.size(net)
  e_obs <- network::network.edgecount(net)
  deg_obs <- sna::degree(net, gmode = "graph")
  mean_deg_obs <- mean(deg_obs[deg_obs > 0], na.rm = TRUE)
  if (is.na(mean_deg_obs) || !is.finite(mean_deg_obs)) mean_deg_obs <- 0
  n_sim_ok <- sum(!sapply(sims_test, function(s) is.null(s$net)))
  n_edges_sim <- vapply(sims_test, function(s) {
    if (is.null(s$net)) NA_integer_ else network::network.edgecount(s$net)
  }, integer(1))
  n_nodes_sim <- vapply(sims_test, function(s) {
    if (is.null(s$net)) NA_integer_ else network::network.size(s$net)
  }, integer(1))
  mean_deg_sim <- vapply(sims_test, function(s) {
    if (is.null(s$net)) NA_real_
    else {
      d <- sna::degree(s$net, gmode = "graph")
      m <- mean(d[d > 0], na.rm = TRUE)
      if (is.na(m) || !is.finite(m)) 0 else m
    }
  }, numeric(1))
  cat("  Observed:  n_edges =", e_obs, "| n_nodes =", n_obs, "| mean_deg(>0) =", round(mean_deg_obs, 2), "\n")
  cat("  Simulated:", n_sim_ok, "/", N_SIM_TEST, "succeeded\n")
  if (n_sim_ok > 0) {
    cat("    n_edges:  ", paste(na.omit(n_edges_sim), collapse = ", "), "\n")
    cat("    n_nodes:  ", paste(na.omit(n_nodes_sim), collapse = ", "), "\n")
    cat("    mean_deg: ", paste(round(na.omit(mean_deg_sim), 2), collapse = ", "), "\n")
  }
  sim_test_res <- list(sims = sims_test, n_obs = n_obs, e_obs = e_obs, mean_deg_obs = mean_deg_obs)
}

# -----------------------------------------------------------------------------
# GOF (fast checks only)
# -----------------------------------------------------------------------------
gof_res <- NULL
if (RUN_GOF) {
  cat("\n--- GOF ---\n")
  # Version-tolerant GOF call: older hawkesNet installs may not accept newer args
  # (e.g., seed_events, growth_only). Also force-evaluate script-level symbols
  # so PSOCK workers don't see unevaluated promises like `GROWTH_ONLY`.
  gof_fun <- if (exists("gof", mode = "function")) get("gof", mode = "function") else hawkesNet::gof
  gof_formals <- names(formals(gof_fun))
  gof_args <- list(
    fit = fit,
    net_obs = net,
    params_init = params_init,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = FORMULA_RHS,
    time_window = time_window,
    truncation = as.integer(TRUNCATION),
    mark_decay = as.character(MARK_DECAY),
    growth_only = isTRUE(GROWTH_ONLY),
    max_node_time = max(get_times(net)$node_times),
    inhom_bg = inhom_bg,
    n_sim = as.integer(N_GOF),
    cores = as.integer(N_CORES_GOF),
    cores_outer = if (N_GOF_OUTER > 0L) as.integer(N_GOF_OUTER) else NULL,
    max_deg = 30L,
    k_esp = 30L,
    degree = 0L,
    esp = 0L,
    mu_multiplier = 5,
    seed_events = as.integer(SEED_EVENTS_GOF),
    verbose = TRUE
  )
  gof_args <- gof_args[names(gof_args) %in% gof_formals]
  gof_res <- do.call(gof_fun, gof_args)
  if (!is.null(gof_res$plots) && requireNamespace("ggplot2", quietly = TRUE)) {
    if (!is.null(gof_res$plots$degree_plot)) print(gof_res$plots$degree_plot)
    if (!is.null(gof_res$plots$esp_plot)) print(gof_res$plots$esp_plot)
    if (!is.null(gof_res$plots$geodist_plot)) print(gof_res$plots$geodist_plot)
    if (!is.null(gof_res$plots$waiting_times_plot)) print(gof_res$plots$waiting_times_plot)
  }
}

# -----------------------------------------------------------------------------
# Save results (cluster_output when on SLURM, same structure as openalex)
# -----------------------------------------------------------------------------
# GOF for star+esp fit (same settings, degree/ESP up to 30)
gof_star_esp <- NULL
gof_decay001 <- NULL
if (RUN_GOF && !is.null(fit_star_esp$fit)) {
  cat("\n--- GOF (star+esp fit) ---\n")
  gof_args_alt <- list(
    fit = fit_star_esp,
    net_obs = net,
    params_init = params_init_star_esp,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = FORMULA_RHS_STAR_ESP,
    time_window = time_window,
    truncation = as.integer(TRUNCATION),
    mark_decay = as.character(MARK_DECAY),
    growth_only = isTRUE(GROWTH_ONLY),
    max_node_time = max(get_times(net)$node_times),
    inhom_bg = inhom_bg,
    n_sim = as.integer(N_GOF),
    cores = as.integer(N_CORES_GOF),
    cores_outer = if (N_GOF_OUTER > 0L) as.integer(N_GOF_OUTER) else NULL,
    max_deg = 30L,
    k_esp = 30L,
    degree = 0L,
    esp = 0L,
    mu_multiplier = 5,
    seed_events = as.integer(SEED_EVENTS_GOF),
    verbose = TRUE
  )
  gof_args_alt <- gof_args_alt[names(gof_args_alt) %in% gof_formals]
  gof_star_esp <- do.call(gof_fun, gof_args_alt)
}

# GOF for decay=0.01 fit
if (RUN_GOF && !is.null(fit_decay001$fit)) {
  cat("\n--- GOF (decay=0.01 fit) ---\n")
  gof_args_decay <- list(
    fit = fit_decay001,
    net_obs = net,
    params_init = params_init_decay001,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = FORMULA_RHS_DECAY001,
    time_window = time_window,
    truncation = as.integer(TRUNCATION),
    mark_decay = as.character(MARK_DECAY),
    growth_only = isTRUE(GROWTH_ONLY),
    max_node_time = max(get_times(net)$node_times),
    inhom_bg = inhom_bg,
    n_sim = as.integer(N_GOF),
    cores = as.integer(N_CORES_GOF),
    cores_outer = if (N_GOF_OUTER > 0L) as.integer(N_GOF_OUTER) else NULL,
    max_deg = 30L,
    k_esp = 30L,
    degree = 0L,
    esp = 0L,
    mu_multiplier = 5,
    seed_events = as.integer(SEED_EVENTS_GOF),
    verbose = TRUE
  )
  gof_args_decay <- gof_args_decay[names(gof_args_decay) %in% gof_formals]
  gof_decay001 <- do.call(gof_fun, gof_args_decay)
}

save_list <- list(
  fit = fit,
  fit_star_esp = fit_star_esp,
  fit_decay001 = fit_decay001,
  gof = gof_res,
  gof_star_esp = gof_star_esp,
  gof_decay001 = gof_decay001,
  GOF = gof_res,  # alias for consistency with openalex
  sim_test = sim_test_res,
  net = net,
  edges = obj$edges,
  inhom_bg = inhom_bg,
  params_init = params_init,
  params_init_star_esp = params_init_star_esp,
  params_init_decay001 = params_init_decay001,
  formula_rhs = FORMULA_RHS,
  formula_rhs_star_esp = FORMULA_RHS_STAR_ESP,
  formula_rhs_decay001 = FORMULA_RHS_DECAY001,
  truncation = TRUNCATION,
  time_window = time_window,
  N_GOF = N_GOF,
  SEED_EVENTS_GOF = SEED_EVENTS_GOF
)

cluster_output_dir <- file.path(PKG_ROOT, "cluster_output")
rds_path_primary <- file.path(cluster_output_dir, "results_hypertext_full.RDS")

# Use cluster_output when running under SLURM or when it exists
use_cluster_output <- nzchar(Sys.getenv("SLURM_JOB_ID")) || dir.exists(cluster_output_dir)

if (use_cluster_output) {
  dir.create(cluster_output_dir, showWarnings = FALSE, recursive = TRUE)
  save_ok <- tryCatch({
    saveRDS(save_list, rds_path_primary)
    cat("\nSaved to cluster_output:", rds_path_primary, "\n")
    TRUE
  }, error = function(e) {
    cat("  Save to ", rds_path_primary, " failed: ", conditionMessage(e), "\n")
    FALSE
  })
  if (!save_ok) {
    out_path <- file.path(PKG_ROOT, "inst", "hypertext_conference", "hypertext_conference_modern_results.rds")
    tryCatch({
      saveRDS(save_list, out_path)
      cat("  Fallback save:", out_path, "\n")
    }, error = function(e) cat("  Fallback save failed:", conditionMessage(e), "\n"))
  }
} else {
  out_path <- file.path(PKG_ROOT, "inst", "hypertext_conference", "hypertext_conference_modern_results.rds")
  saveRDS(save_list, out_path)
  cat("\nSaved:", out_path, "\n")
}
cat("Done.\n")

