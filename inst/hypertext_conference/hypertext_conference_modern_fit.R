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
## Two fits:
##   Fit 1 (primary) — First session only (before first overnight gap).
##                      Clean single-session data, ~946 edges, ~100 nodes.
##   Fit 2 (full)    — All three days, with mu forced to zero during overnight
##                      gaps (>1 h between events). The KDE background is estimated
##                      from active-period events only and zeroed out in gap intervals.
##
## Full run (default under SLURM): both fits, GOF for Fit 1 only.
## Cluster mode knobs:
##   SLURM_CPUS_PER_TASK=100 MAX_ITER=5000 N_GOF=100 GOF_CORES=100
##   RUN_GOF_FULL=TRUE  (to also run GOF on the full-data fit)
##
## Notes:
## - This models *edge formation* (simple graph): we collapse repeated contacts to first contact per dyad.
## - We do NOT use growth_only; transitivity is included via gwesp.
## - Inhomogeneous background is estimated via KDE.
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
# (Removed USE_REPO_CODE block - always use installed package)

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
default_local_quick <- if (nzchar(Sys.getenv("SLURM_JOB_ID")) || nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) "FALSE" else "TRUE"
LOCAL_QUICK <- isTRUE(as.logical(Sys.getenv("LOCAL_QUICK", default_local_quick)))

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 7L))
N_CORES <- max(1L, N_CORES)

MAX_ITER <- as.integer(Sys.getenv("MAX_ITER", if (LOCAL_QUICK) 200L else 5000L))
trunc_env <- Sys.getenv("TRUNCATION", "")
TRUNCATION <- if (nzchar(trunc_env)) suppressWarnings(as.integer(trunc_env)) else NA_integer_
if (length(TRUNCATION) != 1L || !is.finite(TRUNCATION)) TRUNCATION <- NA_integer_

# GOF controls: separate flags for day-1 and full fits.
# Day-1 GOF runs by default under SLURM; full GOF is opt-in.
RUN_GOF_DAY1 <- isTRUE(as.logical(Sys.getenv("RUN_GOF_DAY1",
                        Sys.getenv("RUN_GOF", if (LOCAL_QUICK) "FALSE" else "TRUE"))))
RUN_GOF_FULL <- isTRUE(as.logical(Sys.getenv("RUN_GOF_FULL", "FALSE")))
N_GOF <- as.integer(Sys.getenv("N_GOF", if (LOCAL_QUICK) 2L else 100L))
N_GOF <- max(1L, N_GOF)
N_CORES_GOF <- as.integer(Sys.getenv("GOF_CORES", N_CORES))
N_CORES_GOF <- max(1L, min(N_CORES_GOF, N_GOF))
N_GOF_OUTER <- as.integer(Sys.getenv("GOF_CORES_OUTER", if (N_CORES >= 32L) min(50L, N_CORES) else 0L))
SEED_EVENTS_GOF <- as.integer(Sys.getenv("SEED_EVENTS_GOF", 20L))

RUN_SIM_TEST <- isTRUE(as.logical(Sys.getenv("RUN_SIM_TEST", if (LOCAL_QUICK) "TRUE" else "FALSE")))
N_SIM_TEST <- as.integer(Sys.getenv("N_SIM_TEST", 3L))

# Run Fit 2 (full data, mu=0 in gaps)? Default: TRUE under SLURM, FALSE locally.
RUN_FULL_FIT <- isTRUE(as.logical(Sys.getenv("RUN_FULL_FIT", if (LOCAL_QUICK) "FALSE" else "TRUE")))

# Data shaping
USE_FIRST_CONTACT_ONLY <- isTRUE(as.logical(Sys.getenv("USE_FIRST_CONTACT_ONLY", "TRUE")))
MAX_EDGES <- as.integer(Sys.getenv("MAX_EDGES", 0L))  # 0 = no cap

# Background KDE
GRID_N <- as.integer(Sys.getenv("KDE_GRID_N", if (LOCAL_QUICK) 1024L else 4096L))
BW <- suppressWarnings(as.numeric(Sys.getenv("KDE_BW", NA_real_)))

# Gap threshold (hours) for identifying overnight breaks
GAP_THRESHOLD <- 1.0

# Model specification
FORMULA_RHS <- Sys.getenv("FORMULA_RHS", "edges + degree(0) + gwdegree(0.1) + gwesp(0.1)")
MARK_DECAY <- Sys.getenv("MARK_DECAY", "activity")
GROWTH_ONLY <- FALSE

cat("=== Hypertext conference (modern) ===\n")
cat("  LOCAL_QUICK:", LOCAL_QUICK, "\n")
cat("  cores (fit):", N_CORES, "| cores (gof):", N_CORES_GOF, "| gof_outer:", N_GOF_OUTER, "\n")
cat("  MAX_ITER:", MAX_ITER, "\n")
cat("  RUN_GOF_DAY1:", RUN_GOF_DAY1, "| RUN_GOF_FULL:", RUN_GOF_FULL, "| N_GOF:", N_GOF, "\n")
cat("  RUN_FULL_FIT:", RUN_FULL_FIT, "| RUN_SIM_TEST:", RUN_SIM_TEST, "\n")
cat("  FORMULA_RHS:", FORMULA_RHS, "\n")
cat("  USE_FIRST_CONTACT_ONLY:", USE_FIRST_CONTACT_ONLY, "| MAX_EDGES:", if (MAX_EDGES > 0) MAX_EDGES else "no cap", "\n\n")

# =============================================================================
# Helpers
# =============================================================================

#' Build a network object from hypertext edge data.
#' No timeline compression — times are in hours from first event.
make_hypertext_net <- function(df, use_first_contact_only = TRUE, max_edges = 0L) {
  stopifnot(all(c("from", "to", "time") %in% names(df)))

  # Undirected: canonical ordering (from < to) and drop duplicates
  swap <- df$from > df$to
  df[swap, c("from", "to")] <- df[swap, c("to", "from")]
  df <- df %>% distinct()

  # Time: already in hours and shifted to 0
  df$time <- as.numeric(df$time)
  df <- df[is.finite(df$time), ]
  df <- df[order(df$time), , drop = FALSE]

  if (use_first_contact_only) {
    df <- df %>%
      dplyr::group_by(.data$from, .data$to) %>%
      dplyr::summarise(time = min(.data$time), .groups = "drop")
  }

  df <- df %>% arrange(time)

  if (max_edges > 0L && nrow(df) > max_edges) {
    df <- df[seq_len(max_edges), , drop = FALSE]
  }

  # Map vertex ids to 1..n in *entry-time* order
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

  node_time <- as.numeric(entry_time_by_node)

  el <- as.matrix(df[, c("tail", "head")])
  net <- network::network(el, matrix.type = "edgelist", directed = FALSE)
  network::set.edge.attribute(net, "time", df$time)
  network::set.vertex.attribute(net, "time", node_time)
  network::delete.vertex.attribute(net, "vertex.names")

  list(net = net, edges = df, nodes = nodes, id_map = id_map)
}

#' Identify gap intervals (overnight breaks) from event times.
#' Returns a data.frame with columns: start, end (in hours).
find_gap_intervals <- function(event_times, gap_threshold = 1.0) {
  event_times <- sort(event_times)
  gaps <- diff(event_times)
  big <- which(gaps > gap_threshold)
  if (length(big) == 0) return(data.frame(start = numeric(0), end = numeric(0)))
  data.frame(
    start = event_times[big],
    end   = event_times[big + 1L]
  )
}

#' Check whether a time falls inside any gap interval.
in_gap <- function(t, gap_intervals) {
  if (nrow(gap_intervals) == 0) return(rep(FALSE, length(t)))
  out <- logical(length(t))
  for (i in seq_len(nrow(gap_intervals))) {
    out <- out | (t > gap_intervals$start[i] & t < gap_intervals$end[i])
  }
  out
}

#' Zero out the KDE background during gap intervals and recompute integral_bg.
#'
#' Takes the output of prepare_inhomogeneous_background() and a gap_intervals
#' data.frame (from find_gap_intervals), zeroes the KDE rate in those intervals,
#' and returns a modified inhom_bg with corrected mu_vec, integral_bg, mu_fun.
zero_gaps_inhom_bg <- function(inhom_bg, gap_intervals) {
  if (nrow(gap_intervals) == 0) return(inhom_bg)

  mu_fit <- inhom_bg$mu_fit
  grid   <- mu_fit$grid
  mu_grid <- mu_fit$mu_grid

  # Zero out grid points inside gaps
  is_gap <- in_gap(grid, gap_intervals)
  mu_grid_zeroed <- mu_grid
  mu_grid_zeroed[is_gap] <- 0

  # Rebuild mu_fun with zeroed grid
  mu_fun_zeroed <- approxfun(grid, mu_grid_zeroed, rule = 2)

  # Recompute integral_bg via trapezoidal rule on the zeroed grid
  dx <- diff(grid)
  area <- dx * (head(mu_grid_zeroed, -1) + tail(mu_grid_zeroed, -1)) / 2
  integral_bg_zeroed <- sum(area)

  # Recompute mu_vec at event times
  mu_vec_zeroed <- mu_fun_zeroed(inhom_bg$times)
  mu_vec_zeroed <- pmax(mu_vec_zeroed, 1e-12)

  # Update mu_fit in-place
  mu_fit_zeroed <- mu_fit
  mu_fit_zeroed$mu_grid <- mu_grid_zeroed
  mu_fit_zeroed$mu_fun  <- mu_fun_zeroed

  list(
    mu_vec      = mu_vec_zeroed,
    integral_bg = integral_bg_zeroed,
    times       = inhom_bg$times,
    mu_fit      = mu_fit_zeroed,
    Lambda_fun  = inhom_bg$Lambda_fun,  # note: Lambda_fun is not recomputed; use integral_bg_zeroed
    gap_intervals = gap_intervals
  )
}

#' Subset edges to the first session (before the first gap > gap_threshold).
subset_first_session <- function(df, gap_threshold = 1.0) {
  df <- df[order(df$time), , drop = FALSE]
  gaps <- diff(df$time)
  first_gap <- which(gaps > gap_threshold)[1]
  if (is.na(first_gap)) return(df)  # no gaps — all data is one session
  df[seq_len(first_gap), , drop = FALSE]
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

  impossible_any <- edges_df$time < pmax(node_entry_time[edges_df$tail], node_entry_time[edges_df$head])
  u <- pmin(edges_df$tail, edges_df$head)
  v <- pmax(edges_df$tail, edges_df$head)
  edge_key <- u + n_all * v
  impossible_any <- impossible_any | duplicated(edge_key)
  n_impossible_any <- sum(impossible_any, na.rm = TRUE)

  last_activity <- node_entry_time
  n_edges <- nrow(edges_df)
  n_impossible <- integer(length(trunc_grid))
  n_candidate_edges <- numeric(length(trunc_grid))

  entry_sorted <- node_entry_time
  if (!isTRUE(all(diff(entry_sorted) >= -1e-12))) {
    entry_sorted <- sort(entry_sorted)
  }
  present_count_at_edge <- integer(n_edges)
  ptr <- 0L
  for (i in seq_len(n_edges)) {
    t_i <- edges_df$time[i]
    while (ptr < n_all && entry_sorted[ptr + 1L] <= t_i) ptr <- ptr + 1L
    present_count_at_edge[i] <- ptr
  }

  get_active_nodes <- function(n_present, truncation) {
    if (n_present <= 1L) return(integer(0))
    if (mark_decay == "node_entrance") {
      window_start <- max(1L, n_present - truncation + 1L)
      return(window_start:n_present)
    }
    ids <- seq_len(n_present)
    ord <- order(last_activity[ids], ids, decreasing = TRUE)
    sort(ord[seq_len(min(truncation, n_present))])
  }

  for (k in seq_along(trunc_grid)) {
    tr <- as.integer(trunc_grid[k])
    if (!is.finite(tr) || tr < 2L) tr <- 2L
    last_activity <- node_entry_time
    imp <- 0L
    cand_total <- 0
    for (i in seq_len(n_edges)) {
      if (impossible_any[i]) { imp <- imp + 1L; next }
      n_present <- present_count_at_edge[i]
      active <- get_active_nodes(n_present, tr)
      cand_total <- cand_total + (length(active) * (length(active) - 1L) / 2)
      tail_i <- edges_df$tail[i]
      head_i <- edges_df$head[i]
      if (!(tail_i %in% active && head_i %in% active)) imp <- imp + 1L
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

#' Run a fit + optional sim test + optional GOF.
#' Returns a list with fit, sim_test, gof components.
run_fit_block <- function(net, inhom_bg, time_window, label,
                          params_init, p_scale,
                          formula_rhs, truncation, mark_decay, growth_only,
                          max_iter, n_cores,
                          run_sim_test = FALSE, n_sim_test = 3L,
                          run_gof = FALSE, n_gof = 100L, n_cores_gof = n_cores,
                          n_gof_outer = 0L, seed_events_gof = 20L) {
  cat("\n=== ", label, " ===\n")
  cat("  Network:", network.edgecount(net), "edges |", network.size(net), "nodes\n")
  cat("  Time window:", sprintf("[%.3f, %.3f]", time_window[1], time_window[2]), "\n\n")

  # Fit
  cat("--- Fitting hawkesNet (", label, ") ---\n")
  t_fit <- proc.time()
  fit <- fit_hawkesNet(
    params_init = params_init,
    time_window = time_window,
    mark_filtration = net,
    PMF_mark = PMF_mark_CS,
    mu_vec = inhom_bg$mu_vec,
    integral_bg = inhom_bg$integral_bg,
    formula_RHS = formula_rhs,
    truncation = truncation,
    mark_decay = mark_decay,
    growth_only = growth_only,
    max_node_time = max(get_times(net)$node_times),
    method = "Nelder-Mead",
    maxit = max_iter,
    reltol = 1e-8,
    trace = 1,  # increased trace for more progress info
    verbose = TRUE,
    fixed_params = if (is.null(inhom_bg$mu_vec)) "K" else c("mu", "K"),
    parscale = p_scale,
    cache_intensity = TRUE,
    combine_intensity = TRUE,
    cores = n_cores
  )
  cat("  Fit completed in", round((proc.time() - t_fit)[3], 2), "s\n")
  cat("  Final log-likelihood:", if (!is.null(fit$fit$value)) -fit$fit$value else "failed", "\n")
  cat("  Convergence status:", fit$fit$convergence, "(0 = success)\n")

  if (!is.null(fit$fit_table)) {
    cat("\n--- Fit results (", label, ") ---\n")
    print(fit$fit_table, max = NULL)
  }

  # Optional sim test
  sim_test_res <- NULL
  if (run_sim_test && !is.null(fit$params)) {
    cat("\n--- Simulation test (", label, ") ---\n")
    pfit <- fit$params
    seed_net_test <- NULL
    seed_times_test <- NULL
    if (seed_events_gof > 0) {
      all_times <- get_times(net)$times
      if (length(all_times) >= seed_events_gof) {
        t_seed <- all_times[seed_events_gof]
        seed_net_test <- filtration_to_net(net, t_seed, equals = TRUE)
        seed_times_test <- all_times[1:seed_events_gof]
        cat("  Seeding with first", seed_events_gof, "events (up to t =", round(t_seed, 4), ")\n")
      }
    }
    sims_test <- vector("list", n_sim_test)
    for (i in seq_len(n_sim_test)) {
      sims_test[[i]] <- tryCatch({
        sim_hawkesNet(
          params = pfit,
          time_window = time_window,
          PMF_mark = PMF_mark_CS,
          cond_intensity = cond_intensity,
          formula_RHS = formula_rhs,
          truncation = truncation,
          mark_decay = mark_decay,
          growth_only = growth_only,
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
    cat("  Simulated:", n_sim_ok, "/", n_sim_test, "succeeded\n")
    if (n_sim_ok > 0) {
      cat("    n_edges:  ", paste(na.omit(n_edges_sim), collapse = ", "), "\n")
      cat("    n_nodes:  ", paste(na.omit(n_nodes_sim), collapse = ", "), "\n")
      cat("    mean_deg: ", paste(round(na.omit(mean_deg_sim), 2), collapse = ", "), "\n")
    }
    sim_test_res <- list(sims = sims_test, n_obs = n_obs, e_obs = e_obs, mean_deg_obs = mean_deg_obs)
  }

  # Optional GOF
  gof_res <- NULL
  if (run_gof) {
    cat("\n--- GOF (", label, ") ---\n")
    cat("  Running", n_gof, "simulations on", n_cores_gof, "cores...\n")
    t_gof <- proc.time()
    gof_fun <- if (exists("gof", mode = "function")) get("gof", mode = "function") else hawkesNet::gof
    gof_formals <- names(formals(gof_fun))
    gof_args <- list(
      fit = fit,
      net_obs = net,
      params_init = params_init,
      PMF_mark = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS = formula_rhs,
      time_window = time_window,
      truncation = as.integer(truncation),
      mark_decay = as.character(mark_decay),
      growth_only = isTRUE(growth_only),
      max_node_time = max(get_times(net)$node_times),
      inhom_bg = inhom_bg,
      n_sim = as.integer(n_gof),
      cores = as.integer(n_cores_gof),
      cores_outer = if (n_gof_outer > 0L) as.integer(n_gof_outer) else NULL,
      max_deg = 30L,
      k_esp = 30L,
      degree = 0L,
      esp = 0L,
      mu_multiplier = 5,
      seed_events = as.integer(seed_events_gof),
      verbose = TRUE
    )
    gof_args <- gof_args[names(gof_args) %in% gof_formals]
    gof_res <- tryCatch({
      do.call(gof_fun, gof_args)
    }, error = function(e) {
      cat("  GOF failed for ", label, ": ", conditionMessage(e), "\n")
      NULL
    })
    cat("  GOF completed in", round((proc.time() - t_gof)[3], 2), "s\n")
    if (!is.null(gof_res$plots) && requireNamespace("ggplot2", quietly = TRUE)) {
      if (!is.null(gof_res$plots$degree_plot)) print(gof_res$plots$degree_plot)
      if (!is.null(gof_res$plots$esp_plot)) print(gof_res$plots$esp_plot)
      if (!is.null(gof_res$plots$geodist_plot)) print(gof_res$plots$geodist_plot)
      if (!is.null(gof_res$plots$waiting_times_plot)) print(gof_res$plots$waiting_times_plot)
    }
  }

  list(fit = fit, sim_test = sim_test_res, gof = gof_res)
}

# =============================================================================
# Load data
# =============================================================================
raw <- read.table(system.file("extdata", "ht09_contact_list.dat", package = "hawkesNet"))
df <- data.frame(
  time = raw$V1 / 20 / 3600,  # convert to hours immediately
  from = raw$V2,
  to   = raw$V3
)

set.seed(1)
df$time <- df$time + rnorm(nrow(df), 0, 0.01 / 3600) # jitter in hours (0.01s)
df$time <- df$time - min(df$time) # shift to 0
df <- df[order(df$time), ] # ensure sorted

# Build the FULL network (all days, no compression, no subsetting)
obj_full <- make_hypertext_net(df, use_first_contact_only = USE_FIRST_CONTACT_ONLY, max_edges = MAX_EDGES)

# Identify gap intervals in the raw edge times (before first-contact collapsing)
# Normalize raw times to [0, 1] based on full window for gap identification
all_edge_times_sorted <- sort(obj_full$edges$time)
full_range <- range(all_edge_times_sorted)
norm_times <- (all_edge_times_sorted - full_range[1]) / (full_range[2] - full_range[1])
gap_intervals <- find_gap_intervals(norm_times, gap_threshold = GAP_THRESHOLD / (full_range[2] - full_range[1]))
cat("  Gap intervals found (normalized):", nrow(gap_intervals), "\n")
if (nrow(gap_intervals) > 0) {
  for (g in seq_len(nrow(gap_intervals))) {
    cat(sprintf("    Gap %d: [%.4f, %.4f] (duration %.1f h)\n",
                g, gap_intervals$start[g], gap_intervals$end[g],
                (gap_intervals$end[g] - gap_intervals$start[g]) * (full_range[2] - full_range[1])))
  }
}

# =============================================================================
# Fit 1: First session only (before first overnight gap)
# =============================================================================
cat("\n######################################################################\n")
cat("## FIT 1: First session only\n")
cat("######################################################################\n")

# Subset edges to first session from the raw data, then rebuild network
df_day1 <- subset_first_session(df, gap_threshold = GAP_THRESHOLD)
obj_day1 <- make_hypertext_net(df_day1, use_first_contact_only = USE_FIRST_CONTACT_ONLY, max_edges = MAX_EDGES)
net_day1 <- normalize_times_01(obj_day1$net) # Normalize to [0, 1]

times_day1 <- get_times(net_day1)$times
time_window_day1 <- c(min(times_day1), max(times_day1))
cat("  Day-1 network:", length(times_day1), "events |", network.size(net_day1), "nodes |",
    network.edgecount(net_day1), "edges\n")
cat("  Time window:", sprintf("[%.3f, %.3f]", time_window_day1[1], time_window_day1[2]), "\n")

# Truncation for day-1
n_nodes_day1 <- network.size(net_day1)
TRUNCATION_DAY1 <- if (is.na(TRUNCATION)) n_nodes_day1 else TRUNCATION
cat("  TRUNCATION:", TRUNCATION_DAY1, "\n")
verify_truncation(obj_day1$edges, node_entry_time = net_day1 %v% "time",
                  truncation = TRUNCATION_DAY1, mark_decay = MARK_DECAY)

# KDE background for day-1
cat("--- Using homogeneous background (day-1) ---\n")
inhom_bg_day1 <- list(
  mu_vec = NULL,
  integral_bg = NULL,
  times = times_day1,
  mu_fit = NULL,
  Lambda_fun = NULL
)

# Model init
exp_cs <- expected_params_PMF_mark_CS(net_day1, FORMULA_RHS)
n_cs <- if (!is.na(exp_cs$CS_params_length)) exp_cs$CS_params_length else 3L
cat("  CS_params length:", n_cs, "\n")

mu_init_day1 <- length(times_day1) / (time_window_day1[2] - time_window_day1[1])
params_init_day1 <- list(
  mu = mu_init_day1,
  beta_overall = 0.3,
  K = 0.5,
  beta_edges = 0.3,
  node_lambda = 0.1,
  CS_params = c(-8, -5, rep(0, n_cs - 2))
)

  p_scale_day1 <- c(
    mu = 1,
    beta_overall = 0.1,
  beta_edges = 0.1,
  node_lambda = 0.5,
  setNames(rep(0.1, n_cs), paste0("CS_params", seq_len(n_cs)))
)

res_day1 <- run_fit_block(
  net = net_day1, inhom_bg = inhom_bg_day1, time_window = time_window_day1,
  label = "Fit 1: First session",
  params_init = params_init_day1, p_scale = p_scale_day1,
  formula_rhs = FORMULA_RHS, truncation = TRUNCATION_DAY1,
  mark_decay = MARK_DECAY, growth_only = GROWTH_ONLY,
  max_iter = MAX_ITER, n_cores = N_CORES,
  run_sim_test = RUN_SIM_TEST, n_sim_test = N_SIM_TEST,
  run_gof = RUN_GOF_DAY1, n_gof = N_GOF, n_cores_gof = N_CORES_GOF,
  n_gof_outer = N_GOF_OUTER, seed_events_gof = SEED_EVENTS_GOF
)

# =============================================================================
# Fit 2: Full data with mu=0 during gaps (optional)
# =============================================================================
res_full <- NULL
inhom_bg_full <- NULL
net_full <- NULL
time_window_full <- NULL
TRUNCATION_FULL <- NULL

if (RUN_FULL_FIT) {
  cat("\n######################################################################\n")
  cat("## FIT 2: Full data (mu=0 in overnight gaps)\n")
  cat("######################################################################\n")

  net_full <- normalize_times_01(obj_full$net) # Normalize to [0, 1]
  times_full <- get_times(net_full)$times
  time_window_full <- c(min(times_full), max(times_full))
  cat("  Full network:", length(times_full), "events |", network.size(net_full), "nodes |",
      network.edgecount(net_full), "edges\n")
  cat("  Time window:", sprintf("[%.3f, %.3f]", time_window_full[1], time_window_full[2]), "\n")

  # Truncation for full data
  n_nodes_full <- network.size(net_full)
  TRUNCATION_FULL <- if (is.na(TRUNCATION)) n_nodes_full else TRUNCATION
  cat("  TRUNCATION:", TRUNCATION_FULL, "\n")
  verify_truncation(obj_full$edges, node_entry_time = net_full %v% "time",
                    truncation = TRUNCATION_FULL, mark_decay = MARK_DECAY)

  # KDE background for full data, then zero out gaps
  cat("--- Estimating KDE background (full, then zeroing gaps) ---\n")
  t_bg <- proc.time()
  inhom_bg_full_raw <- prepare_inhomogeneous_background(
    net_full, time_attr = "time",
    bw = if (is.finite(BW)) BW else NULL,
    grid_n = GRID_N
  )
  inhom_bg_full <- zero_gaps_inhom_bg(inhom_bg_full_raw, gap_intervals)
  cat("  KDE + gap-zeroing done in", round((proc.time() - t_bg)[3], 2), "s\n")
  cat("  integral_bg (raw):", round(inhom_bg_full_raw$integral_bg, 2),
      "| integral_bg (gap-zeroed):", round(inhom_bg_full$integral_bg, 2), "\n")

  # Model init (use same formula; may want different starting values)
  exp_cs_full <- expected_params_PMF_mark_CS(net_full, FORMULA_RHS)
  n_cs_full <- if (!is.na(exp_cs_full$CS_params_length)) exp_cs_full$CS_params_length else 3L

  mu_init_full <- inhom_bg_full$integral_bg / (time_window_full[2] - time_window_full[1])
  params_init_full <- list(
    mu = mu_init_full,
    beta_overall = 0.3,
    K = 0.5,
    beta_edges = 0.3,
    node_lambda = 0.1,
    CS_params = c(-8, -5, rep(0, n_cs_full - 2))
  )

  p_scale_full <- c(
    beta_overall = 0.1,
    beta_edges = 0.1,
    node_lambda = 0.5,
    setNames(rep(0.1, n_cs_full), paste0("CS_params", seq_len(n_cs_full)))
  )

  res_full <- run_fit_block(
    net = net_full, inhom_bg = inhom_bg_full, time_window = time_window_full,
    label = "Fit 2: Full data (mu=0 in gaps)",
    params_init = params_init_full, p_scale = p_scale_full,
    formula_rhs = FORMULA_RHS, truncation = TRUNCATION_FULL,
    mark_decay = MARK_DECAY, growth_only = GROWTH_ONLY,
    max_iter = MAX_ITER, n_cores = N_CORES,
    run_sim_test = RUN_SIM_TEST, n_sim_test = N_SIM_TEST,
    run_gof = RUN_GOF_FULL, n_gof = N_GOF, n_cores_gof = N_CORES_GOF,
    n_gof_outer = N_GOF_OUTER, seed_events_gof = SEED_EVENTS_GOF
  )
}

# =============================================================================
# Save results
# =============================================================================
save_list <- list(
  # Fit 1: first session
  fit_day1      = res_day1$fit,
  gof_day1      = res_day1$gof,
  sim_test_day1 = res_day1$sim_test,
  net_day1      = net_day1,
  edges_day1    = obj_day1$edges,
  inhom_bg_day1 = inhom_bg_day1,
  params_init_day1 = params_init_day1,
  time_window_day1 = time_window_day1,
  truncation_day1  = TRUNCATION_DAY1,
  # Fit 2: full data (mu=0 in gaps)
  fit_full      = if (!is.null(res_full)) res_full$fit else NULL,
  gof_full      = if (!is.null(res_full)) res_full$gof else NULL,
  sim_test_full = if (!is.null(res_full)) res_full$sim_test else NULL,
  net_full      = net_full,
  edges_full    = obj_full$edges,
  inhom_bg_full = inhom_bg_full,
  params_init_full = if (exists("params_init_full")) params_init_full else NULL,
  time_window_full = time_window_full,
  truncation_full  = TRUNCATION_FULL,
  # Shared metadata
  gap_intervals = gap_intervals,
  gap_threshold = GAP_THRESHOLD,
  formula_rhs   = FORMULA_RHS,
  N_GOF         = N_GOF,
  SEED_EVENTS_GOF = SEED_EVENTS_GOF,
  RUN_FULL_FIT  = RUN_FULL_FIT
)

cluster_output_dir <- file.path(PKG_ROOT, "cluster_output")
rds_path_primary <- file.path(cluster_output_dir, "results_hypertext_full.RDS")

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
cat("\n######################################################################\n")
cat("## Hypertext study finished at", as.character(Sys.time()), "\n")
cat("######################################################################\n")
cat("Done.\n")
