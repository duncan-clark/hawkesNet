## =============================================================================
## Hypertext 2009 conference: consolidated hawkesNet fits
## =============================================================================
## Run from package root:
##   Rscript inst/hypertext_conference/hypertext_fits.R
##
## Quick local mode:
##   LOCAL_QUICK=TRUE Rscript inst/hypertext_conference/hypertext_fits.R
##
## Cluster mode (SLURM):
##   sbatch inst/hypertext_conference/run_hypertext.slurm
##
## Publication fit (simple data only, Formula A = tri+star,
## mark_decay="activity", GROWTH_ONLY=FALSE):
##
##   Fit B — K fixed (default 1), edges free
##
## The fit fixes node_lambda at the observed nodes/events ratio.
## GOF is run for the fit (controllable via env vars).
## =============================================================================

library(hawkesNet)
library(dplyr)
library(network)
library(ernm)
library(sna)

# =============================================================================
# Config
# =============================================================================
PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) {
  Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  getwd()
}
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) PKG_ROOT <- getwd()

OUTPUT_DIR <- file.path(PKG_ROOT, "cluster_output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

ON_SLURM <- nzchar(Sys.getenv("SLURM_JOB_ID")) || nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))
LOCAL_QUICK <- isTRUE(as.logical(Sys.getenv("LOCAL_QUICK", if (ON_SLURM) "FALSE" else "TRUE")))

N_CORES  <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 25L)))
MAX_ITER <- as.integer(Sys.getenv("MAX_ITER", if (LOCAL_QUICK) 200L else 5000L))

N_GOF       <- max(1L, as.integer(Sys.getenv("N_GOF", N_CORES)))
# Optional PSOCK parallelism for GOF. When set, gof() ignores `cores` and uses
# `cores_outer` workers (capped at 60 inside gof()).
N_GOF_OUTER <- as.integer(Sys.getenv("GOF_CORES_OUTER", N_CORES))
SEED_EVENTS_GOF <- as.integer(Sys.getenv("SEED_EVENTS_GOF", 20L))

RUN_FIT_A <- isTRUE(as.logical(Sys.getenv("RUN_FIT_A", "TRUE")))
RUN_FIT_B <- isTRUE(as.logical(Sys.getenv("RUN_FIT_B", "TRUE")))
RUN_FIT_C <- isTRUE(as.logical(Sys.getenv("RUN_FIT_C", "TRUE")))
RUN_FIT_D <- isTRUE(as.logical(Sys.getenv("RUN_FIT_D", "TRUE")))

RUN_GOF_A <- isTRUE(as.logical(Sys.getenv("RUN_GOF_A", "TRUE")))
RUN_GOF_B <- isTRUE(as.logical(Sys.getenv("RUN_GOF_B", "TRUE")))
RUN_GOF_C <- isTRUE(as.logical(Sys.getenv("RUN_GOF_C", "TRUE")))
RUN_GOF_D <- isTRUE(as.logical(Sys.getenv("RUN_GOF_D", "TRUE")))

# Fixed model settings
MARK_DECAY  <- "activity"
GROWTH_ONLY <- FALSE
EDGES_INIT  <- -5
K_FIXED     <- as.numeric(Sys.getenv("K_FIXED", 0.5))
GAP_THRESHOLD <- 1.0
USE_FIRST_CONTACT_ONLY <- TRUE

FORMULA_A <- "edges + triangles + star(c(2,3))"
FORMULA_B <- "edges + gwesp(0.5) + gwdegree(0.5)"
FORMULA_C <- "edges + degree(2:3) + esp(1:2)"
FORMULA_D <- "edges + triangles + star(c(2,3)) + gwesp(0.5) + gwdegree(0.5)"

cat("=== Hypertext Conference Fits (simple data only) ===\n")
cat("  Mode:", if (ON_SLURM) "SLURM" else if (LOCAL_QUICK) "Local (quick)" else "Local", "\n")
cat("  N_CORES:", N_CORES, "| MAX_ITER:", MAX_ITER, "| N_GOF:", N_GOF, "\n")
fits_abcd <- c(A = RUN_FIT_A, B = RUN_FIT_B, C = RUN_FIT_C, D = RUN_FIT_D)
cat("  Fits to run:", paste(names(fits_abcd)[fits_abcd], collapse = ", "), "\n")
gofs_abcd <- c(A = RUN_GOF_A, B = RUN_GOF_B, C = RUN_GOF_C, D = RUN_GOF_D)
cat("  GOF to run:", paste(names(gofs_abcd)[gofs_abcd], collapse = ", "), "\n")
cat("  Formula:", FORMULA_A, "\n")
cat("  mark_decay:", MARK_DECAY, "\n")
cat("  Fit B: K=", K_FIXED, " fixed, edges free\n", sep = "")
cat("  node_lambda fixed at observed nodes/events\n\n")

# =============================================================================
# Helpers
# =============================================================================

make_hypertext_net <- function(df, use_first_contact_only = TRUE) {
  stopifnot(all(c("from", "to", "time") %in% names(df)))
  swap <- df$from > df$to
  df[swap, c("from", "to")] <- df[swap, c("to", "from")]
  df <- df %>% distinct()
  df$time <- as.numeric(df$time)
  df <- df[is.finite(df$time), ]
  df <- df[order(df$time), , drop = FALSE]

  if (use_first_contact_only) {
    df <- df %>%
      dplyr::group_by(.data$from, .data$to) %>%
      dplyr::summarise(time = min(.data$time), .groups = "drop") %>%
      arrange(time)
  }

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

  el <- as.matrix(df[, c("tail", "head")])
  net <- network::network(el, matrix.type = "edgelist", directed = FALSE)
  network::set.edge.attribute(net, "time", df$time)
  network::set.vertex.attribute(net, "time", as.numeric(entry_time_by_node))
  network::delete.vertex.attribute(net, "vertex.names")

  list(net = net, edges = df, nodes = nodes, id_map = id_map)
}

subset_first_session <- function(df, gap_threshold = 1.0) {
  df <- df[order(df$time), , drop = FALSE]
  gaps <- diff(df$time)
  first_gap <- which(gaps > gap_threshold)[1]
  if (is.na(first_gap)) return(df)
  df[seq_len(first_gap), , drop = FALSE]
}

run_single_fit <- function(net, time_window, label,
                           formula_rhs, mark_decay, growth_only,
                           max_iter, n_cores,
                           fixed_params, params_override = NULL,
                           run_gof, n_gof, n_gof_outer, seed_events_gof) {
  cat("\n######################################################################\n")
  cat("## ", label, "\n")
  cat("######################################################################\n")

  times <- get_times(net)$times
  n_events <- length(times)
  n_edges  <- network.edgecount(net)
  n_nodes  <- network.size(net)
  trunc    <- n_nodes

  cat("  Network:", n_events, "events |", n_nodes, "nodes |", n_edges, "edges\n")
  cat("  Time window:", sprintf("[%.4f, %.4f]", time_window[1], time_window[2]), "\n")
  cat("  Edges/event:", round(n_edges / n_events, 2), "\n")
  cat("  Nodes/event:", round(n_nodes / n_events, 4), "(fixed as node_lambda)\n")
  cat("  Formula:", formula_rhs, "\n")

  exp_cs <- expected_params_PMF_mark_CS(net, formula_rhs)
  n_cs <- if (!is.na(exp_cs$CS_params_length)) exp_cs$CS_params_length else 4L
  cat("  CS_params length:", n_cs, "\n")

  mu_init <- n_events / (time_window[2] - time_window[1])
  m_init  <- n_edges / n_events
  node_lambda_init <- n_nodes / n_events

  cs_init <- rep(0, n_cs)
  cs_init[1] <- EDGES_INIT

  params_init <- list(
    mu = mu_init,
    beta_overall = 0.3,
    K = K_FIXED,
    beta_edges = 0.3,
    node_lambda = node_lambda_init,
    m = m_init,
    CS_params = cs_init
  )
  if (!is.null(params_override)) {
    for (nm in names(params_override)) params_init[[nm]] <- params_override[[nm]]
  }

  p_scale <- c(
    mu = 1, beta_overall = 0.1, K = 0.1, beta_edges = 0.1,
    m = 0.5,
    setNames(rep(0.1, n_cs), paste0("CS_params", seq_len(n_cs)))
  )
  p_scale <- p_scale[!names(p_scale) %in% fixed_params]

  cat("  Fixed params:", paste(fixed_params, collapse = ", "), "\n")
  cat("  params_init:\n")
  for (nm in names(params_init)) {
    val <- params_init[[nm]]
    if (length(val) > 1) {
      cat(sprintf("    %s = [%s]\n", nm, paste(round(val, 4), collapse = ", ")))
    } else {
      cat(sprintf("    %s = %.4f\n", nm, val))
    }
  }

  # Fit
  cat("\n--- Fitting ---\n")
  t_fit <- proc.time()
  fit <- tryCatch({
    fit_hawkesNet(
      params_init = params_init,
      time_window = time_window,
      mark_filtration = net,
      PMF_mark = PMF_mark_CS,
      formula_RHS = formula_rhs,
      truncation = trunc,
      mark_decay = mark_decay,
      growth_only = growth_only,
      max_node_time = max(get_times(net)$node_times),
      method = "Nelder-Mead",
      maxit = max_iter,
      trace = 0,
      verbose = FALSE,
      fixed_params = fixed_params,
      parscale = p_scale,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      cores = n_cores
    )
  }, error = function(e) {
    cat(sprintf("  FIT FAILED: %s\n", e$message))
    NULL
  })
  t_fit_elapsed <- (proc.time() - t_fit)[3]
  cat(sprintf("  Fit completed in %.1f s\n", t_fit_elapsed))

  if (!is.null(fit) && !is.null(fit$fit_table)) {
    cat("\n--- Fit results ---\n")
    print(fit$fit_table, max = NULL)
  }

  # GOF
  gof_res <- NULL
  t_gof_elapsed <- 0
  if (run_gof && !is.null(fit) && !is.null(fit$params)) {
    cat(sprintf("\n--- GOF (%d sims) ---\n", n_gof))
    t_gof <- proc.time()
    gof_fun <- if (existsFunction("gof")) get("gof") else if (existsFunction("gof_hawkesNet")) get("gof_hawkesNet") else NULL
    if (!is.null(gof_fun)) {
      gof_formals <- names(formals(gof_fun))
      gof_args <- list(
        fit = fit,
        net_obs = net,
        params_init = params_init,
        PMF_mark = PMF_mark_CS,
        cond_intensity = cond_intensity,
        formula_RHS = formula_rhs,
        time_window = time_window,
        truncation = as.integer(trunc),
        mark_decay = mark_decay,
        growth_only = growth_only,
        max_node_time = max(get_times(net)$node_times),
        n_sim = as.integer(n_gof),
        cores = as.integer(n_cores),
        cores_outer = if (n_gof_outer > 0L) as.integer(n_gof_outer) else NULL,
        max_deg = 30L,
        k_esp = 30L,
        degree = 0L,
        esp = 0L,
        mu_multiplier = 5,
        seed_events = as.integer(seed_events_gof),
        verbose = FALSE
      )
      gof_args <- gof_args[names(gof_args) %in% gof_formals]
      gof_res <- tryCatch({
        do.call(gof_fun, gof_args)
      }, error = function(e) {
        cat(sprintf("  GOF FAILED: %s\n", e$message))
        NULL
      })
    } else {
      cat("  GOF function not found; skipping.\n")
    }
    t_gof_elapsed <- (proc.time() - t_gof)[3]
    cat(sprintf("  GOF completed in %.1f s\n", t_gof_elapsed))
  }

  t_total <- t_fit_elapsed + t_gof_elapsed
  cat(sprintf("\n  === %s: fit %.1f s + GOF %.1f s = total %.1f s ===\n",
              label, t_fit_elapsed, t_gof_elapsed, t_total))

  # Free the intensity closure cache to avoid memory bloat across sequential fits
  if (!is.null(fit)) {
    fit$intens_funcs <- NULL
  }
  gc()

  list(
    fit = fit,
    gof = gof_res,
    label = label,
    formula_rhs = formula_rhs,
    params_init = params_init,
    net = net,
    time_window = time_window,
    truncation = trunc,
    time_fit = t_fit_elapsed,
    time_gof = t_gof_elapsed
  )
}

# =============================================================================
# Load data
# =============================================================================
t_wall_start <- proc.time()
cat("\n=== Loading Data ===\n")
raw <- read.table(system.file("extdata", "ht09_contact_list.dat", package = "hawkesNet"))

# --- Simple data (jittered so each event = 1 edge) ---
df_simple <- data.frame(
  time = raw$V1 / 3600,
  from = raw$V2,
  to   = raw$V3
)
set.seed(1)
df_simple$time <- df_simple$time + rnorm(nrow(df_simple), 0, 0.01 / 3600)
df_simple$time <- df_simple$time - min(df_simple$time)
df_simple <- df_simple[order(df_simple$time), ]

df_simple_day1 <- subset_first_session(df_simple, gap_threshold = GAP_THRESHOLD)
obj_simple <- make_hypertext_net(df_simple_day1, use_first_contact_only = USE_FIRST_CONTACT_ONLY)
net_simple <- normalize_times_01(obj_simple$net)
times_simple <- get_times(net_simple)$times
tw_simple <- c(min(times_simple), max(times_simple))

cat(sprintf("  Simple day-1: %d events | %d nodes | %d edges\n",
            length(times_simple), network.size(net_simple), network.edgecount(net_simple)))


# =============================================================================
# Run fits
# =============================================================================

results <- list()

save_incremental <- function() {
  rds_tmp <- file.path(OUTPUT_DIR, "results_hypertext_partial.RDS")
  tryCatch({
    saveRDS(list(
      results = results,
      net_simple = net_simple,
      edges_simple = obj_simple$edges,
      tw_simple = tw_simple,
      formula = FORMULA_A,
      mark_decay = MARK_DECAY, edges_init = EDGES_INIT, K_fixed = K_FIXED,
      N_GOF = N_GOF, MAX_ITER = MAX_ITER, SEED_EVENTS_GOF = SEED_EVENTS_GOF
    ), rds_tmp)
    cat(sprintf("  [checkpoint] Saved partial results: %s\n", rds_tmp))
  }, error = function(e) cat("  [checkpoint] Save failed:", e$message, "\n"))
}

# --- Fit A: K fixed, edges free ---
if (RUN_FIT_A) {
  results$fitA <- run_single_fit(
    net = net_simple,
    time_window = tw_simple,
    label = "Fit A: K fixed, edges free",
    formula_rhs = FORMULA_A,
    mark_decay = MARK_DECAY,
    growth_only = GROWTH_ONLY,
    max_iter = MAX_ITER,
    n_cores = N_CORES,
    fixed_params = c("K", "node_lambda"),
    params_override = list(K = K_FIXED),
    run_gof = RUN_GOF_A,
    n_gof = N_GOF,
    n_gof_outer = N_GOF_OUTER,
    seed_events_gof = SEED_EVENTS_GOF
  )
  save_incremental()
}

# --- Fit B: K fixed, edges free ---
if (RUN_FIT_B) {
  results$fitB <- run_single_fit(
    net = net_simple,
    time_window = tw_simple,
    label = "Fit B: K fixed, edges free",
    formula_rhs = FORMULA_B,
    mark_decay = MARK_DECAY,
    growth_only = GROWTH_ONLY,
    max_iter = MAX_ITER,
    n_cores = N_CORES,
    fixed_params = c("K", "node_lambda"),
    params_override = list(K = K_FIXED),
    run_gof = RUN_GOF_B,
    n_gof = N_GOF,
    n_gof_outer = N_GOF_OUTER,
    seed_events_gof = SEED_EVENTS_GOF
  )
  save_incremental()
}

# --- Fit C: K fixed, edges free ---
if (RUN_FIT_C) {
  results$fitC <- run_single_fit(
    net = net_simple,
    time_window = tw_simple,
    label = "Fit C: K fixed, edges free",
    formula_rhs = FORMULA_C,
    mark_decay = MARK_DECAY,
    growth_only = GROWTH_ONLY,
    max_iter = MAX_ITER,
    n_cores = N_CORES,
    fixed_params = c("K", "node_lambda"),
    params_override = list(K = K_FIXED),
    run_gof = RUN_GOF_C,
    n_gof = N_GOF,
    n_gof_outer = N_GOF_OUTER,
    seed_events_gof = SEED_EVENTS_GOF
  )
  save_incremental()
}

# --- Fit D: K fixed, edges free ---
if (RUN_FIT_D) {
  results$fitD <- run_single_fit(
    net = net_simple,
    time_window = tw_simple,
    label = "Fit D: K fixed, edges free",
    formula_rhs = FORMULA_D,
    mark_decay = MARK_DECAY,
    growth_only = GROWTH_ONLY,
    max_iter = MAX_ITER,
    n_cores = N_CORES,
    fixed_params = c("K", "node_lambda"),
    params_override = list(K = K_FIXED),
    run_gof = RUN_GOF_D,
    n_gof = N_GOF,
    n_gof_outer = N_GOF_OUTER,
    seed_events_gof = SEED_EVENTS_GOF
  )
  save_incremental()
}

# results$ernm_fit <- ernm(net_simple ~ edges + gwesp(0.5) + gwdegree(0.5))

# =============================================================================
# Summary table
# =============================================================================
cat("\n######################################################################\n")
cat("## SUMMARY\n")
cat("######################################################################\n")

for (nm in names(results)) {
  res <- results[[nm]]
  if (is.null(res)) next
  cat(sprintf("\n--- %s ---\n", res$label))
  cat(sprintf("  Fit time: %.1f s | GOF time: %.1f s\n", res$time_fit, res$time_gof))
  if (!is.null(res$fit) && !is.null(res$fit$fit_table)) {
    print(res$fit$fit_table)
  }
}

# =============================================================================
# Save results
# =============================================================================
save_list <- list(
  results = results,
  net_simple = net_simple,
  edges_simple = obj_simple$edges,
  tw_simple = tw_simple,
  formula = FORMULA_A,
  mark_decay = MARK_DECAY,
  edges_init = EDGES_INIT,
  K_fixed = K_FIXED,
  N_GOF = N_GOF,
  MAX_ITER = MAX_ITER,
  SEED_EVENTS_GOF = SEED_EVENTS_GOF
)

rds_path <- file.path(OUTPUT_DIR, "results_hypertext.RDS")
tryCatch({
  saveRDS(save_list, rds_path)
  cat(sprintf("\nResults saved: %s\n", rds_path))
}, error = function(e) {
  cat(sprintf("  Save failed: %s\n", e$message))
  fallback <- file.path(PKG_ROOT, "inst", "hypertext_conference", "results_hypertext.rds")
  tryCatch({
    saveRDS(save_list, fallback)
    cat(sprintf("  Fallback save: %s\n", fallback))
  }, error = function(e2) cat("  Fallback also failed:", e2$message, "\n"))
})

# GOF plots
if (requireNamespace("ggplot2", quietly = TRUE)) {
  for (nm in names(results)) {
    res <- results[[nm]]
    if (is.null(res) || is.null(res$gof) || is.null(res$gof$plots)) next
    safe_label <- gsub("[^a-zA-Z0-9_]", "_", tolower(res$label))
    for (plot_name in names(res$gof$plots)) {
      p <- res$gof$plots[[plot_name]]
      if (!is.null(p)) {
        fname <- file.path(OUTPUT_DIR, sprintf("gof_%s_%s.pdf", safe_label, plot_name))
        tryCatch(
          ggplot2::ggsave(filename = fname, plot = p, width = 10, height = 8),
          error = function(e) cat("  Plot save failed:", fname, e$message, "\n")
        )
      }
    }
  }
}

t_wall_elapsed <- (proc.time() - t_wall_start)[3]
cat("\n######################################################################\n")
cat("## Finished at", as.character(Sys.time()), "\n")
cat(sprintf("## Total wall-clock time: %.1f s (%.1f min)\n", t_wall_elapsed, t_wall_elapsed / 60))
cat("######################################################################\n")
