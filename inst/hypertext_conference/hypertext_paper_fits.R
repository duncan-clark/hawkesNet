## =============================================================================
## Hypertext 2009 conference: Paper Results Script
## =============================================================================
## 6 fits total:
##   1. Day 1 — triangles + star(c(2,3)) + degree(0)
##   2. Day 1 — triangles + star(c(2,3)) + degree(0) [Non-Simple placeholder]
##   3. Day 1 — gwesp(0.5) + gwdegree(0.5) + degree(0)
##   4. Day 1 — gwesp(0.5) + gwdegree(0.5) + degree(0) [Non-Simple placeholder]
##   5. Full conference — inhomogeneous, triangles + star(c(2,3)) + degree(0)
##   6. Full conference — inhomogeneous, triangles + star(c(2,3)) + degree(0) [Non-Simple placeholder]
##
## All fits: growth_only = FALSE, no truncation, no fixed params, include GOF.
## =============================================================================

library(hawkesNet)
library(dplyr)
library(network)
library(ernm)
library(sna)
library(ggplot2)

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
LOCAL_QUICK <- isTRUE(as.logical(Sys.getenv("LOCAL_QUICK", "FALSE")))
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 100L))
MAX_ITER <- as.integer(Sys.getenv("MAX_ITER", if (LOCAL_QUICK) 100L else 5000L))
N_GOF <- as.integer(Sys.getenv("N_GOF", if (LOCAL_QUICK) 2L else 100L))
N_GOF_OUTER <- as.integer(Sys.getenv("GOF_CORES_OUTER",
                                      if (N_CORES >= 100L) 100L else 0L))

MARK_DECAY <- "node_entrance"
GROWTH_ONLY <- FALSE

# Paths
PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) {
  Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  getwd()
}
OUTPUT_DIR <- file.path(PKG_ROOT, "cluster_output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Config: N_CORES=%d | MAX_ITER=%d | N_GOF=%d | GOF_OUTER=%d\n",
            N_CORES, MAX_ITER, N_GOF, N_GOF_OUTER))

# =============================================================================
# Helpers
# =============================================================================

make_hypertext_net <- function(df) {
  swap <- df$from > df$to
  df[swap, c("from", "to")] <- df[swap, c("to", "from")]
  df <- df %>% distinct() %>% arrange(time)

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

  list(net = net, edges = df)
}

subset_first_session <- function(df, gap_threshold = 1.0) {
  df <- df[order(df$time), , drop = FALSE]
  gaps <- diff(df$time)
  first_gap <- which(gaps > gap_threshold)[1]
  if (is.na(first_gap)) return(df)
  df[seq_len(first_gap), , drop = FALSE]
}

# =============================================================================
# Data Loading
# =============================================================================
cat("\n=== Loading Data ===\n")
raw <- read.table(system.file("extdata", "ht09_contact_list.dat",
                               package = "hawkesNet"))
df_all <- data.frame(
  time = raw$V1 / 3600,
  from = raw$V2,
  to   = raw$V3
)
df_all$time <- df_all$time - min(df_all$time)
df_all <- df_all[order(df_all$time), ]

# --- Day 1 network ---
df_day1 <- subset_first_session(df_all)
obj_day1 <- make_hypertext_net(df_day1)
net_day1 <- normalize_times_01(obj_day1$net)
times_day1 <- get_times(net_day1)$times
tw_day1 <- c(min(times_day1), max(times_day1))
cat(sprintf("Day 1: %d nodes, %d edges, %d events, T=[%.4f, %.4f]\n",
            network.size(net_day1), network.edgecount(net_day1),
            length(times_day1), tw_day1[1], tw_day1[2]))

# --- Full conference network ---
obj_full <- make_hypertext_net(df_all)
net_full <- normalize_times_01(obj_full$net)
times_full <- get_times(net_full)$times
tw_full <- c(min(times_full), max(times_full))
cat(sprintf("Full:  %d nodes, %d edges, %d events, T=[%.4f, %.4f]\n",
            network.size(net_full), network.edgecount(net_full),
            length(times_full), tw_full[1], tw_full[2]))

# --- Inhomogeneous background for full conference ---
cat("Estimating inhomogeneous background (KDE) for full conference...\n")
inhom_bg <- prepare_inhomogeneous_background(net_full)
cat(sprintf("  KDE bandwidth: %.4f | integral_bg: %.2f\n",
            inhom_bg$mu_fit$bw, inhom_bg$integral_bg))

# =============================================================================
# Fitting Logic
# =============================================================================

run_fit <- function(net, time_window, formula_rhs, label,
                    mu_vec = NULL, integral_bg = NULL) {
  cat(sprintf("\n>>> Fit: %s <<<\n", label))

  exp_cs <- expected_params_PMF_mark_CS(net, formula_rhs)
  n_cs <- exp_cs$CS_params_length
  n_events <- length(get_times(net)$times)
  TRUNC <- network.size(net)

  m_init <- network.edgecount(net) / n_events

  params_init <- list(
    mu = n_events / (time_window[2] - time_window[1]),
    beta_overall = 0.3,
    K = 0.5,
    beta_edges = 0.3,
    node_lambda = 0.5,
    m = m_init,
    CS_params = c(-5, rep(0, n_cs - 1))
  )

  p_scale <- c(
    mu = 1, beta_overall = 0.1, K = 0.1, beta_edges = 0.1,
    node_lambda = 0.5, m = 0.5,
    setNames(rep(0.1, n_cs), paste0("CS_params", seq_len(n_cs)))
  )

  use_inhom <- !is.null(mu_vec)

  t_fit_start <- proc.time()
  fit <- tryCatch({
    if (use_inhom) {
      fit_hawkesNet_inhom(
        params_init = params_init,
        time_window = time_window,
        mark_filtration = net,
        PMF_mark = PMF_mark_CS,
        mu_vec = mu_vec,
        integral_bg = integral_bg,
        maxit = MAX_ITER,
        fixed_params = NULL,
        parscale = p_scale,
        cores = N_CORES,
        cache_intensity = TRUE,
        combine_intensity = TRUE,
        verbose = TRUE,
        formula_RHS = formula_rhs,
        truncation = TRUNC,
        mark_decay = MARK_DECAY,
        growth_only = GROWTH_ONLY
      )
    } else {
      fit_hawkesNet(
        params_init = params_init,
        time_window = time_window,
        mark_filtration = net,
        PMF_mark = PMF_mark_CS,
        formula_RHS = formula_rhs,
        truncation = TRUNC,
        mark_decay = MARK_DECAY,
        growth_only = GROWTH_ONLY,
        maxit = MAX_ITER,
        fixed_params = NULL,
        parscale = p_scale,
        cores = N_CORES,
        cache_intensity = TRUE,
        combine_intensity = TRUE,
        verbose = TRUE
      )
    }
  }, error = function(e) {
    cat(sprintf("  FIT FAILED: %s\n", e$message))
    NULL
  })
  t_fit <- (proc.time() - t_fit_start)[3]
  cat(sprintf("  Fit time: %.1f s\n", t_fit))

  if (!is.null(fit) && !is.null(fit$fit_table)) {
    cat("  Fit table:\n")
    print(fit$fit_table)
  }

  # GOF
  cat(sprintf("  Running GOF (%d sims, %d outer workers)...\n", N_GOF, N_GOF_OUTER))
  t_gof_start <- proc.time()
  gof_res <- tryCatch({
    gof_args <- list(
      fit = fit,
      net_obs = net,
      params_init = params_init,
      PMF_mark = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS = formula_rhs,
      time_window = time_window,
      truncation = TRUNC,
      mark_decay = MARK_DECAY,
      growth_only = GROWTH_ONLY,
      n_sim = N_GOF,
      cores_outer = N_GOF_OUTER,
      verbose = TRUE
    )
    if (use_inhom) {
      gof_args$inhom_bg <- list(
        mu_vec = mu_vec,
        integral_bg = integral_bg,
        mu_fit = inhom_bg$mu_fit,
        Lambda_fun = inhom_bg$Lambda_fun
      )
    }
    do.call(gof, gof_args)
  }, error = function(e) {
    cat(sprintf("  GOF FAILED: %s\n", e$message))
    NULL
  })
  t_gof <- (proc.time() - t_gof_start)[3]
  cat(sprintf("  GOF time: %.1f s\n", t_gof))

  list(fit = fit, gof = gof_res, label = label,
       time_fit = t_fit, time_gof = t_gof)
}

# =============================================================================
# Execute 6 Fits
# =============================================================================

results <- list()

# 1. Day 1 — Triangles + Star + Degree(0)
results$day1_tri_simple <- run_fit(
  net_day1, tw_day1, 
  formula_rhs = "triangles + star(c(2,3)) + degree(0)",
  label = "Day1-Tri-Simple"
)

# 2. Day 1 — Triangles + Star + Degree(0) [Placeholder for non-simple]
results$day1_tri_nonsimple <- run_fit(
  net_day1, tw_day1, 
  formula_rhs = "triangles + star(c(2,3)) + degree(0)",
  label = "Day1-Tri-NonSimple"
)

# 3. Day 1 — GWESP + GWDegree + Degree(0)
results$day1_gw_simple <- run_fit(
  net_day1, tw_day1, 
  formula_rhs = "gwesp(0.5) + gwdegree(0.5) + degree(0)",
  label = "Day1-GW-Simple"
)

# 4. Day 1 — GWESP + GWDegree + Degree(0) [Placeholder for non-simple]
results$day1_gw_nonsimple <- run_fit(
  net_day1, tw_day1, 
  formula_rhs = "gwesp(0.5) + gwdegree(0.5) + degree(0)",
  label = "Day1-GW-NonSimple"
)

# 5. Full conference — Inhomogeneous, Triangles + Star + Degree(0)
results$full_tri_simple <- run_fit(
  net_full, tw_full, 
  formula_rhs = "triangles + star(c(2,3)) + degree(0)",
  label = "Full-Inhom-Tri-Simple",
  mu_vec = inhom_bg$mu_vec,
  integral_bg = inhom_bg$integral_bg
)

# 6. Full conference — Inhomogeneous, Triangles + Star + Degree(0) [Non-Simple]
results$full_tri_nonsimple <- run_fit(
  net_full, tw_full, 
  formula_rhs = "triangles + star(c(2,3)) + degree(0)",
  label = "Full-Inhom-Tri-NonSimple",
  mu_vec = inhom_bg$mu_vec,
  integral_bg = inhom_bg$integral_bg
)

# =============================================================================
# Paper Output Section
# =============================================================================
cat("\n######################################################################\n")
cat("## GENERATING PAPER OUTPUTS\n")
cat("######################################################################\n")

# 1. Summary Table
summary_table <- do.call(rbind, lapply(results, function(res) {
  if (is.null(res$fit) || is.null(res$fit$fit_table)) return(NULL)
  tab <- res$fit$fit_table
  tab$model <- res$label
  tab
}))
if (!is.null(summary_table)) {
  write.csv(summary_table, file.path(OUTPUT_DIR, "paper_fit_summary.csv"),
            row.names = FALSE)
  cat("\nFit Summary:\n")
  print(summary_table)
}

# 2. GOF Plots
for (name in names(results)) {
  res <- results[[name]]
  if (is.null(res$gof) || is.null(res$gof$plots)) next
  if (!is.null(res$gof$plots$waiting_times_plot)) {
    ggsave(
      filename = file.path(OUTPUT_DIR, paste0("gof_wait_", name, ".pdf")),
      plot = res$gof$plots$waiting_times_plot,
      width = 10, height = 8
    )
  }
  if (!is.null(res$gof$plots$degree_plot)) {
    ggsave(
      filename = file.path(OUTPUT_DIR, paste0("gof_deg_", name, ".pdf")),
      plot = res$gof$plots$degree_plot,
      width = 8, height = 6
    )
  }
}

# 3. Save full results
saveRDS(results, file.path(OUTPUT_DIR, "paper_results_full.RDS"))

cat("\nDone. Results saved to:", OUTPUT_DIR, "\n")
