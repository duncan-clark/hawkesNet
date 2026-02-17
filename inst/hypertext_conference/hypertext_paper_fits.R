## =============================================================================
## Hypertext 2009 conference: Cleaned Paper Results Script
## =============================================================================
## This script performs 4 specific m-parameter fits:
##   1. GWESP(0.5) + GWDegree(0.5) | Activity Decay
##   2. GWESP(0.5) + GWDegree(0.5) | Node Entrance Decay
##   3. Triangles + Star(c(2,3))   | Activity Decay
##   4. Triangles + Star(c(2,3))   | Node Entrance Decay
##
## No truncation is used. All fits include GOF.
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
N_CORES_GOF <- as.integer(Sys.getenv("GOF_CORES", N_CORES))
# Outer parallelism for GOF: distribute simulations across workers.
# PSOCK limit is 128; we use 50 outer workers to be safe and efficient.
N_GOF_OUTER <- as.integer(Sys.getenv("GOF_CORES_OUTER", if (N_CORES >= 100L) 50L else 0L))

# Paths
PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) Sys.getenv("SLURM_SUBMIT_DIR") else getwd()
OUTPUT_DIR <- file.path(PKG_ROOT, "cluster_output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Helpers
# =============================================================================

make_hypertext_net <- function(df) {
  # Undirected: canonical ordering (from < to) and drop duplicates
  swap <- df$from > df$to
  df[swap, c("from", "to")] <- df[swap, c("to", "from")]
  df <- df %>% distinct() %>% arrange(time)

  # Map vertex ids to 1..n in entry-time order
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
raw <- read.table(system.file("extdata", "ht09_contact_list.dat", package = "hawkesNet"))
df_ns <- data.frame(
  time = raw$V1 / 3600,
  from = raw$V2,
  to   = raw$V3
)
df_ns$time <- df_ns$time - min(df_ns$time)
df_ns <- df_ns[order(df_ns$time), ]

# Use first session only for clean results
df_day1 <- subset_first_session(df_ns)
obj_ns <- make_hypertext_net(df_day1)
net <- normalize_times_01(obj_ns$net)
times <- get_times(net)$times
time_window <- c(min(times), max(times))

# No truncation
TRUNCATION <- network.size(net)

# =============================================================================
# Fitting Logic
# =============================================================================

run_paper_fit <- function(formula_rhs, mark_decay, label) {
  cat("\n>>> Running Fit:", label, "<<<\n")
  
  exp_cs <- expected_params_PMF_mark_CS(net, formula_rhs)
  n_cs <- exp_cs$CS_params_length
  
  # Initial m from data
  m_init <- network.edgecount(net) / length(times)
  
  params_init <- list(
    mu = length(times) / (time_window[2] - time_window[1]),
    beta_overall = 0.3,
    K = 0.5,
    beta_edges = 0.3,
    node_lambda = 0.5,
    m = m_init,
    CS_params = c(-5, -3, rep(0, n_cs - 2))
  )
  
  p_scale <- c(
    mu = 1, beta_overall = 0.1, K = 0.1, beta_edges = 0.1, node_lambda = 0.5, m = 0.5,
    setNames(rep(0.1, n_cs), paste0("CS_params", seq_len(n_cs)))
  )
  
  fit <- fit_hawkesNet(
    params_init = params_init,
    time_window = time_window,
    mark_filtration = net,
    PMF_mark = PMF_mark_CS,
    formula_RHS = formula_rhs,
    truncation = TRUNCATION,
    mark_decay = mark_decay,
    growth_only = FALSE,
    maxit = MAX_ITER,
    fixed_params = NULL, # m-parameter model: all params free
    parscale = p_scale,
    cores = N_CORES,
    cache_intensity = TRUE,
    combine_intensity = TRUE,
    verbose = TRUE
  )
  
  gof_res <- gof(
    fit = fit,
    net_obs = net,
    params_init = params_init,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = formula_rhs,
    time_window = time_window,
    truncation = TRUNCATION,
    mark_decay = mark_decay,
    n_sim = N_GOF,
    cores_outer = N_GOF_OUTER,
    verbose = TRUE
  )
  
  list(fit = fit, gof = gof_res, label = label)
}

# =============================================================================
# Execute 4 Fits
# =============================================================================

results <- list()

# 1. GWESP + GWDegree | Activity
results$gw_activity <- run_paper_fit(
  formula_rhs = "gwesp(0.5, fixed=TRUE) + gwdegree(0.5, fixed=TRUE)",
  mark_decay = "activity",
  label = "GW-Activity"
)

# 2. GWESP + GWDegree | Node Entrance
results$gw_ne <- run_paper_fit(
  formula_rhs = "gwesp(0.5, fixed=TRUE) + gwdegree(0.5, fixed=TRUE)",
  mark_decay = "node_entrance",
  label = "GW-NodeEntrance"
)

# 3. Triangles + Star | Activity
results$tri_activity <- run_paper_fit(
  formula_rhs = "triangles + star(c(2,3))",
  mark_decay = "activity",
  label = "Tri-Activity"
)

# 4. Triangles + Star | Node Entrance
results$tri_ne <- run_paper_fit(
  formula_rhs = "triangles + star(c(2,3))",
  mark_decay = "node_entrance",
  label = "Tri-NodeEntrance"
)

# =============================================================================
# Paper Output Section
# =============================================================================
cat("\n######################################################################\n")
cat("## GENERATING PAPER OUTPUTS\n")
cat("######################################################################\n")

# 1. Summary Table
summary_table <- do.call(rbind, lapply(results, function(res) {
  tab <- res$fit$fit_table
  tab$model <- res$label
  tab
}))
write.csv(summary_table, file.path(OUTPUT_DIR, "paper_fit_summary.csv"), row.names = FALSE)

# 2. GOF Plots
for (name in names(results)) {
  res <- results[[name]]
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

saveRDS(results, file.path(OUTPUT_DIR, "paper_results_full.RDS"))

cat("\nDone. Results saved to:", OUTPUT_DIR, "\n")
