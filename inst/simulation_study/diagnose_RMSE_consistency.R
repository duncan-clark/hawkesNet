#!/usr/bin/env Rscript
# ==============================================================================
# RMSE Consistency Diagnostic Script
# ==============================================================================
# Purpose: Investigate why RMSE does not decrease monotonically with T.
# Run a minimal study with T = c(2, 4, 6) and verify RMSE(T=6) < RMSE(T=4) < RMSE(T=2).
# Iterate through hypotheses (filter by convergence, init, etc.) to find root cause.
# ==============================================================================

library(hawkesNet)
library(ernm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

# Run from package root
PKG_ROOT <- getwd()
if (!file.exists(file.path(PKG_ROOT, "inst"))) {
  PKG_ROOT <- dirname(dirname(sys.frame(1)$ofile))
  if (nzchar(PKG_ROOT)) setwd(PKG_ROOT)
}

# Minimal settings for fast iteration (match CS study design: T=2,4,6,10,20)
N_CORES <- 7L  # Use 7 of 8 cores
time_windows <- c(2, 4, 6, 10, 20)
N_REPS <- 12L  # More reps for less noise
TRUNCATION <- 500
MAX_ITER <- 3000
p_scale <- c(beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.1,
             CS_params1 = 1, CS_params2 = 0.1, CS_params3 = 0.1, CS_params4 = 0.1)

params_true <- list(
  mu = 10,
  beta_overall = 2,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-7, 3, 0.1, -0.1)
)

# Map true values to fit par names (net optional for formula-name fallback)
map_true_vals <- function(par_names, params, net = NULL) {
  true_vals <- setNames(numeric(length(par_names)), par_names)
  if ("mu" %in% par_names) true_vals["mu"] <- params$mu
  true_vals["beta_overall"] <- params$beta_overall
  if ("K" %in% par_names) true_vals["K"] <- params$K
  true_vals["beta_edges"] <- params$beta_edges
  true_vals["node_lambda"] <- params$node_lambda
  for (k in seq_along(params$CS_params)) {
    nm <- paste0("CS_params", k)
    if (nm %in% par_names) true_vals[nm] <- params$CS_params[k]
  }
  if (!is.null(net)) {
    exp_cs <- tryCatch(expected_params_PMF_mark_CS(net, "edges + triangles + star(c(2,3))"), error = function(e) NULL)
    if (!is.null(exp_cs$CS_params_names) && length(params$CS_params) >= length(exp_cs$CS_params_names)) {
      for (k in seq_along(exp_cs$CS_params_names)) {
        if (exp_cs$CS_params_names[k] %in% par_names) true_vals[exp_cs$CS_params_names[k]] <- params$CS_params[k]
      }
    }
  }
  true_vals
}

run_one_simfit <- function(i, curr_time, use_near_true_init = TRUE) {
  sim_res <- tryCatch({
    sim_hawkesNet(params = params_true,
                  time_window = c(0, curr_time),
                  PMF_mark = PMF_mark_CS,
                  cond_intensity = cond_intensity,
                  hashed_edges = TRUE,
                  mu_multiplier = 3,
                  verbose = FALSE,
                  truncation = TRUNCATION,
                  formula_RHS = "edges + triangles + star(c(2,3))",
                  mark_decay = "node_entrance",
                  growth_only = FALSE)
  }, error = function(e) return(NULL))
  if (is.null(sim_res)) return(NULL)

  if (use_near_true_init) {
    params_init <- list(
      mu = params_true$mu,
      beta_overall = max(0.1, params_true$beta_overall * exp(rnorm(1, 0, 0.2))),
      K = params_true$K,
      beta_edges = max(0.1, params_true$beta_edges * exp(rnorm(1, 0, 0.2))),
      node_lambda = max(0.1, params_true$node_lambda * exp(rnorm(1, 0, 0.2))),
      CS_params = params_true$CS_params + rnorm(length(params_true$CS_params), 0, 0.5)
    )
    params_init$CS_params[!is.finite(params_init$CS_params)] <- params_true$CS_params[!is.finite(params_init$CS_params)]
  } else {
    params_init <- list(
      mu = 10, beta_overall = 1, K = 0.5, beta_edges = 1, node_lambda = 1,
      CS_params = c(-7, 0, 0, 0)
    )
  }

  fit_res <- tryCatch({
    fit_hawkesNet(params_init = params_init,
                  time_window = c(0, curr_time),
                  mark_filtration = sim_res$net,
                  PMF_mark = PMF_mark_CS,
                  formula_RHS = "edges + triangles + star(c(2,3))",
                  maxit = MAX_ITER,
                  truncation = TRUNCATION,
                  mark_decay = "node_entrance",
                  growth_only = FALSE,
                  cache_intensity = TRUE,
                  combine_intensity = TRUE,
                  verbose = FALSE,
                  fixed_params = c("K", "mu"),
                  parscale = p_scale,
                  cores = 1L,
                  method = "Nelder-Mead")
  }, error = function(e) return(NULL))
  if (is.null(fit_res) || is.null(fit_res$fit)) return(NULL)

  keep <- (!is.null(fit_res$fit) &&
           length(fit_res$fit) > 0 &&
           fit_res$fit$convergence == 0 &&
           all(is.finite(fit_res$fit$par)) &&
           !any(fit_res$fit$par > 100) &&
           length(fit_res$fit$par) >= 2)
  # beta_overall check (par[2] may be beta_overall depending on order)
  par_names <- names(fit_res$fit$par)
  idx_bo <- which(par_names == "beta_overall")
  if (length(idx_bo) > 0 && fit_res$fit$par[idx_bo] > 10) keep <- FALSE

  true_vals <- map_true_vals(par_names, params_true, sim_res$net)
  n_events <- length(sim_res$events$t)

  data.frame(
    keep = keep,
    sim_id = i,
    time_window = curr_time,
    n_events = n_events,
    convergence = fit_res$fit$convergence,
    param = par_names,
    estimate = as.numeric(fit_res$fit$par),
    true_value = as.numeric(true_vals)
  )
}

cat("=== RMSE Consistency Diagnostic ===\n")
cat("Time windows:", paste(time_windows, collapse = ", "), "| Reps per window:", N_REPS, "| Cores:", N_CORES, "\n\n")

t_start <- proc.time()
all_results <- data.frame()
for (curr_time in time_windows) {
  cat("T =", curr_time, "... ")
  res_list <- mclapply(seq_len(N_REPS), function(i) run_one_simfit(i, curr_time, use_near_true_init = TRUE),
                      mc.cores = N_CORES, mc.preschedule = TRUE)
  res_list <- res_list[!sapply(res_list, is.null)]
  if (length(res_list) == 0) {
    cat("no successful fits\n")
    next
  }
  res_df <- do.call(rbind, res_list)
  res_df$time_window <- curr_time
  all_results <- rbind(all_results, res_df)
  n_conv <- sum(res_df %>% filter(!duplicated(paste(sim_id, time_window))) %>% pull(keep))
  cat(length(unique(res_df$sim_id)), "sims,", n_conv, "converged\n")
}

if (nrow(all_results) == 0) {
  stop("No results. Check sim/fit for errors.")
}

# Summary: RMSE with all runs vs converged only
summary_all_full <- all_results %>%
  group_by(time_window, param) %>%
  summarise(rmse_all = sqrt(mean((estimate - true_value)^2)), .groups = "drop")

# Per-sim keep for correct rmse_conv (converged fits only)
sim_keep <- all_results %>%
  group_by(time_window, sim_id) %>%
  summarise(keep = first(keep), .groups = "drop")
results_conv <- all_results %>%
  inner_join(sim_keep %>% filter(keep) %>% select(time_window, sim_id), by = c("time_window", "sim_id"))
summary_conv_correct <- results_conv %>%
  group_by(time_window, param) %>%
  summarise(rmse_conv = sqrt(mean((estimate - true_value)^2)), n_conv = n(), .groups = "drop")

cat("\n--- RMSE by T (all runs) ---\n")
print(summary_all_full %>% select(time_window, param, rmse_all) %>% pivot_wider(names_from = time_window, values_from = rmse_all))

cat("\n--- RMSE by T (converged only) ---\n")
print(summary_conv_correct %>% select(time_window, param, rmse_conv) %>% pivot_wider(names_from = time_window, values_from = rmse_conv))

# Check monotonicity: RMSE should decrease as T increases
check_monotonic <- function(rmse_df, rmse_col = "rmse_conv") {
  params <- unique(rmse_df$param)
  tw <- sort(unique(rmse_df$time_window))
  ok <- TRUE
  for (p in params) {
    vals <- rmse_df %>% filter(param == p) %>% arrange(time_window) %>% pull(!!sym(rmse_col))
    if (length(vals) < 2) next
    for (j in 2:length(vals)) {
      if (vals[j] > vals[j-1] * 1.01) {  # allow 1% tolerance
        cat("  NON-MONOTONIC:", p, "T=", tw[j], "rmse=", round(vals[j], 4), "> T=", tw[j-1], "rmse=", round(vals[j-1], 4), "\n")
        ok <- FALSE
      }
    }
  }
  if (ok) cat("  All params monotonic (or near)\n")
  ok
}

cat("\n--- Monotonicity check (converged) ---\n")
check_monotonic(summary_conv_correct, "rmse_conv")

cat("\n--- Monotonicity check (all runs) ---\n")
check_monotonic(summary_all_full, "rmse_all")

# Convergence rate by T
cat("\n--- Convergence rate by T ---\n")
conv_rate <- all_results %>%
  group_by(time_window, sim_id) %>%
  summarise(keep = first(keep), .groups = "drop") %>%
  group_by(time_window) %>%
  summarise(n_total = n(), n_conv = sum(keep), rate = mean(keep), .groups = "drop")
print(conv_rate)

# Event count by T
cat("\n--- Mean events per T ---\n")
evt_summary <- all_results %>%
  group_by(time_window, sim_id) %>%
  summarise(n_events = first(n_events), .groups = "drop") %>%
  group_by(time_window) %>%
  summarise(mean_events = mean(n_events), sd_events = sd(n_events), .groups = "drop")
print(evt_summary)

# Elapsed time for title
elapsed_sec <- (proc.time() - t_start)[3]
elapsed_str <- if (elapsed_sec >= 3600) {
  sprintf("%.1f h", elapsed_sec / 3600)
} else if (elapsed_sec >= 60) {
  sprintf("%.1f min", elapsed_sec / 60)
} else {
  sprintf("%.0f s", elapsed_sec)
}

# Plot (same format as CS study for direct comparison)
summary_for_plot <- summary_conv_correct %>% rename(rmse = rmse_conv)
p_rmse <- ggplot(summary_for_plot, aes(x = time_window, y = rmse)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~param, scales = "free_y") +
  labs(title = sprintf("RMSE Decay as Data Increases (CS model) — run time: %s", elapsed_str),
       x = "Time Window Length (T)",
       y = "Root Mean Squared Error") +
  theme_bw()
print(p_rmse)

# Save for inspection
out_dir <- file.path(PKG_ROOT, "cluster_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(
  all_results = all_results,
  summary_conv = summary_conv_correct,
  summary_all = summary_all_full,
  conv_rate = conv_rate
), file.path(out_dir, "diagnose_RMSE_results.RDS"))
ggsave(file.path(out_dir, "diagnose_RMSE_plot.png"), plot = p_rmse, width = 10, height = 8, dpi = 150)
cat("\nSaved to", file.path(out_dir, "diagnose_RMSE_results.RDS"), "\n")
cat("Saved RMSE plot to", file.path(out_dir, "diagnose_RMSE_plot.png"), "\n")

# ==============================================================================
# SAME-REALIZATION TRUNCATION DESIGN (proper consistency test)
# ==============================================================================
# Simulate long processes (T=20), then fit on truncated data at T=2,4,6,10,20.
# Same underlying process -> more data should give lower RMSE.
# ==============================================================================
cat("\n\n=== Same-Realization Truncation Design ===\n")
cat("Simulate at T=20, fit on truncated data at T=2,4,6,10,20 (same process, more data)\n")

T_MAX <- 20
N_TRUNC_REPS <- 8L
trunc_windows <- c(2, 4, 6, 10, 20)

run_one_trunc_realization <- function(i) {
  sim_full <- tryCatch({
    sim_hawkesNet(params = params_true,
                  time_window = c(0, T_MAX),
                  PMF_mark = PMF_mark_CS,
                  cond_intensity = cond_intensity,
                  hashed_edges = TRUE,
                  mu_multiplier = 3,
                  verbose = FALSE,
                  truncation = TRUNCATION,
                  formula_RHS = "edges + triangles + star(c(2,3))",
                  mark_decay = "node_entrance",
                  growth_only = FALSE)
  }, error = function(e) NULL)
  if (is.null(sim_full)) return(NULL)
  rows <- list()
  for (curr_t in trunc_windows) {
    net_t <- filtration_to_net(sim_full$net, curr_t, equals = TRUE)
    n_evt <- length(get_times(net_t)$times)
    if (n_evt < 5) next
    params_init <- list(
      mu = params_true$mu * exp(rnorm(1, 0, 0.2)),
      beta_overall = max(0.1, params_true$beta_overall * exp(rnorm(1, 0, 0.2))),
      K = params_true$K,
      beta_edges = max(0.1, params_true$beta_edges * exp(rnorm(1, 0, 0.2))),
      node_lambda = max(0.1, params_true$node_lambda * exp(rnorm(1, 0, 0.2))),
      CS_params = params_true$CS_params + rnorm(4, 0, 0.5)
    )
    fit_t <- tryCatch({
      fit_hawkesNet(params_init = params_init,
                    time_window = c(0, curr_t),
                    mark_filtration = net_t,
                    PMF_mark = PMF_mark_CS,
                    formula_RHS = "edges + triangles + star(c(2,3))",
                    maxit = MAX_ITER,
                    truncation = TRUNCATION,
                    mark_decay = "node_entrance",
                    growth_only = FALSE,
                    cache_intensity = TRUE,
                    combine_intensity = TRUE,
                    verbose = FALSE,
                    fixed_params = c("K", "mu"),
                    parscale = p_scale,
                    cores = 1L,
                    method = "Nelder-Mead")
    }, error = function(e) NULL)
    if (is.null(fit_t) || is.null(fit_t$fit)) next
    keep_t <- (fit_t$fit$convergence == 0 && all(is.finite(fit_t$fit$par)) && !any(fit_t$fit$par > 100))
    par_names <- names(fit_t$fit$par)
    true_vals <- map_true_vals(par_names, params_true, net_t)
    rows[[length(rows) + 1]] <- data.frame(
      rep_id = i, time_window = curr_t, keep = keep_t, n_events = n_evt,
      param = par_names, estimate = as.numeric(fit_t$fit$par),
      true_value = as.numeric(true_vals)
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

cat("  Running", N_TRUNC_REPS, "realizations in parallel (", N_CORES, "cores)...\n")
trunc_list <- mclapply(seq_len(N_TRUNC_REPS), run_one_trunc_realization,
                      mc.cores = N_CORES, mc.preschedule = TRUE)
trunc_valid <- trunc_list[!sapply(trunc_list, is.null)]
trunc_results <- if (length(trunc_valid) > 0) do.call(rbind, trunc_valid) else data.frame()

if (nrow(trunc_results) > 0) {
  trunc_conv <- trunc_results %>% filter(keep == TRUE)
  if (nrow(trunc_conv) > 0) {
    trunc_summary <- trunc_conv %>%
      group_by(time_window, param) %>%
      summarise(rmse = sqrt(mean((estimate - true_value)^2)), n = n(), .groups = "drop")
    cat("\n--- RMSE by T (same realization, converged) ---\n")
    print(trunc_summary %>% select(time_window, param, rmse) %>%
            pivot_wider(names_from = time_window, values_from = rmse))
    cat("\n--- Monotonicity (same-realization) ---\n")
    check_monotonic(trunc_summary, "rmse")
  }
}
