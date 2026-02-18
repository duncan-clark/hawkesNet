dat <- results_hypertext_full

t0 <- get_times(dat$net)
t1 <- get_times(dat$gof$nets_sim[[1]])
t2 <- get_times(dat$gof_decay05$nets_sim[[1]])
t3 <- get_times(dat$gof_star_esp$nets_sim[[1]])

hist(t0$times)
hist(t1$times)
hist(t2$times)
hist(t3$times)

# lets simulate some networks with the fitted params:
res <- results_hypertext_full

# ------------------------------------------------------------------
# Simulate from the fitted hypertext model — all 5 fits
# ------------------------------------------------------------------
# After loading the results:
#   res <- readRDS("cluster_output/results_hypertext_full.RDS")
#   # or: res <- results_hypertext_full  (if loaded interactively)

library(hawkesNet)
library(network)
library(sna)

# Load results
res <- readRDS("/nesi/project/uoo04008/Duncan/hawkes_net/hawkesNet/cluster_output/results_hypertext_full.RDS")

formula_rhs <- res$formula_rhs
SEED_EVENTS <- 50
tw          <- c(0, 1)

# =============================================================================
# Helper: seed + simulate
# =============================================================================
sim_block <- function(fit, net_obs, inhom_bg, trunc, mark_decay, label,
                      seed_events = SEED_EVENTS, mu_mult = 2) {
  cat(sprintf("\n--- Simulating from %s ---\n", label))
  pfit <- fit$params

  seed_net   <- NULL
  seed_times <- NULL
  all_times  <- get_times(net_obs)$times
  if (seed_events > 0 && length(all_times) >= seed_events) {
    t_seed     <- all_times[seed_events]
    seed_net   <- filtration_to_net(net_obs, t_seed, equals = TRUE)
    seed_times <- all_times[1:seed_events]
    cat("  Seeding with first", seed_events, "events (up to t =", round(t_seed, 4), ")\n")
  }

  sim <- tryCatch({
    sim_hawkesNet(
      params         = pfit,
      time_window    = tw,
      PMF_mark       = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS    = formula_rhs,
      truncation     = trunc,
      mark_decay     = mark_decay,
      growth_only    = FALSE,
      max_node_time  = Inf,
      hashed_edges   = TRUE,
      verbose        = 2,
      mu_multiplier  = mu_mult,
      inhom_bg       = inhom_bg,
      seed_net       = seed_net,
      seed_times     = seed_times
    )
  }, error = function(e) { cat("  FAILED:", e$message, "\n"); NULL })

  if (!is.null(sim) && !is.null(sim$net)) {
    deg <- sna::degree(sim$net, gmode = "graph")
    cat(sprintf("  %s: %d nodes, %d edges, %d events, mean deg(>0) = %.2f\n",
                label, network.size(sim$net), network.edgecount(sim$net),
                length(sim$events$t), mean(deg[deg > 0], na.rm = TRUE)))
  }
  sim
}

# =============================================================================
# BLOCK 1: Fit 1 — Simple, activity decay
# =============================================================================
sim_1 <- sim_block(
  fit       = res$fit_day1,
  net_obs   = res$net_day1,
  inhom_bg  = res$inhom_bg_day1,
  trunc     = res$truncation_day1,
  mark_decay = "activity",
  label     = "Fit 1 (Simple, activity)"
)

# =============================================================================
# BLOCK 2: Fit 1b — Non-Simple, activity decay
# =============================================================================
sim_1b <- sim_block(
  fit       = res$fit_nonsimple,
  net_obs   = res$net_nonsimple,
  inhom_bg  = res$inhom_bg_nonsimple,
  trunc     = res$truncation_nonsimple,
  mark_decay = "activity",
  label     = "Fit 1b (Non-Simple, activity)",
  mu_mult   = 5
)

# =============================================================================
# BLOCK 3: Fit 3 — Simple, node_entrance decay
# =============================================================================
sim_3 <- sim_block(
  fit       = res$fit_ne_simple,
  net_obs   = res$net_day1,
  inhom_bg  = res$inhom_bg_day1,
  trunc     = res$truncation_day1,
  mark_decay = "node_entrance",
  label     = "Fit 3 (Simple, node_entrance)"
)

# =============================================================================
# BLOCK 4: Fit 4 — Non-Simple, node_entrance decay
# =============================================================================
sim_4 <- sim_block(
  fit       = res$fit_ne_nonsimple,
  net_obs   = res$net_nonsimple,
  inhom_bg  = res$inhom_bg_nonsimple,
  trunc     = res$truncation_nonsimple,
  mark_decay = "node_entrance",
  label     = "Fit 4 (Non-Simple, node_entrance)",
  mu_mult   = 5
)

# =============================================================================
# BLOCK 5: Fit 5 — Non-Simple, activity decay, m-parameter (Poisson edge count)
# =============================================================================
sim_5 <- sim_block(
  fit       = res$fit_m,
  net_obs   = res$net_nonsimple,
  inhom_bg  = res$inhom_bg_nonsimple,
  trunc     = res$truncation_nonsimple,
  mark_decay = "activity",
  label     = "Fit 5 (Non-Simple, activity+m)",
  mu_mult   = 5
)

# =============================================================================
# SUMMARY TABLE
# =============================================================================
cat("\n======================================================================\n")
cat("                        OBSERVED vs SIMULATED\n")
cat("======================================================================\n")

print_row <- function(label, net_obs, times_obs, sim_obj) {
  n_obs  <- network.size(net_obs)
  e_obs  <- network.edgecount(net_obs)
  ev_obs <- length(times_obs)
  if (is.null(sim_obj) || is.null(sim_obj$net)) {
    cat(sprintf("  %-35s  Obs: %3d n %4d e %4d ev  |  Sim: FAILED\n",
                label, n_obs, e_obs, ev_obs))
  } else {
    n_sim  <- network.size(sim_obj$net)
    e_sim  <- network.edgecount(sim_obj$net)
    ev_sim <- length(sim_obj$events$t)
    cat(sprintf("  %-35s  Obs: %3d n %4d e %4d ev  |  Sim: %3d n %4d e %4d ev\n",
                label, n_obs, e_obs, ev_obs, n_sim, e_sim, ev_sim))
  }
}

times_1  <- get_times(res$net_day1)$times
times_1b <- get_times(res$net_nonsimple)$times

print_row("Fit 1  (simple, activity)",       res$net_day1,      times_1,  sim_1)
print_row("Fit 1b (non-simple, activity)",   res$net_nonsimple, times_1b, sim_1b)
print_row("Fit 3  (simple, node_entrance)",  res$net_day1,      times_1,  sim_3)
print_row("Fit 4  (non-simple, node_entr.)", res$net_nonsimple, times_1b, sim_4)
print_row("Fit 5  (non-simple, activity+m)", res$net_nonsimple, times_1b, sim_5)
cat("======================================================================\n")

# Print fitted m value if available
if (!is.null(res$fit_m)) {
  cat(sprintf("\n  Fitted m (Poisson edge rate): %.3f\n", res$fit_m$params$m))
  cat(sprintf("  Observed edges/event:         %.3f\n",
              network.edgecount(res$net_nonsimple) / length(times_1b)))
}
