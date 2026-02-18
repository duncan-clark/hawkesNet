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
# Simulate from the fitted hypertext model
# ------------------------------------------------------------------
# After loading the results:
#   res <- readRDS("cluster_output/results_hypertext_full.RDS")
#   # or: res <- results_hypertext_full  (if loaded interactively)

library(hawkesNet)
library(network)
library(sna)

# --- Choose which fit to simulate from: day-1 or full ---
# Day-1 (first session, homogeneous background):
fit       <- res$fit_day1
net_obs   <- res$net_day1
inhom_bg  <- res$inhom_bg_day1
tw        <- res$time_window_day1
trunc     <- res$truncation_day1

# # Full data (all 3 days, inhomogeneous KDE background):
# fit       <- res$fit_full
# net_obs   <- res$net_full
# inhom_bg  <- res$inhom_bg_full
# tw        <- res$time_window_full
# trunc     <- res$truncation_full

formula_rhs <- res$formula_rhs   # "edges + degree(0) + gwdegree(0.1) + gwesp(0.1)"

# --- Extract fitted parameters ---
pfit <- fit$params

# --- (Optional) Seed the simulation with observed data ---
# Seeding conditions the simulation on the first N observed events,
# then simulates forward. Set to 0 for unconditional simulation.
SEED_EVENTS <- 50
seed_net   <- NULL
seed_times <- NULL

if (SEED_EVENTS > 0) {
  all_times <- get_times(net_obs)$times
  if (length(all_times) >= SEED_EVENTS) {
    t_seed     <- all_times[SEED_EVENTS]
    seed_net   <- filtration_to_net(net_obs, t_seed, equals = TRUE)
    seed_times <- all_times[1:SEED_EVENTS]
    cat("Seeding with first", SEED_EVENTS, "events (up to t =", round(t_seed, 4), ")\n")
  }
}

# --- Simulate ---
N_SIMS <- 1
sims <- vector("list", N_SIMS)

for (i in seq_len(N_SIMS)) {
  cat("Simulation", i, "/", N_SIMS, "... ")
  sims[[i]] <- tryCatch({
    sim_hawkesNet(
      params         = pfit,
      time_window    = tw,
      PMF_mark       = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS    = formula_rhs,
      truncation     = trunc,
      mark_decay     = "activity",
      growth_only    = FALSE,
      max_node_time  = max(get_times(net_obs)$node_times),
      hashed_edges   = TRUE,
      verbose        = TRUE,
      mu_multiplier  = 2,
      stop_on_full_network = FALSE,
      inhom_bg       = inhom_bg,
      seed_net       = seed_net,
      seed_times     = seed_times
    )
  }, error = function(e) {
    cat("FAILED:", e$message, "\n")
    list(net = NULL, events = list(t = numeric(0)), error = e$message)
  })
  if (!is.null(sims[[i]]$net)) {
    cat("OK —",
        network.size(sims[[i]]$net), "nodes,",
        network.edgecount(sims[[i]]$net), "edges,",
        length(sims[[i]]$events$t), "events\n")
  }
}

# --- Quick comparison to observed ---
cat("\n--- Observed vs Simulated ---\n")
cat("Observed:", network.size(net_obs), "nodes,",
    network.edgecount(net_obs), "edges ",
    length(get_times(net_obs)$times), "events\n")

for (i in seq_len(N_SIMS)) {
  s <- sims[[i]]
  if (!is.null(s$net)) {
    deg <- sna::degree(s$net, gmode = "graph")
    cat(sprintf("Sim %d: %d nodes, %d edges, mean deg(>0) = %.2f\n",
                i, network.size(s$net), network.edgecount(s$net),
                mean(deg[deg > 0], na.rm = TRUE)))
  }
}
