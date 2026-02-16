## Quick verification that the hypertext pipeline works end-to-end
## with the [0,1] time normalization. Uses maxit=3 and cores=1 to keep
## runtime under ~5 min on a laptop.
library(hawkesNet)
library(dplyr)
library(network)

cat("=== Hypertext pipeline verification ===\n")

# ---- 1. Load data (V1 is in seconds with 20-s resolution) ----
raw <- read.table(system.file("extdata", "ht09_contact_list.dat", package = "hawkesNet"))
df <- data.frame(time = raw$V1 / 3600, from = raw$V2, to = raw$V3)

set.seed(1)
df$time <- df$time + rnorm(nrow(df), 0, 0.01 / 3600)
df$time <- df$time - min(df$time)
df <- df[order(df$time), ]
cat("  Raw data: ", nrow(df), " rows, time range", round(range(df$time), 3), "hours\n")

# ---- 2. Helpers (mirroring the main script) ----
subset_first_session <- function(df, gap_threshold = 1.0) {
  df <- df[order(df$time), , drop = FALSE]
  gaps <- diff(df$time)
  first_gap <- which(gaps > gap_threshold)[1]
  if (is.na(first_gap)) return(df)
  df[seq_len(first_gap), , drop = FALSE]
}

make_hypertext_net <- function(df, use_first_contact_only = TRUE) {
  swap <- df$from > df$to
  df[swap, c("from", "to")] <- df[swap, c("to", "from")]
  df <- df %>% distinct()
  df$time <- as.numeric(df$time)
  df <- df[is.finite(df$time), ]
  df <- df[order(df$time), , drop = FALSE]
  if (use_first_contact_only) {
    df <- df %>% group_by(.data$from, .data$to) %>%
      summarise(time = min(.data$time), .groups = "drop")
  }
  df <- df %>% arrange(time)
  nodes_raw <- sort(unique(c(df$from, df$to)))
  entry <- vapply(nodes_raw, function(v) min(df$time[df$from == v | df$to == v], na.rm = TRUE), numeric(1))
  ord <- order(entry, nodes_raw)
  nodes <- nodes_raw[ord]
  id_map <- setNames(seq_along(nodes), nodes)
  df$tail <- unname(id_map[as.character(df$from)])
  df$head <- unname(id_map[as.character(df$to)])
  el <- as.matrix(df[, c("tail", "head")])
  net <- network::network(el, matrix.type = "edgelist", directed = FALSE)
  network::set.edge.attribute(net, "time", df$time)
  network::set.vertex.attribute(net, "time", entry)
  list(net = net, edges = df)
}

find_gap_intervals <- function(event_times, gap_threshold = 1.0) {
  event_times <- sort(event_times)
  gaps <- diff(event_times)
  big <- which(gaps > gap_threshold)
  if (length(big) == 0) return(data.frame(start = numeric(0), end = numeric(0)))
  data.frame(start = event_times[big], end = event_times[big + 1L])
}

in_gap <- function(t, gap_intervals) {
  if (nrow(gap_intervals) == 0) return(rep(FALSE, length(t)))
  out <- logical(length(t))
  for (i in seq_len(nrow(gap_intervals))) {
    out <- out | (t > gap_intervals$start[i] & t < gap_intervals$end[i])
  }
  out
}

zero_gaps_inhom_bg <- function(inhom_bg, gap_intervals) {
  if (nrow(gap_intervals) == 0) return(inhom_bg)
  mu_fit <- inhom_bg$mu_fit; grid <- mu_fit$grid; mu_grid <- mu_fit$mu_grid
  is_gap <- in_gap(grid, gap_intervals)
  mu_grid[is_gap] <- 0
  mu_fun_z <- approxfun(grid, mu_grid, rule = 2)
  dx <- diff(grid)
  integral_z <- sum(dx * (head(mu_grid, -1) + tail(mu_grid, -1)) / 2)
  mu_vec_z <- pmax(mu_fun_z(inhom_bg$times), 1e-12)
  mu_fit$mu_grid <- mu_grid; mu_fit$mu_fun <- mu_fun_z
  list(mu_vec = mu_vec_z, integral_bg = integral_z, times = inhom_bg$times,
       mu_fit = mu_fit, Lambda_fun = inhom_bg$Lambda_fun, gap_intervals = gap_intervals)
}

# ---- 3. Build networks ----
FORMULA_RHS <- "edges + degree(0) + gwdegree(0.1) + gwesp(0.1)"

# Day-1
df_day1 <- subset_first_session(df, gap_threshold = 1.0)
obj_day1 <- make_hypertext_net(df_day1, use_first_contact_only = TRUE)
net_day1_raw <- obj_day1$net
net_day1 <- normalize_times_01(net_day1_raw)
times_d1 <- get_times(net_day1)$times
cat("  Day-1: ", network.edgecount(net_day1), " edges, ",
    network.size(net_day1), " nodes, time [",
    round(min(times_d1), 4), ",", round(max(times_d1), 4), "]\n")
stopifnot(abs(min(times_d1)) < 1e-10, abs(max(times_d1) - 1) < 1e-10)

# Full
obj_full <- make_hypertext_net(df, use_first_contact_only = TRUE)
net_full_raw <- obj_full$net

# Gap detection on normalized scale
all_et <- sort(obj_full$edges$time)
full_rng <- range(all_et)
norm_et <- (all_et - full_rng[1]) / (full_rng[2] - full_rng[1])
gap_intervals <- find_gap_intervals(norm_et, gap_threshold = 1.0 / (full_rng[2] - full_rng[1]))

net_full <- normalize_times_01(net_full_raw)
times_f <- get_times(net_full)$times
cat("  Full:  ", network.edgecount(net_full), " edges, ",
    network.size(net_full), " nodes, time [",
    round(min(times_f), 4), ",", round(max(times_f), 4), "]\n")
cat("  Gaps found: ", nrow(gap_intervals), "\n")
stopifnot(abs(min(times_f)) < 1e-10, abs(max(times_f) - 1) < 1e-10)
stopifnot(nrow(gap_intervals) >= 2) # should have overnight gaps

# ---- 4. KDE background for full fit ----
cat("  Computing KDE background... ")
inhom_raw <- prepare_inhomogeneous_background(net_full, time_attr = "time", grid_n = 512)
inhom_full <- zero_gaps_inhom_bg(inhom_raw, gap_intervals)
cat("done. integral_bg raw=", round(inhom_raw$integral_bg, 1),
    " zeroed=", round(inhom_full$integral_bg, 1), "\n")
stopifnot(inhom_full$integral_bg < inhom_raw$integral_bg)  # zeroing must reduce

# ---- 5. Fit 1 — Day-1 (homogeneous bg, fix K only) ----
cat("\n--- Fit 1: Day-1 (maxit=3) ---\n")
exp_cs <- expected_params_PMF_mark_CS(net_day1, FORMULA_RHS)
n_cs <- exp_cs$CS_params_length
tw_d1 <- c(0, 1)
params_d1 <- list(mu = length(times_d1), beta_overall = 0.3, K = 0.5,
                  beta_edges = 0.3, node_lambda = 0.1,
                  CS_params = c(-8, -5, rep(0, n_cs - 2)))
ps_d1 <- c(mu = 1, beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.5,
            setNames(rep(0.1, n_cs), paste0("CS_params", seq_len(n_cs))))

fit1 <- fit_hawkesNet(
  params_init = params_d1, time_window = tw_d1, mark_filtration = net_day1,
  PMF_mark = PMF_mark_CS, formula_RHS = FORMULA_RHS,
  truncation = network.size(net_day1), mark_decay = "activity",
  growth_only = FALSE, max_node_time = max(get_times(net_day1)$node_times),
  method = "Nelder-Mead", maxit = 3, verbose = TRUE,
  fixed_params = "K", parscale = ps_d1, cores = 1,
  cache_intensity = TRUE, combine_intensity = TRUE
)
cat("  Fit 1 value (neg-loglik):", fit1$fit$value, "\n")
stopifnot(is.finite(fit1$fit$value))
cat("  PASS: Fit 1 completed without non-finite errors.\n")

# ---- 6. Fit 2 — Full (inhom bg, fix mu+K) ----
cat("\n--- Fit 2: Full (maxit=3) ---\n")
exp_cs2 <- expected_params_PMF_mark_CS(net_full, FORMULA_RHS)
n_cs2 <- exp_cs2$CS_params_length
tw_f <- c(0, 1)
params_f <- list(mu = inhom_full$integral_bg, beta_overall = 0.3, K = 0.5,
                 beta_edges = 0.3, node_lambda = 0.1,
                 CS_params = c(-8, -5, rep(0, n_cs2 - 2)))
ps_f <- c(beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.5,
           setNames(rep(0.1, n_cs2), paste0("CS_params", seq_len(n_cs2))))

fit2 <- fit_hawkesNet(
  params_init = params_f, time_window = tw_f, mark_filtration = net_full,
  PMF_mark = PMF_mark_CS, mu_vec = inhom_full$mu_vec,
  integral_bg = inhom_full$integral_bg,
  formula_RHS = FORMULA_RHS, truncation = network.size(net_full),
  mark_decay = "activity", growth_only = FALSE,
  max_node_time = max(get_times(net_full)$node_times),
  method = "Nelder-Mead", maxit = 3, verbose = TRUE,
  fixed_params = c("mu", "K"), parscale = ps_f, cores = 1,
  cache_intensity = TRUE, combine_intensity = TRUE
)
cat("  Fit 2 value (neg-loglik):", fit2$fit$value, "\n")
stopifnot(is.finite(fit2$fit$value))
cat("  PASS: Fit 2 completed without non-finite errors.\n")

# ---- 7. Fit 1b — Non-simple day-1 (no jitter, multi-edge events) ----
cat("\n--- Fit 1b: Non-simple day-1 (maxit=3) ---\n")

# Reload raw data WITHOUT jitter
raw_ns <- read.table(system.file("extdata", "ht09_contact_list.dat", package = "hawkesNet"))
df_ns <- data.frame(time = raw_ns$V1 / 3600, from = raw_ns$V2, to = raw_ns$V3)
df_ns$time <- df_ns$time - min(df_ns$time)
df_ns <- df_ns[order(df_ns$time), ]

df_ns_day1 <- subset_first_session(df_ns, gap_threshold = 1.0)
obj_ns <- make_hypertext_net(df_ns_day1, use_first_contact_only = TRUE)
net_ns <- normalize_times_01(obj_ns$net)
times_ns <- get_times(net_ns)$times
n_events_ns <- length(times_ns)
n_edges_ns <- network.edgecount(net_ns)
n_nodes_ns <- network.size(net_ns)
cat("  Non-simple day-1:", n_events_ns, "events |",
    n_nodes_ns, "nodes |", n_edges_ns, "edges\n")
cat("  Edges per event (mean):", round(n_edges_ns / n_events_ns, 2), "\n")

exp_cs_ns <- expected_params_PMF_mark_CS(net_ns, FORMULA_RHS)
n_cs_ns <- exp_cs_ns$CS_params_length
tw_ns <- c(0, 1)
params_ns <- list(mu = n_events_ns, beta_overall = 0.3, K = 0.5,
                  beta_edges = 0.3, node_lambda = 0.5,
                  CS_params = c(-5, -3, rep(0, n_cs_ns - 2)))
ps_ns <- c(mu = 1, beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.5,
            setNames(rep(0.1, n_cs_ns), paste0("CS_params", seq_len(n_cs_ns))))

fit1b <- fit_hawkesNet(
  params_init = params_ns, time_window = tw_ns, mark_filtration = net_ns,
  PMF_mark = PMF_mark_CS, formula_RHS = FORMULA_RHS,
  truncation = n_nodes_ns, mark_decay = "activity",
  growth_only = FALSE, max_node_time = max(get_times(net_ns)$node_times),
  method = "Nelder-Mead", maxit = 3, verbose = TRUE,
  fixed_params = "K", parscale = ps_ns, cores = 1,
  cache_intensity = TRUE, combine_intensity = TRUE
)
cat("  Fit 1b value (neg-loglik):", fit1b$fit$value, "\n")
stopifnot(is.finite(fit1b$fit$value))
cat("  PASS: Fit 1b completed without non-finite errors.\n")

# ---- 8. Simulation comparison: simple vs non-simple ----
cat("\n--- Simulation comparison ---\n")

sim_from_fit <- function(fit_obj, net_obs, inhom_bg, formula_rhs,
                         truncation, mark_decay, seed_events = 20L) {
  pfit <- fit_obj$params
  tw <- c(0, 1)
  seed_net <- NULL
  seed_times <- NULL
  all_times <- get_times(net_obs)$times
  if (length(all_times) >= seed_events) {
    t_seed <- all_times[seed_events]
    seed_net <- filtration_to_net(net_obs, t_seed, equals = TRUE)
    seed_times <- all_times[1:seed_events]
  }
  sim_hawkesNet(
    params = pfit,
    time_window = tw,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    formula_RHS = formula_rhs,
    truncation = truncation,
    mark_decay = mark_decay,
    growth_only = FALSE,
    max_node_time = Inf,
    hashed_edges = TRUE,
    verbose = FALSE,
    mu_multiplier = 5,
    stop_on_full_network = FALSE,
    inhom_bg = inhom_bg,
    seed_net = seed_net,
    seed_times = seed_times
  )
}

cat("  Simulating from Fit 1 (simple)... ")
sim_simple <- tryCatch(
  sim_from_fit(fit1, net_day1,
               list(mu_vec = NULL, integral_bg = NULL, times = times_d1),
               FORMULA_RHS, network.size(net_day1), "activity"),
  error = function(e) { cat("FAILED:", e$message, "\n"); NULL }
)
if (!is.null(sim_simple) && !is.null(sim_simple$net)) {
  cat("OK\n")
  cat("    Observed (simple):  ", network.size(net_day1), "nodes,",
      network.edgecount(net_day1), "edges,", length(times_d1), "events\n")
  cat("    Simulated (simple): ", network.size(sim_simple$net), "nodes,",
      network.edgecount(sim_simple$net), "edges,",
      length(sim_simple$events$t), "events\n")
}

cat("  Simulating from Fit 1b (non-simple)... ")
sim_nonsimple <- tryCatch(
  sim_from_fit(fit1b, net_ns,
               list(mu_vec = NULL, integral_bg = NULL, times = times_ns),
               FORMULA_RHS, n_nodes_ns, "activity"),
  error = function(e) { cat("FAILED:", e$message, "\n"); NULL }
)
if (!is.null(sim_nonsimple) && !is.null(sim_nonsimple$net)) {
  cat("OK\n")
  cat("    Observed (non-simple):  ", n_nodes_ns, "nodes,",
      n_edges_ns, "edges,", n_events_ns, "events\n")
  cat("    Simulated (non-simple): ", network.size(sim_nonsimple$net), "nodes,",
      network.edgecount(sim_nonsimple$net), "edges,",
      length(sim_nonsimple$events$t), "events\n")
}

cat("\n=== ALL CHECKS PASSED ===\n")
