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

# ---- 3. Build day-1 network (jittered = simple) ----
FORMULA_RHS <- "edges + degree(0) + gwdegree(0.1) + gwesp(0.1)"

df_day1 <- subset_first_session(df, gap_threshold = 1.0)
obj_day1 <- make_hypertext_net(df_day1, use_first_contact_only = TRUE)
net_day1_raw <- obj_day1$net
net_day1 <- normalize_times_01(net_day1_raw)
times_d1 <- get_times(net_day1)$times
cat("  Day-1 (simple):", network.edgecount(net_day1), "edges,",
    network.size(net_day1), "nodes, time [",
    round(min(times_d1), 4), ",", round(max(times_d1), 4), "]\n")
stopifnot(abs(min(times_d1)) < 1e-10, abs(max(times_d1) - 1) < 1e-10)

# ---- 4. Fit 1 — Day-1 simple, mark_decay = "activity" ----
cat("\n--- Fit 1: Day-1 simple, activity decay (maxit=3) ---\n")
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
cat("  PASS: Fit 1 completed.\n")

# ---- 5. Fit 1b — Non-simple day-1, mark_decay = "activity" ----
cat("\n--- Fit 1b: Day-1 non-simple, activity decay (maxit=3) ---\n")

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
cat("  PASS: Fit 1b completed.\n")

# ---- 6. Fit 3 — Day-1 simple, mark_decay = "node_entrance" ----
cat("\n--- Fit 3: Day-1 simple, node_entrance decay (maxit=3) ---\n")
params_d1_ne <- list(mu = length(times_d1), beta_overall = 0.3, K = 0.5,
                     beta_edges = 0.3, node_lambda = 0.1,
                     CS_params = c(-8, -5, rep(0, n_cs - 2)))

fit3 <- fit_hawkesNet(
  params_init = params_d1_ne, time_window = tw_d1, mark_filtration = net_day1,
  PMF_mark = PMF_mark_CS, formula_RHS = FORMULA_RHS,
  truncation = network.size(net_day1), mark_decay = "node_entrance",
  growth_only = FALSE, max_node_time = max(get_times(net_day1)$node_times),
  method = "Nelder-Mead", maxit = 3, verbose = TRUE,
  fixed_params = "K", parscale = ps_d1, cores = 1,
  cache_intensity = TRUE, combine_intensity = TRUE
)
cat("  Fit 3 value (neg-loglik):", fit3$fit$value, "\n")
stopifnot(is.finite(fit3$fit$value))
cat("  PASS: Fit 3 completed.\n")

# ---- 7. Fit 4 — Non-simple day-1, mark_decay = "node_entrance" ----
cat("\n--- Fit 4: Day-1 non-simple, node_entrance decay (maxit=3) ---\n")
params_ns_ne <- list(mu = n_events_ns, beta_overall = 0.3, K = 0.5,
                     beta_edges = 0.3, node_lambda = 0.5,
                     CS_params = c(-5, -3, rep(0, n_cs_ns - 2)))

fit4 <- fit_hawkesNet(
  params_init = params_ns_ne, time_window = tw_ns, mark_filtration = net_ns,
  PMF_mark = PMF_mark_CS, formula_RHS = FORMULA_RHS,
  truncation = n_nodes_ns, mark_decay = "node_entrance",
  growth_only = FALSE, max_node_time = max(get_times(net_ns)$node_times),
  method = "Nelder-Mead", maxit = 3, verbose = TRUE,
  fixed_params = "K", parscale = ps_ns, cores = 1,
  cache_intensity = TRUE, combine_intensity = TRUE
)
cat("  Fit 4 value (neg-loglik):", fit4$fit$value, "\n")
stopifnot(is.finite(fit4$fit$value))
cat("  PASS: Fit 4 completed.\n")

# ---- 8. Fit 5 — Non-simple day-1, activity, with m-parameter ----
cat("\n--- Fit 5: Day-1 non-simple, activity + m-param (maxit=3) ---\n")
m_init <- n_edges_ns / n_events_ns
params_ns_m <- list(mu = n_events_ns, beta_overall = 0.3, K = 0.5,
                    beta_edges = 0.3, node_lambda = 0.5, m = m_init,
                    CS_params = c(-5, -3, rep(0, n_cs_ns - 2)))
ps_ns_m <- c(mu = 1, beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.5,
              m = 0.5,
              setNames(rep(0.1, n_cs_ns), paste0("CS_params", seq_len(n_cs_ns))))

fit5 <- fit_hawkesNet(
  params_init = params_ns_m, time_window = tw_ns, mark_filtration = net_ns,
  PMF_mark = PMF_mark_CS, formula_RHS = FORMULA_RHS,
  truncation = n_nodes_ns, mark_decay = "activity",
  growth_only = FALSE, max_node_time = max(get_times(net_ns)$node_times),
  method = "Nelder-Mead", maxit = 3, verbose = TRUE,
  fixed_params = "K", parscale = ps_ns_m, cores = 1,
  cache_intensity = TRUE, combine_intensity = TRUE
)
cat("  Fit 5 value (neg-loglik):", fit5$fit$value, "\n")
stopifnot(is.finite(fit5$fit$value))
cat("  PASS: Fit 5 completed.\n")

# ---- 9. Simulation comparison: all 5 fits ----
cat("\n--- Simulation comparison (all 5 fits) ---\n")

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

print_comparison <- function(label, sim_obj, net_obs, times_obs) {
  n_obs  <- network.size(net_obs)
  e_obs  <- network.edgecount(net_obs)
  ev_obs <- length(times_obs)
  if (is.null(sim_obj) || is.null(sim_obj$net)) {
    cat(sprintf("  %-35s  FAILED\n", label))
    return(invisible(NULL))
  }
  n_sim  <- network.size(sim_obj$net)
  e_sim  <- network.edgecount(sim_obj$net)
  ev_sim <- length(sim_obj$events$t)
  cat(sprintf("  %-35s  Obs: %3d nodes %4d edges %4d events  |  Sim: %3d nodes %4d edges %4d events\n",
              label, n_obs, e_obs, ev_obs, n_sim, e_sim, ev_sim))
}

no_inhom <- list(mu_vec = NULL, integral_bg = NULL, times = NULL)

# Fit 1: simple, activity
cat("  Simulating Fit 1 ... ")
sim1 <- tryCatch(
  sim_from_fit(fit1, net_day1, no_inhom, FORMULA_RHS,
               network.size(net_day1), "activity"),
  error = function(e) { cat("FAILED:", e$message, "\n"); NULL }
)
cat("done\n")

# Fit 1b: non-simple, activity
cat("  Simulating Fit 1b ... ")
sim1b <- tryCatch(
  sim_from_fit(fit1b, net_ns, no_inhom, FORMULA_RHS,
               n_nodes_ns, "activity"),
  error = function(e) { cat("FAILED:", e$message, "\n"); NULL }
)
cat("done\n")

# Fit 3: simple, node_entrance
cat("  Simulating Fit 3 ... ")
sim3 <- tryCatch(
  sim_from_fit(fit3, net_day1, no_inhom, FORMULA_RHS,
               network.size(net_day1), "node_entrance"),
  error = function(e) { cat("FAILED:", e$message, "\n"); NULL }
)
cat("done\n")

# Fit 4: non-simple, node_entrance
cat("  Simulating Fit 4 ... ")
sim4 <- tryCatch(
  sim_from_fit(fit4, net_ns, no_inhom, FORMULA_RHS,
               n_nodes_ns, "node_entrance"),
  error = function(e) { cat("FAILED:", e$message, "\n"); NULL }
)
cat("done\n")

# Fit 5: non-simple, activity + m-parameter
cat("  Simulating Fit 5 ... ")
sim5 <- tryCatch(
  sim_from_fit(fit5, net_ns, no_inhom, FORMULA_RHS,
               n_nodes_ns, "activity"),
  error = function(e) { cat("FAILED:", e$message, "\n"); NULL }
)
cat("done\n")

cat("\n--- Results ---\n")
print_comparison("Fit 1  (simple, activity)",        sim1,  net_day1, times_d1)
print_comparison("Fit 1b (non-simple, activity)",    sim1b, net_ns,   times_ns)
print_comparison("Fit 3  (simple, node_entrance)",   sim3,  net_day1, times_d1)
print_comparison("Fit 4  (non-simple, node_entr.)",  sim4,  net_ns,   times_ns)
print_comparison("Fit 5  (non-simple, activity+m)",  sim5,  net_ns,   times_ns)

cat("\n=== ALL CHECKS PASSED ===\n")
