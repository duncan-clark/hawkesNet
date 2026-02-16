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

cat("\n=== ALL CHECKS PASSED ===\n")
