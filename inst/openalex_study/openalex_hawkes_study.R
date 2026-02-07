# =============================================================================
# OpenAlex Hawkes study: inhomogeneous fit, temporal fit, KS tests, GOF
# =============================================================================
# Run from package root: Rscript inst/openalex_study/openalex_hawkes_study.R
# Or submit via SLURM: sbatch inst/openalex_study/run_openalex.slurm
#
# Requires: hawkesGrowthNet package (includes inhomogeneous fit and KDE background).
# =============================================================================

library(hawkesGrowthNet)
library(pbmcapply)
library(ggplot2)
library(dplyr)
library(network)
library(sna)
library(ernm)
library(parallel)

# Paths: run from package root (directory containing inst/)
PKG_ROOT <- getwd()
source(file.path(PKG_ROOT, "inst", "openalex_study", "get_network_openalex.R"))

# =============================================================================
# Config (override via env or edit)
# =============================================================================
EMAIL <- Sys.getenv("OPENALEX_EMAIL", "duncan-clark@outlook.com")
SEARCH_STRING <- Sys.getenv("OPENALEX_STRING", "Hawkes")
PAGES <- as.integer(Sys.getenv("OPENALEX_PAGES", 100))
PER_PAGE <- 100L
MIN_DATE <- "1971-04-01"
MAX_DATE <- "2020-01-01"
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 7))
MAX_ITER <- 1000
TRUNCATION <- 1000L
N_GOF <- 50L   # number of simulated networks for goodness-of-fit
PAPER_OUTPUT <- TRUE
RUN_GOF <- TRUE

# =============================================================================
# Helpers: waiting times *between* formations (triangle, 2-star, 3-star)
# =============================================================================
#' Waiting times between consecutive structure formations (triangle / 2-star / 3-star).
#'
#' Replays the network event-by-event (grouped by event time). An "event" is all
#' edges sharing the same time — typically a new node joining with its edges.
#' Records the event time whenever the triangle count, 2-star count, or 3-star
#' count increases, then returns diff() of those times.
#'
#' @param net Network with vertex/edge time attributes.
#' @param time_attr Name of the time attribute (default "time").
#' @return list(triangle = numeric(), star2 = numeric(), star3 = numeric())
waiting_times_between_formations <- function(net, time_attr = "time") {
  out <- list(triangle = numeric(), star2 = numeric(), star3 = numeric())

  # Get edge list and edge times
  el <- network::as.edgelist(net)
  if (nrow(el) == 0) return(out)

  # Try edge times first; fall back to vertex times for edge ordering
  edge_times <- NULL
  if (time_attr %in% network::list.edge.attributes(net)) {
    edge_times <- network::get.edge.attribute(net, time_attr)
  }
  if (is.null(edge_times) || length(edge_times) != nrow(el)) {
    # Use vertex times: assign each edge the max time of its endpoints
    vtimes <- network::get.vertex.attribute(net, time_attr)
    if (is.null(vtimes) || all(is.na(vtimes))) return(out)
    edge_times <- pmax(vtimes[el[, 1]], vtimes[el[, 2]], na.rm = TRUE)
  }

  # Sort edges by time
  ord <- order(edge_times)
  el <- el[ord, , drop = FALSE]
  edge_times <- edge_times[ord]

  # Group edges by unique event times
  event_times <- unique(edge_times)
  n <- network::network.size(net)
  g <- network::network.initialize(n, directed = network::is.directed(net))

  t_tri <- numeric()
  t_2s <- numeric()
  t_3s <- numeric()
  n_tri_prev <- 0
  n_2s_prev <- 0
  n_3s_prev <- 0

  for (t_cur in event_times) {
    # Add all edges for this event at once
    idx <- which(edge_times == t_cur)
    for (j in idx) {
      network::add.edges(g, tail = el[j, 1], head = el[j, 2])
    }

    # Check 2-star and 3-star counts
    degs <- sna::degree(g, gmode = "graph")
    n_2s <- sum(choose(degs, 2))
    n_3s <- sum(choose(degs, 3))
    if (n_2s > n_2s_prev) t_2s <- c(t_2s, t_cur)
    if (n_3s > n_3s_prev) t_3s <- c(t_3s, t_cur)
    n_2s_prev <- n_2s
    n_3s_prev <- n_3s

    # Triangle count
    A <- as.matrix(g, matrix.type = "adjacency")
    n_tri <- 0
    if (nrow(A) >= 3 && ncol(A) >= 3) {
      n_tri <- sum(diag(A %*% A %*% A)) / 6
    }
    if (n_tri > n_tri_prev) t_tri <- c(t_tri, t_cur)
    n_tri_prev <- n_tri
  }

  out$triangle <- if (length(t_tri) >= 2) diff(t_tri) else numeric()
  out$star2    <- if (length(t_2s) >= 2) diff(t_2s) else numeric()
  out$star3    <- if (length(t_3s) >= 2) diff(t_3s) else numeric()
  out
}

#' Degree distribution as vector of counts (degree 0, 1, 2, ... up to max_deg).
degree_dist <- function(net, max_deg = 20) {
  degs <- sna::degree(net, gmode = "graph")
  tab <- table(factor(degs, levels = 0:max_deg))
  as.vector(tab)
}

#' ESP distribution via ernm (edge-wise shared partners 0, 1, ... k).
esp_dist <- function(net, k_max = 15) {
  tryCatch({
    as.vector(ernm::calculateStatistics(net ~ esp(0:k_max)))
  }, error = function(e) rep(NA_real_, k_max + 1))
}

#' Geodesic distance distribution (upper triangle of distance matrix, excluding Inf).
geodist_dist <- function(net) {
  d <- sna::geodist(net, inf.replace = NA)
  if (is.list(d)) d <- d$gdist
  d <- as.vector(d)
  d <- d[!is.na(d) & d > 0]
  if (length(d) == 0) return(numeric(0))
  d
}

# =============================================================================
# 1. Fetch data and prepare network
# =============================================================================
t_total <- proc.time()
cat("=== OpenAlex Hawkes Study ===\n")
cat("  Search:", SEARCH_STRING, "| Pages:", PAGES, "| Cores:", N_CORES, "\n")
cat("  Date range:", MIN_DATE, "to", MAX_DATE, "\n\n")

cat("--- Step 1: Fetch data and prepare network ---\n")
t_step <- proc.time()
out <- get_network(email = EMAIL, pages = PAGES, per_page = PER_PAGE,
                   string = SEARCH_STRING, min_date = MIN_DATE, max_date = MAX_DATE)
net_raw <- out$net
edges <- out$edges
network::set.vertex.attribute(net_raw, "time", net_raw %v% "time_scaled")
network::set.edge.attribute(net_raw, "time", net_raw %e% "time_scaled")
net_raw <- hawkesGrowthNet::normalize_times_01(net_raw, attr = "time", keep_na = TRUE)
n_events <- length(hawkesGrowthNet::get_times(net_raw)$times)
n_nodes <- network::network.size(net_raw)
cat("  Network:", n_events, "events,", n_nodes, "nodes\n")
cat("  Step 1 took:", round((proc.time() - t_step)[3], 1), "s\n\n")

# =============================================================================
# 2. Inhomogeneous (KDE) + CS fit with vertex_categorical
# =============================================================================
cat("--- Step 2: Inhomogeneous (KDE) + CS fit ---\n")
t_step <- proc.time()
FORMULA_RHS <- "edges + triangles + star(c(2,3)) + nodeMix('gender')"
cat("  Formula:", FORMULA_RHS, "\n")
time_window_01 <- c(0, 1)
cat("  Preparing inhomogeneous background (KDE)...\n")
t_kde <- proc.time()
inhom_bg <- tryCatch(
  prepare_inhomogeneous_background(net_raw, time_attr = "time", bw = NULL, grid_n = 2048),
  error = function(e) { cat("  ERROR: prepare_inhomogeneous_background failed:", e$message, "\n"); NULL }
)
cat("  KDE background:", round((proc.time() - t_kde)[3], 1), "s\n")
fit_inhom <- NULL
# CS_params length must match number of change statistics from formula (ernm)
exp_cs <- expected_params_PMF_mark_CS(net_raw, FORMULA_RHS)
n_cs <- if (!is.na(exp_cs$CS_params_length)) exp_cs$CS_params_length else 5L
params_init_inhom <- list(
  mu = 1,
  beta_overall = 1,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-10, rep(0, n_cs - 1)),
  vertex_categorical = list(gender = c(female = 0.1, male = 0.5)),
  vertex_categorical_levels = list(gender = c("female", "male", "unknown"))
)
p_scale_inhom <- c(
  beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
  setNames(rep(0.1, n_cs), paste0("CS_params", seq_len(n_cs))),
  vertex_categorical.gender.female = 0.1, vertex_categorical.gender.male = 0.1
)

if (!is.null(inhom_bg)) {
  cat("  Fitting CS model (inhomogeneous + vertex_categorical)...\n")
  t_fit <- proc.time()
  fit_inhom <- tryCatch(
  fit_hawkesGrowthNet_inhom(
      params_init = params_init_inhom,
      time_window = time_window_01,
      mark_filtration = net_raw,
      PMF_mark = PMF_mark_CS,
      mu_vec = inhom_bg$mu_vec,
      integral_bg = inhom_bg$integral_bg,
      formula_RHS = FORMULA_RHS,
      truncation = TRUNCATION,
      mark_decay = "node_entrance",
      max_node_time = 1,
      method = "L-BFGS-B",
      maxit = MAX_ITER,
      trace = 1,
      reltol = 1e-8,
      verbose = FALSE,
      fixed_params = c("K", "mu"),
      parscale = p_scale_inhom,
      cache_intensity = TRUE,
      cores = N_CORES
    ),
    error = function(e) { cat("  ERROR: fit_hawkesGrowthNet_inhom failed:", e$message, "\n"); NULL }
  )
  elapsed_fit <- (proc.time() - t_fit)[3]
  if (!is.null(fit_inhom)) {
    cat("  Inhomogeneous fit completed:", round(elapsed_fit, 1), "s (",
        round(elapsed_fit / 60, 1), "min)\n")
    cat("  Fitted par:", paste(round(fit_inhom$fit$par, 4), collapse = ", "), "\n")
    cat("  Convergence code:", fit_inhom$fit$convergence, "\n")
  } else {
    cat("  Inhomogeneous fit FAILED after", round(elapsed_fit, 1), "s\n")
  }
}
cat("  Step 2 total:", round((proc.time() - t_step)[3], 1), "s\n\n")

# =============================================================================
# 3. Temporal Hawkes fit and KS test
# =============================================================================
cat("--- Step 3: Temporal Hawkes fit + KS test ---\n")
t_step <- proc.time()
t_events <- sort(unique(c(hawkesGrowthNet::get_times(net_raw)$node_times,
                         hawkesGrowthNet::get_times(net_raw)$edge_times)))
t_events <- t_events[!is.na(t_events)]
windowT <- c(min(t_events), max(t_events))
realiz <- data.frame(t = t_events)
fit_temporal <- NULL
init_gamma <- max(length(t_events) * 0.5, 10)
params_init_exp <- list(gamma = init_gamma, beta = 10, K = 0.2)
cat("  Fitting temporal Hawkes (exp kernel)...\n")
t_fit <- proc.time()
fit_temporal <- tryCatch(
  hawkesGrowthNet::fit_temporal_hawkes(
    params_init = params_init_exp,
    realiz = realiz,
    windowT = windowT,
    method = "L-BFGS-B",
    maxit = 500,
    kernel = "exp",
    trace = 0
  ),
  error = function(e) { cat("  ERROR: fit_temporal_hawkes failed:", e$message, "\n"); NULL }
)
cat("  Temporal fit:", round((proc.time() - t_fit)[3], 1), "s\n")
if (!is.null(fit_temporal)) {
  cat("  Temporal par:", paste(names(fit_temporal$par), "=", round(fit_temporal$par, 4), collapse = ", "), "\n")
}
ks_temporal_pval <- NA_real_
if (!is.null(fit_temporal) && exists("ks_test_pval_temporal")) {
  cat("  Computing KS test...\n")
  ks_temporal_pval <- tryCatch(
    hawkesGrowthNet::ks_test_pval_temporal(
      realiz = realiz,
      windowT = windowT,
      hawkes_par = fit_temporal$par,
      kernel = "exp",
      use_kde = TRUE
    ),
    error = function(e) NA_real_
  )
  cat("  Temporal KS p-value:", ks_temporal_pval, "\n")
}
cat("  Step 3 total:", round((proc.time() - t_step)[3], 1), "s\n\n")

# =============================================================================
# 4. Goodness-of-fit: simulate from fitted model, compare degree/ESP/geodesic/waiting times
# =============================================================================
cat("--- Step 4: Goodness-of-fit ---\n")
t_step <- proc.time()
GOF_results <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                    geodist_obs = NULL, geodist_sim = NULL,
                    wait_obs = NULL, wait_sim = NULL)
if (RUN_GOF && !is.null(fit_inhom)) {
  cat("  Simulating", N_GOF, "networks from fitted model...\n")
  # Reconstruct fitted params: strip vertex_categorical_levels from skeleton
  # (the fitter stripped it before unlist, so fit$par doesn't include it)
  skel <- params_init_inhom
  skel$vertex_categorical_levels <- NULL
  pfit <- relist(fit_inhom$fit$par, skeleton = skel)
  # Restore metadata and fixed params
  pfit$vertex_categorical_levels <- params_init_inhom$vertex_categorical_levels
  pfit$K <- params_init_inhom$K
  pfit$mu <- inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1])
  Tval <- time_window_01[2] - time_window_01[1]
  sim_nets <- list()
  
  # fix soome params just for now:
  pfit$K <- min(max(pfit$K, 0.001), 0.999)  # K must be in (0,1) for stability
  pfit$mu <- max(pfit$mu, 0.001)  # mu must be positive
  pfit$node_lambda <- max(pfit$node_lambda, 1)  # node_lambda must be >= 1 for stability
  pfit$beta_overall <- max(pfit$beta_overall, 0.001)
  pfit$beta_edges <- max(pfit$beta_edges, 0.001)
  
  pfit$vertex_categorical$gender <- c(0.1,0.5)
  
  n_success <- 0L
  n_fail <- 0L
  for (i in seq_len(N_GOF)) {
    t_sim_i <- proc.time()
    s <- tryCatch(
      sim_hawkesGrowthNet(
        params = pfit,
        time_window = c(0,0.05),
        PMF_mark = PMF_mark_CS,
        cond_intensity = cond_intensity,
        formula_RHS = FORMULA_RHS,
        truncation = TRUNCATION,
        mark_decay = "activity",
        max_node_time = 1,
        hashed_edges = TRUE,
        verbose = FALSE,
        mu_multiplier = 5,
        stop_on_full_network = FALSE
      ),
      error = function(e) { cat("  GOF sim", i, "FAILED:", e$message, "\n"); NULL }
    )
    elapsed_i <- round((proc.time() - t_sim_i)[3], 1)
    if (!is.null(s)) {
      n_success <- n_success + 1L
      sim_nets[[i]] <- s$net
      n_sim_nodes <- network::network.size(s$net)
      n_sim_edges <- network::network.edgecount(s$net)
      cat("  GOF sim", i, "/", N_GOF, ":", n_sim_nodes, "nodes,",
          n_sim_edges, "edges (", elapsed_i, "s)\n")
    } else {
      n_fail <- n_fail + 1L
    }
  }
  sim_nets <- sim_nets[!sapply(sim_nets, is.null)]
  cat("  GOF simulations:", n_success, "succeeded,", n_fail, "failed\n")
  
  if (length(sim_nets) > 0) {
    cat("  Computing GOF statistics...\n")
    t_stats <- proc.time()
    max_deg <- 15
    k_esp <- 15
    GOF_results$degree_obs <- degree_dist(net_raw, max_deg)
    GOF_results$degree_sim <- do.call(rbind, lapply(sim_nets, function(n) degree_dist(n, max_deg)))
    cat("    Degree: done\n")
    GOF_results$esp_obs <- esp_dist(net_raw, k_esp)
    GOF_results$esp_sim <- do.call(rbind, lapply(sim_nets, function(n) esp_dist(n, k_esp)))
    cat("    ESP: done\n")
    GOF_results$geodist_obs <- geodist_dist(net_raw)
    GOF_results$geodist_sim <- lapply(sim_nets, geodist_dist)
    cat("    Geodesic: done\n")
    GOF_results$wait_obs <- waiting_times_between_formations(net_raw)
    GOF_results$wait_sim <- lapply(sim_nets, waiting_times_between_formations)
    cat("    Waiting times: done\n")
    cat("  GOF statistics:", round((proc.time() - t_stats)[3], 1), "s\n")
  }
} else {
  if (!RUN_GOF) cat("  RUN_GOF = FALSE; skipping\n")
  if (is.null(fit_inhom)) cat("  No fit available; skipping GOF\n")
}
cat("  Step 4 total:", round((proc.time() - t_step)[3], 1), "s\n\n")

# =============================================================================
# 5. Save full state for rehydration
# =============================================================================
cat("--- Step 5: Save full state ---\n")
save_list <- list(
  net_raw = net_raw,
  edges = edges,
  inhom_bg = inhom_bg,
  fit_inhom = fit_inhom,
  params_init_inhom = params_init_inhom,
  fit_temporal = fit_temporal,
  ks_temporal_pval = ks_temporal_pval,
  realiz = realiz,
  windowT = windowT,
  GOF_results = GOF_results,
  N_GOF = N_GOF,
  SEARCH_STRING = SEARCH_STRING,
  time_window_01 = time_window_01
)
saveRDS(save_list, file.path(PKG_ROOT, "inst", "openalex_study", "results_openalex_full.RDS"))
cat("  Saved to inst/openalex_study/results_openalex_full.RDS\n\n")

# =============================================================================
# 6. PAPER_OUTPUT: rehydrate and produce figures/tables
# =============================================================================
if (PAPER_OUTPUT) {
  cat("--- Step 6: Paper output (figures & tables) ---\n")
  t_step <- proc.time()
  dat <- readRDS(file.path(PKG_ROOT, "inst", "openalex_study", "results_openalex_full.RDS"))
  list2env(dat, envir = .GlobalEnv)
  cat("  Rehydrated; producing figures and tables.\n")

  if (!is.null(dat$fit_inhom)) {
    if (!is.null(dat$fit_inhom$fit_table)) {
      cat("  Inhomogeneous fit: parameter estimates and standard errors\n")
      print(dat$fit_inhom$fit_table)
    } else {
      print(fit_inhom$fit$par)
    }
    if (exists("params_init_inhom") && !is.null(params_init_inhom$vertex_categorical)) {
      skel2 <- params_init_inhom
      skel2$vertex_categorical_levels <- NULL
      pfit2 <- relist(fit_inhom$fit$par, skeleton = skel2)
      pfit2$vertex_categorical_levels <- params_init_inhom$vertex_categorical_levels
      if (!is.null(pfit2$vertex_categorical$gender)) {
        levs <- params_init_inhom$vertex_categorical_levels$gender
        pgender <- expand_vertex_categorical_probs(pfit2$vertex_categorical$gender, levs)
        cat("  Fitted gender proportions (n-1 expanded):\n"); print(pgender)
      }
    }
  }
  if (!is.null(dat$fit_temporal)) {
    cat("  Temporal Hawkes par:\n"); print(fit_temporal$par)
    cat("  Temporal KS p-value:", ks_temporal_pval, "\n")
  }

  # GOF plots: degree, ESP, geodesic, waiting times
  if (!is.null(dat$GOF_results) && !is.null(dat$GOF_results$degree_obs)) {
    gof <- dat$GOF_results
    max_deg <- length(gof$degree_obs) - 1
    # Degree: observed vs simulated boxplots
    deg_df <- rbind(
      data.frame(degree = 0:max_deg, count = gof$degree_obs, type = "Observed"),
      data.frame(degree = rep(0:max_deg, each = nrow(gof$degree_sim)),
                 count = as.vector(gof$degree_sim),
                 type = "Simulated")
    )
    p_deg <- ggplot(deg_df, aes(x = factor(degree), y = count, fill = type)) +
      geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
      labs(title = "GOF: Degree distribution", x = "Degree", y = "Count") +
      theme_minimal() + theme(legend.position = "bottom")
    print(p_deg)
    # ESP
    esp_obs <- gof$esp_obs
    esp_sim <- gof$esp_sim
    if (!is.null(esp_sim) && nrow(esp_sim) > 0) {
      esp_df <- rbind(
        data.frame(esp = 0:(length(esp_obs)-1), value = esp_obs, type = "Observed"),
        data.frame(esp = rep(0:(ncol(esp_sim)-1), each = nrow(esp_sim)),
                   value = as.vector(esp_sim), type = "Simulated")
      )
      p_esp <- ggplot(esp_df, aes(x = factor(esp), y = value, fill = type)) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
        labs(title = "GOF: ESP distribution", x = "ESP", y = "Count") +
        theme_minimal() + theme(legend.position = "bottom")
      print(p_esp)
    }
    # Geodesic: boxplot of pair counts at each distance + ECDF
    g_obs <- gof$geodist_obs
    g_sim <- gof$geodist_sim
    if (length(g_obs) > 0 && length(g_sim) > 0) {
      # Boxplot: tabulate counts at each integer distance
      max_geod <- min(max(c(g_obs, unlist(g_sim)), na.rm = TRUE), 20)
      geod_levels <- seq_len(max_geod)
      obs_tab <- table(factor(g_obs, levels = geod_levels))
      sim_tabs <- lapply(g_sim, function(g) {
        as.vector(table(factor(g, levels = geod_levels)))
      })
      sim_mat <- do.call(rbind, sim_tabs)
      geod_box_df <- rbind(
        data.frame(distance = geod_levels, count = as.vector(obs_tab), type = "Observed"),
        data.frame(distance = rep(geod_levels, each = nrow(sim_mat)),
                   count = as.vector(sim_mat), type = "Simulated")
      )
      p_geod_box <- ggplot(geod_box_df, aes(x = factor(distance), y = count, fill = type)) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
        labs(title = "GOF: Geodesic distance distribution",
             x = "Geodesic distance", y = "Number of pairs") +
        theme_minimal() + theme(legend.position = "bottom")
      print(p_geod_box)

      # ECDF version
      max_d <- max(c(g_obs, unlist(g_sim)), na.rm = TRUE)
      x_seq <- seq(0, min(max_d, 20), length.out = 200)
      ecdf_obs <- sapply(x_seq, function(x) mean(g_obs <= x, na.rm = TRUE))
      ecdf_sim <- sapply(x_seq, function(x) mean(unlist(g_sim) <= x, na.rm = TRUE))
      geod_df <- rbind(
        data.frame(dist = x_seq, ecdf = ecdf_obs, type = "Observed"),
        data.frame(dist = x_seq, ecdf = ecdf_sim, type = "Simulated")
      )
      p_geod <- ggplot(geod_df, aes(x = dist, y = ecdf, color = type)) +
        geom_line(linewidth = 1) +
        labs(title = "GOF: Geodesic distance (ECDF)", x = "Distance", y = "ECDF") +
        theme_minimal() + theme(legend.position = "bottom")
      print(p_geod)
    }
    # Waiting times between formations: triangle, 2-star, 3-star
    w_obs <- gof$wait_obs
    w_sim <- gof$wait_sim
    if (!is.null(w_obs) && is.list(w_obs) && !is.null(w_sim) && length(w_sim) > 0) {
      obs_vec <- c(w_obs$triangle, w_obs$star2, w_obs$star3)
      obs_metric <- rep(c("Triangle", "2-star", "3-star"),
                       c(length(w_obs$triangle), length(w_obs$star2), length(w_obs$star3)))
      sim_vec <- unlist(lapply(w_sim, function(x) c(x$triangle, x$star2, x$star3)))
      sim_metric <- unlist(lapply(w_sim, function(x) rep(c("Triangle", "2-star", "3-star"),
                         c(length(x$triangle), length(x$star2), length(x$star3)))))
      wait_df <- rbind(
        data.frame(metric = obs_metric, value = obs_vec, type = "Observed"),
        data.frame(metric = sim_metric, value = sim_vec, type = "Simulated")
      )
      wait_df <- wait_df[!is.na(wait_df$value), ]
      if (nrow(wait_df) > 0) {
        p_wait <- ggplot(wait_df, aes(x = metric, y = value, fill = type)) +
          geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
          labs(title = "GOF: Waiting time between formations (triangle / 2-star / 3-star)",
               x = "", y = "Waiting time") +
          theme_minimal() + theme(legend.position = "bottom")
        print(p_wait)
      }
    }
  }
  cat("  Step 6 total:", round((proc.time() - t_step)[3], 1), "s\n\n")
}

# =============================================================================
# Total elapsed time
# =============================================================================
cat("=== OpenAlex study complete ===\n")
cat("  Total wall time:", round((proc.time() - t_total)[3], 1), "s (",
    round((proc.time() - t_total)[3] / 60, 1), "min)\n")
