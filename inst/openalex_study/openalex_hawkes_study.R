# =============================================================================
# OpenAlex Hawkes study: inhomogeneous fit, temporal fit, KS tests, GOF
# =============================================================================
# Run from package root: Rscript inst/openalex_study/openalex_hawkes_study.R
# Or submit via SLURM: sbatch inst/openalex_study/run_openalex.slurm
#
# Requires: hawkesGrowthNet package (includes inhomogeneous fit and KDE background).
# =============================================================================

library(hawkesGrowthNet)
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
PAGES <- as.integer(Sys.getenv("OPENALEX_PAGES", 50))
PER_PAGE <- 100L
MIN_DATE <- "1971-04-01"
MAX_DATE <- "2020-01-01"
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 7))
MAX_ITER <- 2000L
TRUNCATION <- 100L
N_GOF <- 50L   # number of simulated networks for goodness-of-fit
PAPER_OUTPUT <- TRUE
RUN_GOF <- TRUE

# =============================================================================
# Helpers: waiting times *between* formations (triangle, 2-star, 3-star)
# =============================================================================
#' From a network with edge attribute "time", return waiting times between
#' consecutive formations: each time the triangle count (or 2-star / 3-star count)
#' increases, record that time; return diff() of those times.
#' Returns list(triangle = numeric(), star2 = numeric(), star3 = numeric()).
waiting_times_between_formations <- function(net, time_attr = "time") {
  out <- list(triangle = numeric(), star2 = numeric(), star3 = numeric())
  if (!time_attr %in% network::list.edge.attributes(net)) return(out)
  el <- network::as.edgelist(net)
  times <- network::get.edge.attribute(net, time_attr)
  if (length(times) != nrow(el)) return(out)
  ord <- order(times)
  el <- el[ord, , drop = FALSE]
  times <- times[ord]
  n <- network::network.size(net)
  g <- network::network.initialize(n, directed = FALSE)
  t_tri <- numeric()
  t_2s <- numeric()
  t_3s <- numeric()
  n_tri_prev <- 0
  n_2s_prev <- 0
  n_3s_prev <- 0
  for (i in seq_len(nrow(el))) {
    network::add.edges(g, tail = el[i, 1], head = el[i, 2])
    t_cur <- times[i]
    degs <- sna::degree(g, gmode = "graph")
    n_2s <- sum(degs >= 2)
    n_3s <- sum(degs >= 3)
    if (n_2s > n_2s_prev) t_2s <- c(t_2s, t_cur)
    if (n_3s > n_3s_prev) t_3s <- c(t_3s, t_cur)
    n_2s_prev <- n_2s
    n_3s_prev <- n_3s
    if (nrow(el) >= 3) {
      A <- as.matrix(g)
      n_tri <- 0
      if (nrow(A) >= 3 && ncol(A) >= 3) {
        n_tri <- sum(diag(A %*% A %*% A)) / 6
      }
      if (n_tri > n_tri_prev) t_tri <- c(t_tri, t_cur)
      n_tri_prev <- n_tri
    }
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
message("Fetching OpenAlex data: ", SEARCH_STRING, " ...")
out <- get_network(email = EMAIL, pages = PAGES, per_page = PER_PAGE,
                   string = SEARCH_STRING, min_date = MIN_DATE, max_date = MAX_DATE)
net_raw <- out$net
edges <- out$edges
network::set.vertex.attribute(net_raw, "time", net_raw %v% "time_scaled")
network::set.edge.attribute(net_raw, "time", net_raw %e% "time_scaled")
net_raw <- hawkesGrowthNet::normalize_times_01(net_raw, attr = "time", keep_na = TRUE)
n_events <- length(hawkesGrowthNet::get_times(net_raw)$times)
n_nodes <- network::network.size(net_raw)
message(sprintf("Network: %d events, %d nodes.", n_events, n_nodes))

# =============================================================================
# 2. Inhomogeneous (KDE) + CS fit with vertex_categorical
# =============================================================================
time_window_01 <- c(0, 1)
inhom_bg <- tryCatch(
  hawkesGrowthNet::prepare_inhomogeneous_background(net_raw, time_attr = "time", bw = NULL, grid_n = 2048),
  error = function(e) { message("prepare_inhomogeneous_background failed: ", e$message); NULL }
)
fit_inhom <- NULL
params_init_inhom <- list(
  mu = 1,
  beta_overall = 1,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-10, rep(0, 9)),
  vertex_categorical = list(gender = c(male = 0.33, female = 0.33, unknown = 0.34))
)
p_scale_inhom <- c(
  beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
  CS_params1 = 1, CS_params2 = 0.1, CS_params3 = 0.1, CS_params4 = 0.01,
  CS_params5 = 0.1, CS_params6 = 0.1, CS_params7 = 0.1, CS_params8 = 0.1,
  CS_params9 = 0.1, CS_params10 = 0.1,
  vertex_categorical.gender.male = 0.1, vertex_categorical.gender.female = 0.1,
  vertex_categorical.gender.unknown = 0.1
)
if (!is.null(inhom_bg)) {
  message("Fitting CS model with inhomogeneous (KDE) background and vertex_categorical (gender)...")
  fit_inhom <- tryCatch(
    hawkesGrowthNet::fit_hawkesGrowthNet_inhom(
      params_init = params_init_inhom,
      time_window = time_window_01,
      mark_filtration = net_raw,
      PMF_mark = hawkesGrowthNet::PMF_mark_CS,
      mu_vec = inhom_bg$mu_vec,
      integral_bg = inhom_bg$integral_bg,
      formula_RHS = "edges + triangles + star(c(2,3)) + nodeMix('gender')",
      truncation = TRUNCATION,
      mark_decay = "activity",
      max_node_time = 1,
      maxit = MAX_ITER,
      trace = 1,
      reltol = 1e-8,
      verbose = FALSE,
      get_hessian = TRUE,
      fixed_params = c("K", "mu"),
      parscale = p_scale_inhom,
      cache_intensity = TRUE,
      cores = N_CORES
    ),
    error = function(e) { message("fit_hawkesGrowthNet_inhom failed: ", e$message); NULL }
  )
  if (!is.null(fit_inhom)) message("Inhomogeneous fit par: ", paste(round(fit_inhom$fit$par, 4), collapse = ", "))
}

# =============================================================================
# 3. Temporal Hawkes fit and KS test
# =============================================================================
t_events <- sort(unique(c(hawkesGrowthNet::get_times(net_raw)$node_times,
                         hawkesGrowthNet::get_times(net_raw)$edge_times)))
t_events <- t_events[!is.na(t_events)]
windowT <- c(min(t_events), max(t_events))
realiz <- data.frame(t = t_events)
fit_temporal <- NULL
init_gamma <- max(length(t_events) * 0.5, 10)
params_init_exp <- list(gamma = init_gamma, beta = 10, K = 0.2)
message("Fitting temporal Hawkes (exp kernel, KDE background)...")
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
  error = function(e) { message("fit_temporal_hawkes failed: ", e$message); NULL }
)
ks_temporal_pval <- NA_real_
if (!is.null(fit_temporal) && exists("ks_test_pval_temporal")) {
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
  message(sprintf("Temporal KS test p-value: %s", ks_temporal_pval))
}

# =============================================================================
# 4. Goodness-of-fit: simulate from fitted model, compare degree/ESP/geodesic/waiting times
# =============================================================================
GOF_results <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                    geodist_obs = NULL, geodist_sim = NULL,
                    wait_obs = NULL, wait_sim = NULL)
if (RUN_GOF && !is.null(fit_inhom)) {
  message("GOF: simulating ", N_GOF, " networks from fitted model...")
  pfit <- relist(fit_inhom$fit$par, skeleton = params_init_inhom)
  pfit$K <- 0.5
  pfit$mu <- inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1])
  Tval <- time_window_01[2] - time_window_01[1]
  sim_nets <- list()
  for (i in seq_len(N_GOF)) {
    s <- tryCatch(
      hawkesGrowthNet::sim_hawkesGrowthNet(
        params = pfit,
        time_window = time_window_01,
        PMF_mark = hawkesGrowthNet::PMF_mark_CS,
        cond_intensity = hawkesGrowthNet::cond_intensity,
        formula_RHS = "edges + triangles + star(c(2,3)) + nodeMix('gender')",
        truncation = TRUNCATION,
        hashed_edges = TRUE,
        verbose = FALSE,
        mu_multiplier = 5
      ),
      error = function(e) NULL
    )
    if (!is.null(s)) sim_nets[[i]] <- s$net
  }
  sim_nets <- sim_nets[!sapply(sim_nets, is.null)]
  message("GOF: ", length(sim_nets), " simulated networks.")
  if (length(sim_nets) > 0) {
    max_deg <- 15
    k_esp <- 15
    GOF_results$degree_obs <- degree_dist(net_raw, max_deg)
    GOF_results$degree_sim <- do.call(rbind, lapply(sim_nets, function(n) degree_dist(n, max_deg)))
    GOF_results$esp_obs <- esp_dist(net_raw, k_esp)
    GOF_results$esp_sim <- do.call(rbind, lapply(sim_nets, function(n) esp_dist(n, k_esp)))
    GOF_results$geodist_obs <- geodist_dist(net_raw)
    GOF_results$geodist_sim <- lapply(sim_nets, geodist_dist)
    GOF_results$wait_obs <- waiting_times_between_formations(net_raw)
    GOF_results$wait_sim <- lapply(sim_nets, waiting_times_between_formations)
  }
}

# =============================================================================
# 5. Save full state for rehydration
# =============================================================================
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
message("Saved full state to inst/openalex_study/results_openalex_full.RDS")

# =============================================================================
# 6. PAPER_OUTPUT: rehydrate and produce figures/tables
# =============================================================================
if (PAPER_OUTPUT) {
  dat <- readRDS(file.path(PKG_ROOT, "inst", "openalex_study", "results_openalex_full.RDS"))
  list2env(dat, envir = .GlobalEnv)
  message("Rehydrated; producing figures and tables.")

  if (!is.null(dat$fit_inhom)) {
    print(fit_inhom$fit$par)
    if (exists("params_init_inhom") && !is.null(params_init_inhom$vertex_categorical)) {
      pfit <- relist(fit_inhom$fit$par, skeleton = params_init_inhom)
      if (!is.null(pfit$vertex_categorical$gender)) {
        pgender <- pfit$vertex_categorical$gender / sum(pfit$vertex_categorical$gender)
        message("Fitted gender proportions: "); print(pgender)
      }
    }
  }
  if (!is.null(dat$fit_temporal)) {
    message("Temporal Hawkes par: "); print(fit_temporal$par)
    message("Temporal KS p-value: ", ks_temporal_pval)
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
    # Geodesic: ECDF observed vs simulated
    g_obs <- gof$geodist_obs
    g_sim <- gof$geodist_sim
    if (length(g_obs) > 0 && length(g_sim) > 0) {
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
  message("PAPER_OUTPUT done.")
}
