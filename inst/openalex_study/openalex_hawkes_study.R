# =============================================================================
# OpenAlex Hawkes study: inhomogeneous fit, temporal fit, KS tests, GOF
# =============================================================================
# Run from package root: Rscript inst/openalex_study/openalex_hawkes_study.R
# Or submit via SLURM: sbatch inst/openalex_study/run_openalex.slurm
#
# Requires: hawkesGrowthNet package (includes inhomogeneous fit and KDE background).
# =============================================================================

# Load hawkesGrowthNet package
library(hawkesGrowthNet)

# Load required libraries
library(dplyr)
library(network)
library(sna)
library(ernm)
library(parallel)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
}

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
MAX_ITER <- 5000
TRUNCATION <- 100L
N_GOF <- 50L   # number of simulated networks for goodness-of-fit
PAPER_OUTPUT <- TRUE
RUN_GOF <- TRUE

# =============================================================================
# GOF functions are now in R/gof.R (exported from package)
# =============================================================================

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
# Check gender distribution
if ("gender" %in% network::list.vertex.attributes(net_raw)) {
  gender_vals <- net_raw %v% "gender"
  gender_counts <- table(gender_vals, useNA = "ifany")
  cat("  Gender distribution:", paste(names(gender_counts), "=", gender_counts, collapse = ", "), "\n")
  if (all(gender_vals == "unknown", na.rm = TRUE)) {
    warning("⚠ WARNING: All genders are 'unknown'. This will cause nodeMix to only detect 1 statistic instead of 6.")
    warning("⚠ Install 'gender' and 'genderdata' packages for gender prediction: install.packages(c('gender', 'genderdata'))")
  }
} else {
  warning("⚠ No 'gender' attribute found on network vertices")
}
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
  cat("  Method: Nelder-Mead (max 500 iterations)\n")
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
      mark_decay = "activity",
      max_node_time = 1,
      method = "Nelder-Mead",
      maxit = MAX_ITER,
      trace = 1,
      reltol = 1e-8,
      verbose = FALSE,
      fixed_params = c("K", "mu"),
      parscale = p_scale_inhom,
      cache_intensity = TRUE,
      cores = N_CORES
    ),
    error = function(e) { cat("  ERROR: Fit failed:", e$message, "\n"); NULL }
  )
  
  elapsed_fit <- (proc.time() - t_fit)[3]
  if (!is.null(fit_inhom)) {
    cat("  Fit completed:", round(elapsed_fit, 1), "s (", round(elapsed_fit / 60, 1), "min)\n")
    cat("  Convergence:", fit_inhom$fit$convergence, "\n")
    cat("  Iterations:", fit_inhom$fit$counts[1], "\n")
  } else {
    cat("  Fit FAILED after", round(elapsed_fit, 1), "s\n")
  }
}

# =============================================================================
# 2b. Inhomogeneous fit with nodeMatch (homophily) instead of nodeMix
# =============================================================================
fit_inhom_nodematch <- NULL
FORMULA_RHS_NODEMATCH <- "edges + triangles + star(c(2,3)) + nodeMatch('gender')"
if (!is.null(inhom_bg) && !is.null(fit_inhom)) {
  cat("--- Step 2b: Inhomogeneous fit with nodeMatch (homophily) ---\n")
  cat("  Formula:", FORMULA_RHS_NODEMATCH, "\n")
  cat("  Using nodeMix fit as starting point...\n")
  t_step_nodematch <- proc.time()
  
  # Get expected parameters for nodeMatch formula
  exp_cs_nodematch <- expected_params_PMF_mark_CS(net_raw, FORMULA_RHS_NODEMATCH)
  n_cs_nodematch <- if (!is.na(exp_cs_nodematch$CS_params_length)) exp_cs_nodematch$CS_params_length else 4L
  
  # Extract fitted parameters from nodeMix fit
  skel <- params_init_inhom
  skel$vertex_categorical_levels <- NULL
  params_fitted <- relist(fit_inhom$fit$par, skeleton = skel)
  params_fitted$vertex_categorical_levels <- params_init_inhom$vertex_categorical_levels
  params_fitted <- hawkesGrowthNet:::reconstruct_vertex_categorical_names(
    params_fitted, params_init_inhom$vertex_categorical_levels)
  params_fitted <- hawkesGrowthNet:::repair_vertex_categorical_params(params_fitted, eps = 1e-6)
  
  # Build nodeMatch params_init: keep structural params, replace nodeMix with nodeMatch
  # Structural CS_params (edges, triangles, stars) are first 4; nodeMix terms are 5-10
  cs_structural <- head(params_fitted$CS_params, 4)  # edges, triangles, star2, star3
  cs_nodemix <- tail(params_fitted$CS_params, -4)     # nodeMix terms (if any)
  cs_nodematch <- if (length(cs_nodemix) > 0) mean(cs_nodemix) else 0
  
  params_init_nodematch <- list(
    mu = inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1]),
    beta_overall = params_fitted$beta_overall,
    K = params_init_inhom$K,
    beta_edges = params_fitted$beta_edges,
    node_lambda = params_fitted$node_lambda,
    CS_params = c(cs_structural, cs_nodematch)[1:n_cs_nodematch],
    vertex_categorical = params_fitted$vertex_categorical,
    vertex_categorical_levels = params_init_inhom$vertex_categorical_levels
  )
  
  # Ensure CS_params has correct length
  if (length(params_init_nodematch$CS_params) < n_cs_nodematch) {
    params_init_nodematch$CS_params <- c(params_init_nodematch$CS_params, 
                                          rep(0, n_cs_nodematch - length(params_init_nodematch$CS_params)))
  }
  
  params_init_nodematch <- hawkesGrowthNet:::repair_vertex_categorical_params(params_init_nodematch, eps = 1e-6)
  
  # Keep vertex_categorical - nodeMatch still needs it to identify node attributes
  # (nodeMatch uses it differently than nodeMix, but it's still required)
  
  # Create parscale for nodeMatch
  p_scale_nodematch <- c(
    beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
    setNames(rep(0.1, n_cs_nodematch), paste0("CS_params", seq_len(n_cs_nodematch))),
    vertex_categorical.gender.female = 0.1, vertex_categorical.gender.male = 0.1
  )
  
  cat("  Method: Nelder-Mead (max 500 iterations)\n")
  t_fit_nodematch <- proc.time()
  
  fit_inhom_nodematch <- tryCatch(
    fit_hawkesGrowthNet_inhom(
      params_init = params_init_nodematch,
      time_window = time_window_01,
      mark_filtration = net_raw,
      PMF_mark = PMF_mark_CS,
      mu_vec = inhom_bg$mu_vec,
      integral_bg = inhom_bg$integral_bg,
      formula_RHS = FORMULA_RHS_NODEMATCH,
      truncation = TRUNCATION,
      mark_decay = "activity",
      max_node_time = 1,
      method = "Nelder-Mead",
      maxit = MAX_ITER,
      trace = 1,
      reltol = 1e-8,
      verbose = FALSE,
      fixed_params = c("K", "mu"),
      parscale = p_scale_nodematch,
      cache_intensity = TRUE,
      cores = N_CORES
    ),
    error = function(e) { cat("  ERROR: Fit failed:", e$message, "\n"); NULL }
  )
  
  elapsed_fit_nodematch <- (proc.time() - t_fit_nodematch)[3]
  if (!is.null(fit_inhom_nodematch)) {
    cat("  Fit completed:", round(elapsed_fit_nodematch, 1), "s (", round(elapsed_fit_nodematch / 60, 1), "min)\n")
    cat("  Convergence:", fit_inhom_nodematch$fit$convergence, "\n")
    cat("  Iterations:", fit_inhom_nodematch$fit$counts[1], "\n")
  } else {
    cat("  Fit FAILED after", round(elapsed_fit_nodematch, 1), "s\n")
  }
  cat("  Step 2b total:", round((proc.time() - t_step_nodematch)[3], 1), "s\n\n")
} else {
  if (is.null(inhom_bg)) cat("  No inhomogeneous background; skipping nodeMatch fit\n")
  if (is.null(fit_inhom)) cat("  No nodeMix fit available; skipping nodeMatch fit\n")
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
                    wait_obs = NULL, wait_sim = NULL,
                    nodemix_obs = NULL, nodemix_sim = NULL)
GOF_results_nodematch <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                               geodist_obs = NULL, geodist_sim = NULL,
                               wait_obs = NULL, wait_sim = NULL,
                               nodemix_obs = NULL, nodemix_sim = NULL)

# GOF for nodeMix model
if (RUN_GOF && !is.null(fit_inhom)) {
  cat("  GOF for nodeMix model...\n")
  # For GOF simulations, use cond_intensity (not cond_intensity_inhom)
  # The average mu from inhom_bg will be used (computed in gof() function)
  GOF_results <- gof(
    fit = fit_inhom,
    net_obs = net_raw,
    params_init = params_init_inhom,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,  # Use homogeneous version for simulations
    formula_RHS = FORMULA_RHS,
    time_window = c(0, 0.05),
    truncation = TRUNCATION,
    mark_decay = "activity",
    max_node_time = 1,
    inhom_bg = inhom_bg,
    n_sim = N_GOF,
    cores = N_CORES,
    max_deg = 15,
    k_esp = 15,
    mu_multiplier = 5,
    verbose = TRUE
  )
} else {
  if (!RUN_GOF) cat("  RUN_GOF = FALSE; skipping\n")
  if (is.null(fit_inhom)) cat("  No nodeMix fit available; skipping GOF\n")
}

# GOF for nodeMatch model
if (RUN_GOF && !is.null(fit_inhom_nodematch)) {
  cat("\n  GOF for nodeMatch model...\n")
  
  # Reconstruct params_init for nodeMatch
  skel_nodematch_gof <- params_init_inhom
  skel_nodematch_gof$vertex_categorical_levels <- NULL
  params_init_nodematch_gof <- relist(fit_inhom_nodematch$fit$par, skeleton = skel_nodematch_gof)
  params_init_nodematch_gof$vertex_categorical_levels <- params_init_inhom$vertex_categorical_levels
  params_init_nodematch_gof$K <- params_init_inhom$K
  params_init_nodematch_gof$mu <- inhom_bg$integral_bg / (time_window_01[2] - time_window_01[1])
  # Restore names and repair parameters before GOF
  params_init_nodematch_gof <- hawkesGrowthNet:::reconstruct_vertex_categorical_names(
    params_init_nodematch_gof, params_init_inhom$vertex_categorical_levels)
  params_init_nodematch_gof <- hawkesGrowthNet:::repair_vertex_categorical_params(params_init_nodematch_gof, eps = 1e-6)
  
  # For GOF simulations, use cond_intensity (not cond_intensity_inhom)
  # The average mu from inhom_bg will be used (computed in gof() function)
  GOF_results_nodematch <- gof(
    fit = fit_inhom_nodematch,
    net_obs = net_raw,
    params_init = params_init_nodematch_gof,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,  # Use homogeneous version for simulations
    formula_RHS = FORMULA_RHS_NODEMATCH,
    time_window = c(0, 0.05),
    truncation = TRUNCATION,
    mark_decay = "activity",
    max_node_time = 1,
    inhom_bg = inhom_bg,
    n_sim = N_GOF,
    cores = N_CORES,
    max_deg = 15,
    k_esp = 15,
    mu_multiplier = 5,
    verbose = TRUE
  )
} else {
  if (!RUN_GOF) cat("  RUN_GOF = FALSE; skipping nodeMatch GOF\n")
  if (is.null(fit_inhom_nodematch)) cat("  No nodeMatch fit available; skipping nodeMatch GOF\n")
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
  fit_inhom_nodematch = fit_inhom_nodematch,
  params_init_inhom = params_init_inhom,
  FORMULA_RHS = FORMULA_RHS,
  FORMULA_RHS_NODEMATCH = FORMULA_RHS_NODEMATCH,
  fit_temporal = fit_temporal,
  ks_temporal_pval = ks_temporal_pval,
  realiz = realiz,
  windowT = windowT,
  GOF_results = GOF_results,
  GOF_results_nodematch = GOF_results_nodematch,
  N_GOF = N_GOF,
  SEARCH_STRING = SEARCH_STRING,
  time_window_01 = time_window_01
)
rds_path <- file.path(PKG_ROOT, "cluster_output", "results_openalex_full.RDS")
dir.create(file.path(PKG_ROOT, "cluster_output"), showWarnings = FALSE, recursive = TRUE)
saveRDS(save_list, rds_path)
cat("  Saved to", rds_path, "\n\n")

# =============================================================================
# 6. PAPER_OUTPUT: rehydrate and produce figures/tables
# =============================================================================
if (PAPER_OUTPUT) {
  cat("--- Step 6: Paper output (figures & tables) ---\n")
  t_step <- proc.time()
  dat <- readRDS(file.path(PKG_ROOT, "cluster_output", "results_openalex_full.RDS"))
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
      # Restore names and repair parameters before expanding
      pfit2 <- hawkesGrowthNet:::reconstruct_vertex_categorical_names(
        pfit2, params_init_inhom$vertex_categorical_levels)
      pfit2 <- hawkesGrowthNet:::repair_vertex_categorical_params(pfit2, eps = 1e-6)
      if (!is.null(pfit2$vertex_categorical$gender)) {
        levs <- params_init_inhom$vertex_categorical_levels$gender
        pgender <- expand_vertex_categorical_probs(pfit2$vertex_categorical$gender, levs)
        if (!is.null(pgender)) {
          cat("  Fitted gender proportions (n-1 expanded):\n"); print(pgender)
        } else {
          cat("  Fitted gender proportions: could not expand (invalid parameters)\n")
        }
      }
    }
  }
  if (!is.null(dat$fit_temporal)) {
    cat("  Temporal Hawkes par:\n"); print(fit_temporal$par)
    cat("  Temporal KS p-value:", ks_temporal_pval, "\n")
  }

  # GOF plots: use plots from gof() function if available, otherwise generate here
  if (!is.null(dat$GOF_results)) {
    # First, try to use plots from gof() function (if available)
    if (!is.null(dat$GOF_results$plots) && length(dat$GOF_results$plots) > 0) {
      cat("  Displaying GOF plots from gof() function...\n")
      if (!is.null(dat$GOF_results$plots$degree_plot)) print(dat$GOF_results$plots$degree_plot)
      if (!is.null(dat$GOF_results$plots$esp_plot)) print(dat$GOF_results$plots$esp_plot)
      if (!is.null(dat$GOF_results$plots$geodist_plot)) print(dat$GOF_results$plots$geodist_plot)
      if (!is.null(dat$GOF_results$plots$nodemix_plot)) print(dat$GOF_results$plots$nodemix_plot)
      if (!is.null(dat$GOF_results$plots$waiting_times_plot)) print(dat$GOF_results$plots$waiting_times_plot)
    } else if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(dat$GOF_results$degree_obs)) {
      # Fallback: generate plots here (legacy code)
      cat("  Generating GOF plots (legacy method)...\n")
      gof <- dat$GOF_results
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
    } else {
      cat("  ggplot2 not available; skipping GOF plots\n")
    }
  } else {
    cat("  No GOF results available for plotting\n")
  }
  cat("  Step 6 total:", round((proc.time() - t_step)[3], 1), "s\n\n")
}

# =============================================================================
# Total elapsed time
# =============================================================================
cat("=== OpenAlex study complete ===\n")
cat("  Total wall time:", round((proc.time() - t_total)[3], 1), "s (",
    round((proc.time() - t_total)[3] / 60, 1), "min)\n")
