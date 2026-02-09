# This script fits a LOLOG style model and a BA model to messaging data:
#
# RUNTIME ESTIMATE (16 cores, cluster):
#   One sim at TIME=50 is ~6 min (gives ~600 events); fit typically similar or longer.
#   Main study (SIMULATE):  N_SIMS=7  -> ~10-20 min total.
#   Consistency:            4 windows x 100 reps = 400 (sim+fit) pairs;
#                            ~8-15 min per run -> 400*10/16 ~ 4-4.5 h.
#   Explosive:               2 sims only, T=5 -> ~2-5 min.
#   Total (all blocks):      ~4.5-5 h. Set N_CORES via SLURM_CPUS_PER_TASK (e.g. 16).

library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(pbapply)
library(parallel)
library(doParallel)
library(R.utils)
library(network)
library(sna)
library(hash)
library(ernm)
library(hawkesGrowthNet)

# Paths: run from package root (directory containing inst/)
PKG_ROOT <- getwd()
CLUSTER_OUTPUT_DIR <- file.path(PKG_ROOT, "cluster_output")
dir.create(CLUSTER_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===================================================
# Change Statistic Mark Generation
# ===================================================
TIME <- 50
# takes ~ 6 minutes to simulate under this setting, gives ~600 events
params <- list(mu = 10,
               beta_overall = 2,
               K = 0.5,
               beta_edges = 1,
               node_lambda = 1,
               CS_params = c(-6.7,2,0.1,-0.1)
               )
TRUNCATION  = 100
SIMULATE = TRUE
PAPER_OUTPUT = TRUE
RUN_EXPLOSIVE <- TRUE
RUN_CONSISTENCY <- FALSE
MAX_ITER = 2000

N_SIMS <- 100 #should take ~ 30 minuts with 2 inner cores per fit
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 50))
# Core allocation: split roughly evenly between inner and outer
# N_CORES_INNER = cores per fit (for intensity cache parallelism via mclapply/fork).
# N_CORES_OUTER = parallel sim+fit workers (PSOCK cluster).
# Total cores used = N_CORES_OUTER * N_CORES_INNER.
# CS model benefits from inner parallelism (ERNM change stats are expensive).
# For 50 cores: split as 7 inner x 7 outer = 49 cores (efficient use)
# Or 5 inner x 10 outer = 50 cores (uses all cores)
N_CORES_INNER <- as.numeric(Sys.getenv("CORES_INNER", 2))
N_CORES_OUTER <- max(1L, floor(N_CORES / N_CORES_INNER))
# Ensure we don't exceed available cores
if (N_CORES_OUTER * N_CORES_INNER > N_CORES) {
  # Rebalance: prefer more outer workers for better parallelization
  N_CORES_OUTER <- max(1L, floor(sqrt(N_CORES)))
  N_CORES_INNER <- max(1L, floor(N_CORES / N_CORES_OUTER))
}
cat("Core allocation:", N_CORES_OUTER, "outer x", N_CORES_INNER, "inner =",
    N_CORES_OUTER * N_CORES_INNER, "total (of", N_CORES, "available)\n")

SEED <- 01267

make_cluster <- function(N_CORES){
  # setup the cluster:
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  # export libraries to cluster:
  clusterEvalQ(cl, {
    library(spatstat)
    library(ggplot2)
    library(dplyr)
    library(data.table)
    library(pbapply)
    library(parallel)
    library(doParallel)
    library(R.utils)
    library(ernm)
    library(network)
    library(sna)
    library(hash)
    library(hawkesGrowthNet)
  })
  clusterExport(cl, c("params",
                      "TIME",
                      "TRUNCATION",
                      "SEED",
                      "MAX_ITER"
                      ))
  return(cl)
}

if(SIMULATE){
  # make the cluster: use N_CORES_OUTER for outer parallelization
  t <- proc.time()
  cl <- make_cluster(N_CORES_OUTER)
  clusterExport(cl, c("N_CORES_INNER"))
  set.seed(SEED)
  
  # Ensure the cluster will be stopped no matter what.
  on.exit({
    if (!is.null(cl)) {
      stopCluster(cl)
    }
  }, add = TRUE)

  print("Commencing simulation:")
  print(paste("Using", N_CORES_OUTER, "outer workers,", N_CORES_INNER, "inner cores each"))
  sims <- parLapply(cl=cl,1:N_SIMS,function(x){
    results <- tryCatch({
      sim_hawkesGrowthNet(params =  params,
                          time_window = c(0,TIME),
                          PMF_mark = PMF_mark_CS,
                          cond_intensity = cond_intensity,
                          hashed_edges = T,
                          verbose = F,
                          mu_multiplier = 3,
                          joint_accept = F,
                          truncation = TRUNCATION,
                          formula_RHS = "edges  + triangles() + star(c(2,3))")
    }, error = function(e) {
      # Already inside parallel worker; just return NULL or partial data
      message("Error in sim_hawkesGrowthNet: ", e$message)
      return(NULL)
    })
    return(results)
  })
  print("Simulation took:")
  print(proc.time()-t)
  
  # only keep non null sims:
  sims <- sims[sapply(sims,length)!=0]
  
  fits <- NULL
  params_init <- list(mu = 10,
                      beta_overall = 1,
                      K = 0.5,
                      beta_edges = 1,
                      node_lambda = 1,
                      CS_params = c(-10,0,0,0)
  )
  clusterExport(cl, c("params_init", "N_CORES_INNER"))
  print("Commencing Fitting")
  print(paste("Using", N_CORES_OUTER, "outer workers,", N_CORES_INNER, "inner cores each"))
  t1 <- proc.time()
  fits <- parLapply(cl=cl,sims,function(x){
    fit <- tryCatch({
      fit_hawkesGrowthNet(
        params_init = params_init,
        time_window = c(0, TIME),
        mark_filtration = x$net,
        PMF_mark = PMF_mark_CS,
        formula_RHS = "edges + triangles + star(c(2,3))",
        trace = 0,
        maxit = MAX_ITER,
        truncation = TRUNCATION,
        fixed_params = c("K"),
        method = "Nelder-Mead",
        cores = N_CORES_INNER,
        cache_intensity = TRUE
      )
    }, error = function(e) {
      # Already inside parallel worker; just return NULL or partial data
      message("Error in fit_hawkesGrowthNet: ", e$message)
      return(e$message)
    })
    return(fit)
  })
  print("Fitting took:")
  print(proc.time()-t1)
  # Save the fits and final network for analysis:
  temp_hawkes_fits <- lapply(sims,function(x){
    fit <- fit_temporal_hawkes(params_init = list(mu = 0.1,
                                                  beta = 1,
                                                  K = 0.1),
                               realiz = data.frame(t = x$events$t,
                                                   n = rep(x$events$n,length(x$events$t))),
                               windowT = c(0,TIME),
                               trace = 0,
                               maxit = 1000
    )
    return(fit)
  })
  stopCluster(cl)
  cl <- NULL
  saveRDS(list(sims=sims,
               fits = fits,
               temp_hawkes_fits = temp_hawkes_fits,
               params = params,
               params_init = params_init),
          file = file.path(CLUSTER_OUTPUT_DIR, "results_CS.RDS"))
  print("Simulating and fitting took:")
  print((proc.time()-t)[3])
}

# ==============================================================================
# STUDY 1: Consistency Analysis (Sliding Window / Increasing T) - CS model
# ==============================================================================
# Goal: Show that as TIME increases, the variance of estimates decreases and
# means converge to truth.
# ==============================================================================

if(RUN_CONSISTENCY){

  # 1. Define Time Windows to test
  time_windows <- c(5, 10, 20, 30, 40, 50)
  N_SIMS_CONSISTENCY <- 50

  # Parameters (Standard/Stable regime) - CS model
  params_true <- list(mu = 10,
                      beta_overall = 2,
                      K = 0.5,
                      beta_edges = 1,
                      node_lambda = 1,
                      CS_params = c(-6.7, 2, 0.1, -0.1))

  # Setup Cluster — use N_CORES_OUTER workers, each fit gets N_CORES_INNER for intensity cache
  t_consistency_total <- proc.time()
  cat("=== Consistency Study (CS) ===\n")
  cat("  Time windows:", paste(time_windows, collapse = ", "), "\n")
  cat("  N_SIMS per window:", N_SIMS_CONSISTENCY, "\n")
  cat("  Core allocation:", N_CORES_OUTER, "outer x", N_CORES_INNER, "inner =",
      N_CORES_OUTER * N_CORES_INNER, "total (of", N_CORES, "available)\n")
  cat("Setting up cluster...\n")
  t_cluster <- proc.time()
  cl <- make_cluster(N_CORES_OUTER)
  clusterExport(cl, c("params_true", "TRUNCATION", "N_CORES_INNER", "MAX_ITER"))
  cat("  Cluster setup:", round((proc.time() - t_cluster)[3], 1), "s\n")

  # Storage for results
  consistency_results <- data.frame()

  for(curr_time in time_windows){
    t_window <- proc.time()
    cat("\n--- T =", curr_time, "(", which(time_windows == curr_time), "/",
        length(time_windows), ") ---\n")

    # Export current time to cluster
    clusterExport(cl, "curr_time", envir = environment())

    # Parallel Simulation & Fitting Loop
    cat("  Running", N_SIMS_CONSISTENCY, "sim+fit pairs:", N_CORES_OUTER, "parallel x",
        N_CORES_INNER, "inner cores...\n")
    t_simfit <- proc.time()
    res_list <- parLapply(cl = cl, X = 1:N_SIMS_CONSISTENCY, fun = function(i){

      # A. Simulate
      sim_res <- tryCatch({
        sim_hawkesGrowthNet(params = params_true,
                            time_window = c(0, curr_time),
                            PMF_mark = PMF_mark_CS,
                            cond_intensity = cond_intensity,
                            hashed_edges = TRUE,
                            mu_multiplier = 3,
                            verbose = FALSE,
                            truncation = TRUNCATION,
                            formula_RHS = "edges + triangles() + star(c(2,3))")
      }, error = function(e) return(NULL))

      if(is.null(sim_res)) return(NULL)

      # B. Fit - CS options: formula_RHS, fixed_params = c("K")
      # Initialize near true params + small noise for better convergence
      # Ensure all parameters are positive and finite
      params_init <- list(
        mu = max(0.1, params_true$mu * exp(rnorm(1, 0, 0.2))),
        beta_overall = max(0.1, params_true$beta_overall * exp(rnorm(1, 0, 0.2))),
        K = params_true$K,
        beta_edges = max(0.1, params_true$beta_edges * exp(rnorm(1, 0, 0.2))),
        node_lambda = max(0.1, params_true$node_lambda * exp(rnorm(1, 0, 0.2))),
        CS_params = params_true$CS_params + rnorm(length(params_true$CS_params), 0, 0.5)
      )
      # Ensure CS_params are finite
      params_init$CS_params[!is.finite(params_init$CS_params)] <- params_true$CS_params[!is.finite(params_init$CS_params)]

      fit_res <- tryCatch({
        fit_hawkesGrowthNet(params_init = params_init,
                            time_window = c(0, curr_time),
                            mark_filtration = sim_res$net,
                            PMF_mark = PMF_mark_CS,
                            formula_RHS = "edges + triangles + star(c(2,3))",
                            maxit = MAX_ITER,
                            truncation = TRUNCATION,
                            cache_intensity = TRUE,
                            verbose = FALSE,
                            fixed_params = c("K"),
                            cores = N_CORES_INNER,
                            method = "Nelder-Mead")
      }, error = function(e) return(NULL))

      if(is.null(fit_res) || is.null(fit_res$fit)) return(NULL)
      
      # Check if fit succeeded
      keep <- (!is.null(fit_res$fit) && 
               length(fit_res$fit) > 0 && 
               fit_res$fit$convergence == 0 &&
               all(is.finite(fit_res$fit$par)) &&
               !any(fit_res$fit$par > 100) && 
               length(fit_res$fit$par) >= 2 &&
               fit_res$fit$par[2] <= 10)

      # Return row
      return(data.frame(
        keep = keep,
        sim_id = i,
        time_window = curr_time,
        param = names(fit_res$fit$par),
        estimate = as.numeric(fit_res$fit$par),
        true_value = as.numeric(unlist(params_true)[names(fit_res$fit$par)])
      ))
    })
    elapsed_simfit <- (proc.time() - t_simfit)[3]

    # Bind results
    res_df <- do.call(rbind, res_list)
    n_success <- length(unique(res_df$sim_id))
    n_fail <- N_SIMS_CONSISTENCY - n_success
    consistency_results <- rbind(consistency_results, res_df)

    elapsed_window <- (proc.time() - t_window)[3]
    cat("  Sim+fit:", round(elapsed_simfit, 1), "s |",
        "Success:", n_success, "/", N_SIMS_CONSISTENCY,
        "(", n_fail, "failed)\n")
    cat("  Window total:", round(elapsed_window, 1), "s (",
        round(elapsed_window / 60, 1), "min)\n")
    elapsed_so_far <- (proc.time() - t_consistency_total)[3]
    remaining_windows <- length(time_windows) - which(time_windows == curr_time)
    cat("  Elapsed so far:", round(elapsed_so_far / 60, 1), "min |",
        "Windows remaining:", remaining_windows, "\n")
  }

  stopCluster(cl)
  elapsed_consistency <- (proc.time() - t_consistency_total)[3]
  cat("\n=== Consistency Study (CS) complete ===\n")
  cat("  Total time:", round(elapsed_consistency / 60, 1), "min (",
      round(elapsed_consistency / 3600, 2), "h)\n")

  # ==========================
  # Visualization
  # ==========================
  prop_keep <- consistency_results %>%
    group_by(time_window, param) %>%
    summarise(prop_keep = sum(keep) / N_SIMS_CONSISTENCY)
  print(prop_keep)

  # Calculate Bias and RMSE (use all runs; keep filter commented)
  summary_stats <- consistency_results %>%
    # filter(keep == TRUE) %>%
    group_by(time_window, param) %>%
    summarise(
      mean_est = mean(estimate),
      sd_est = sd(estimate),
      rmse = sqrt(mean((estimate - true_value)^2)),
      true_val = mean(true_value),
    )

  print(summary_stats,n=100)

  # Plot 1: Boxplots of convergence (use all runs; keep filter commented)
  p_cons <- ggplot(consistency_results, aes(x = factor(time_window), y = estimate)) +  # %>% filter(keep == TRUE) 
    geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightblue") +
    geom_jitter(width = 0.2, alpha = 0.3) +
    geom_hline(aes(yintercept = true_value), color = "red", linetype = "dashed", size = 1) +
    facet_wrap(~param, scales = "free_y") +
    labs(title = "Parameter Consistency vs Time Window (T) - CS model",
         subtitle = "Red dashed line indicates true parameter value",
         x = "Time Window Length (T)",
         y = "Parameter Estimate") +
    theme_minimal()

  print(p_cons)

  # Plot 2: RMSE decay (The "Getting Better" plot)
  p_rmse <- ggplot(summary_stats, aes(x = time_window, y = rmse)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    facet_wrap(~param, scales = "free_y") +
    labs(title = "RMSE Decay as Data Increases (CS model)",
         x = "Time Window Length (T)",
         y = "Root Mean Squared Error") +
    theme_bw()

  print(p_rmse)
}

# ==============================================================================
# STUDY 2: Explosive Regime Analysis - CS model
# ==============================================================================
# Goal: Set beta_overall and beta_edges -> 0.
# 1. beta_overall -> 0 means the memory of the process never decays.
#    If K > 0, the integral of intensity diverges (Explosive / Super-critical).
# 2. beta_edges -> 0 means the preferential attachment logic considers ALL past
#    nodes equally (no time decay on degree relevance).
# ==============================================================================

if(RUN_EXPLOSIVE){

  # Define Explosive Parameters - CS model (include node_lambda, CS_params)
  params_explosive <- list(
    mu = 10,
    beta_overall = 0.1,
    K = 0.99,
    beta_edges = 0.1,
    node_lambda = 1,
    CS_params = c(-6.7, 2, 0.1, -0.1)
  )

  # Compare with Stable Parameters - CS model
  params_stable <- list(
    mu = 10,
    beta_overall = 2.0,
    K = 0.5,
    beta_edges = 1.0,
    node_lambda = 1,
    CS_params = c(-6.7, 2, 0.1, -0.1)
  )

  print("Simulating Explosive Regime (CS model)...")

  T_explode <- 5

  # Simulate Explosive
  sim_exp <- sim_hawkesGrowthNet(params = params_explosive,
                                 time_window = c(0, T_explode),
                                 PMF_mark = PMF_mark_CS,
                                 cond_intensity = cond_intensity,
                                 hashed_edges = TRUE,
                                 verbose = FALSE,
                                 truncation = TRUNCATION,
                                 formula_RHS = "edges + triangles() + star(c(2,3))",
                                 mu_multiplier = 50)

  print("Simulating Stable Regime (CS model)...")
  sim_stable <- sim_hawkesGrowthNet(params = params_stable,
                                    time_window = c(0, T_explode),
                                    PMF_mark = PMF_mark_CS,
                                    cond_intensity = cond_intensity,
                                    hashed_edges = TRUE,
                                    verbose = FALSE,
                                    truncation = TRUNCATION,
                                    formula_RHS = "edges + triangles() + star(c(2,3))",
                                    mu_multiplier = 5)

  # ==========================
  # Visualization: Cumulative Events
  # ==========================
  df_exp <- data.frame(t = sim_exp$events$t,
                       N = seq_len(length(sim_exp$events$t)),
                       Type = "Explosive (K ~ 1)")

  df_stable <- data.frame(t = sim_stable$events$t,
                          N = seq_len(length(sim_stable$events$t)),
                          Type = "Stable (K ~ 0.25)")

  df_compare <- rbind(df_exp, df_stable)

  p_expl <- ggplot(df_compare, aes(x = t, y = N, color = Type)) +
    geom_line(size = 1.2) +
    labs(title = "Explosive vs Stable Process Dynamics (CS model)",
         x = "Time",
         y = "Cumulative Number of Events (N)") +
    theme_minimal() +
    theme(legend.position = "bottom")

  print(p_expl)

  # ==========================
  # Node growth over time (number of nodes at each event time)
  # ==========================
  node_times_exp <- sim_exp$net %v% "time"
  node_times_stable <- sim_stable$net %v% "time"
  n_nodes_exp <- sapply(sim_exp$events$t, function(t) sum(node_times_exp <= t))
  n_nodes_stable <- sapply(sim_stable$events$t, function(t) sum(node_times_stable <= t))
  df_nodes_exp <- data.frame(t = sim_exp$events$t, N_nodes = n_nodes_exp, Type = "Explosive (K ~ 1)")
  df_nodes_stable <- data.frame(t = sim_stable$events$t, N_nodes = n_nodes_stable, Type = "Stable (K ~ 0.25)")
  df_nodes_compare <- rbind(df_nodes_exp, df_nodes_stable)
  p_node_growth <- ggplot(df_nodes_compare, aes(x = t, y = N_nodes, color = Type)) +
    geom_line(size = 1.2) +
    labs(title = "Node Growth Over Time: Explosive vs Stable (CS)",
         x = "Time",
         y = "Number of Nodes") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p_node_growth)

  # ==========================
  # Visualization: Network Structure Impact
  # ==========================
  op <- par(mfrow = c(1, 2))

  plot(sim_stable$net, main = "Stable Network (CS)\n(Recent Activity Matters)",
       vertex.cex = 0.5, edge.col = "gray")

  plot(sim_exp$net, main = "Explosive/Memory Network (CS)\n(History Never Dies)",
       vertex.cex = 0.5, edge.col = "gray")

  par(op)

  max_deg_stable <- max(degree(sim_stable$net))
  max_deg_exp <- max(degree(sim_exp$net))

  print(paste("Max Degree Stable:", max_deg_stable))
  print(paste("Max Degree Explosive:", max_deg_exp))
}

# ==============================================================================
# Save full state at end for re-hydration (main + consistency + explosive if run)
# ==============================================================================
save_list <- list()
if(exists("sims") && !is.null(sims)){
  save_list$sims <- sims
  save_list$fits <- fits
  save_list$temp_hawkes_fits <- temp_hawkes_fits
  save_list$params <- params
  save_list$params_init <- params_init
  save_list$TIME <- TIME
  save_list$N_SIMS <- N_SIMS
}
if(exists("consistency_results")){
  save_list$consistency_results <- consistency_results
  save_list$summary_stats <- summary_stats
  save_list$p_cons <- p_cons
  save_list$p_rmse <- p_rmse
  save_list$N_SIMS_CONSISTENCY <- N_SIMS_CONSISTENCY
  save_list$time_windows <- time_windows
}
if(exists("sim_exp")){
  save_list$sim_exp <- sim_exp
  save_list$sim_stable <- sim_stable
  save_list$df_compare <- df_compare
  save_list$p_expl <- p_expl
  save_list$df_nodes_compare <- df_nodes_compare
  save_list$p_node_growth <- p_node_growth
  save_list$max_deg_stable <- max_deg_stable
  save_list$max_deg_exp <- max_deg_exp
}
saveRDS(save_list, file.path(CLUSTER_OUTPUT_DIR, "results_CS_full.RDS"))
print(paste("Saved full state to", file.path(CLUSTER_OUTPUT_DIR, "results_CS_full.RDS")))

# ==============================================================================
# PAPER OUTPUT (at end: main study + consistency + explosive)
# Re-hydrate from results_CS_full.RDS and produce all figures/tables.
# ==============================================================================
if(PAPER_OUTPUT){
  dat <- readRDS(file.path(CLUSTER_OUTPUT_DIR, "results_CS_full.RDS"))
  list2env(dat, envir = .GlobalEnv)

  # ---------- Main study output ----------
  if(!is.null(dat$sims)){
    net_stats <- do.call(rbind, lapply(sims, function(s){
      net <- s$net
      degs <- ernm::calculateStatistics(net ~ degree(0:20,"in"))
      esps <- ernm::calculateStatistics(net ~ esp(0:20))
      tmp <- as.data.frame(cbind(c(degs,esps), rep(0:20, times = 2), c(rep("degree",21), rep("esp",21))))
      rownames(tmp) <- NULL
      return(tmp)
    }))
    names(net_stats) <- c("value","var","type")
    net_stats$value <- as.numeric(net_stats$value)
    net_stats$var <- as.numeric(net_stats$var)
    net_stats <- net_stats[net_stats$var <= 15,]
    net_stats <- net_stats %>% mutate(var = factor(var, levels = sort(unique(var))))

    deg_plot <- ggplot(net_stats[net_stats$type == "degree",], aes(x = var, y = value)) +
      geom_boxplot() + labs(title = "Degree Distribution", x = "Degree", y = "Value") + theme_minimal()
    print(deg_plot)
    esp_plot <- ggplot(net_stats[net_stats$type == "esp",], aes(x = var, y = value)) +
      geom_boxplot() + labs(title = "ESP Distribution", x = "ESP", y = "Value") + theme_minimal()
    print(esp_plot)

    mean_degs <- sapply(sims, function(s) mean(degree(s$net, gmode = "graph")))
    mean_deg_df <- data.frame(mean_deg = mean_degs)
    hist(mean_deg_df$mean_deg, main = "Histogram of Mean Degrees", xlab = "Mean Degree", breaks = 10)
    abline(v = mean(mean_deg_df$mean_deg), col = "red", lwd = 2)

    keep <- which(sapply(fits, function(x){ length(x) != 0 & x$fit$convergence == 0 & !any(x$fit$par > 100) & !any(x$fit$par[2] > 10) }))
    print(paste0("keeping ", length(keep), " of ", length(fits), " fits"))
    estims <- do.call(rbind, lapply(seq_along(keep), function(i){
      sim_idx <- keep[i]
      est_df <- as.data.frame(t(fits[[sim_idx]]$fit$par), names = names(fits[[sim_idx]]$fit$par))
      est_df$sim_id <- sim_idx
      return(est_df)
    }))
    params_vec <- unlist(params)[names(params) %in% colnames(estims)]
    params_init_vec <- unlist(params_init)[names(params_init) %in% colnames(estims)]
    par_estim <- estims[, -dim(estims)[2]]
    results <- data.frame(mean = colMeans(par_estim), sd = apply(par_estim, 2, sd), true = params_vec[colnames(par_estim)], init = params_init_vec[colnames(par_estim)])
    print(results)

    estim_long <- data.frame(param = rep(colnames(par_estim), each = nrow(par_estim)), estimate = c(as.matrix(par_estim)))
    ref_lines <- data.frame(param = colnames(par_estim), true = as.numeric(results["true",]), init = as.numeric(results["init",]))
    p_est_dist <- ggplot(estim_long, aes(x = estimate)) +
      geom_histogram(bins = 20, fill = "lightblue", alpha = 0.7) +
      geom_vline(data = ref_lines, aes(xintercept = true), color = "red", linetype = "dashed", linewidth = 1) +
      geom_vline(data = ref_lines, aes(xintercept = init), color = "darkgreen", linetype = "dotted", linewidth = 1) +
      facet_wrap(~param, scales = "free") +
      labs(title = "Distribution of parameter estimates (CS)", subtitle = "Red dashed = true; green dotted = init", x = "Estimate") + theme_minimal()
    print(p_est_dist)

    comps <- lapply(sims, function(x) compensators_hawkesGrowthNet(params = params, mark_filtration = x$net, time_window = c(0, TIME)))
    marked_p_vals <- mapply(seq_along(keep), FUN = function(i){
      sim_idx <- keep[i]
      sim <- sims[[sim_idx]]
      fit_par <- fits[[sim_idx]]$fit$par
      times <- get_times(sim$net)$times
      ks_test_pval_temporal(realiz = data.frame(t = times, n = rep(length(times), length(times))), windowT = c(0, TIME), hawkes_par = fit_par)
    })
    temp_p_vals <- mapply(sims, temp_hawkes_fits, FUN = function(x,y){
      ks_test_pval_temporal(realiz = data.frame(t = x$events$t, n = rep(x$events$n, length(x$events$t))), windowT = c(0, TIME), hawkes_par = y$par)
    })
    print(mean(marked_p_vals))
    print(mean(temp_p_vals))
  }

  # ---------- Consistency study output ----------
  if(!is.null(dat$consistency_results)){
    prop_keep <- consistency_results %>% group_by(time_window, param) %>% summarise(prop_keep = sum(keep) / N_SIMS_CONSISTENCY)
    print(prop_keep)
    print(summary_stats)
    print(p_cons)
    print(p_rmse)
  }

  # ---------- Explosive study output ----------
  if(!is.null(dat$sim_exp)){
    print(p_expl)
    if(!is.null(dat$p_node_growth)) print(p_node_growth)
    op <- par(mfrow = c(1, 2))
    plot(sim_stable$net, main = "Stable Network (CS)\n(Recent Activity Matters)", vertex.cex = 0.5, edge.col = "gray")
    plot(sim_exp$net, main = "Explosive/Memory Network (CS)\n(History Never Dies)", vertex.cex = 0.5, edge.col = "gray")
    par(op)
    print(paste("Max Degree Stable:", max_deg_stable))
    print(paste("Max Degree Explosive:", max_deg_exp))
  }
}
