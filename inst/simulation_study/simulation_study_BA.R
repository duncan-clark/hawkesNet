# RUNTIME ESTIMATE (16 cores, cluster):
#   One sim at TIME=50 is ~6 min (gives ~600 events); fit typically similar or longer.
#   Main study (SIMULATE):  N_SIMS=100  -> ~1-2 h total.
#   Consistency:            10 windows x 50 reps = 500 (sim+fit) pairs;
#                            T=1000 window is expensive (~hours per rep).
#                            Expect 12-24+ h total depending on cores.
#   Explosive:               2 sims only, T=5 -> ~2-5 min.
#   Total (all blocks):      ~15-30 h. Set N_CORES via SLURM_CPUS_PER_TASK (e.g. 64+).

library(spatstat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(pbapply)
library(parallel)
library(doParallel)
library(R.utils)
library(network)
library(sna)
library(hash)
library(ernm)
library(hawkesNet)

# Paths: run from package root (directory containing inst/).
# Under SLURM, use submit dir so path stays valid if getwd() breaks later.
PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) Sys.getenv("SLURM_SUBMIT_DIR") else getwd()
source(file.path(PKG_ROOT, "inst", "resolve_output_dir.R"))
CLUSTER_OUTPUT_DIR <- hawkesnet_resolve_output_dir(PKG_ROOT)
cat("CLUSTER_OUTPUT_DIR:", CLUSTER_OUTPUT_DIR, "\n")

# ===================================================
# BA simulation and fitting
# ===================================================
ON_SLURM <- nzchar(Sys.getenv("SLURM_JOB_ID"))

# Defaults differ for interactive (RStudio on cluster) vs SLURM batch jobs
if (ON_SLURM) {
  TIME     <- 50
  N_SIMS   <- 100
  N_CORES  <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 50))
} else {
  # Interactive / RStudio defaults — fast iteration for debugging
  TIME     <- 5
  N_SIMS   <- 8
  N_CORES  <- 64L
}

params <- list(mu = 10,
               beta_overall = 1,
               K = 0.5,
               beta_edges = 1,
               m = 1  # expected edges per event (Poisson); m=1 gives mean degree ~2
)
TRUNCATION  <- 100
SIMULATE <- TRUE
PAPER_OUTPUT <- TRUE
RUN_EXPLOSIVE <- ON_SLURM          # skip explosive in interactive mode
RUN_CONSISTENCY <- ON_SLURM  # Run on cluster batch jobs only (slow)
MAX_ITER <- if (ON_SLURM) 5000 else 2000
cat(sprintf("Mode: %s | TIME=%d | N_SIMS=%d | N_CORES=%d | MAX_ITER=%d\n",
    if (ON_SLURM) "SLURM batch" else "Interactive (RStudio)", TIME, N_SIMS, N_CORES, MAX_ITER))

if (nzchar(Sys.getenv("CORES_OVERRIDE"))) {
  N_CORES <- as.numeric(Sys.getenv("CORES_OVERRIDE"))
}

# --- Core allocation strategy ---
# Simulation is single-threaded: use all cores as outer workers (phase 1).
# Fitting has two phases:
#   (a) Intensity cache build: uses inner cores (mclapply fork). This is ~5-30s.
#   (b) Optimization (optim): single-threaded. This is the majority of runtime.
# Therefore, fitting benefits most from many outer workers with modest inner cores.
# Too many PSOCK workers create overhead and memory bloat. Sweet spot: ~25 outer.
CORES_OUTER_ENV <- Sys.getenv("CORES_OUTER", "")
if (nzchar(CORES_OUTER_ENV)) {
  N_CORES_OUTER <- as.numeric(CORES_OUTER_ENV)
  N_CORES_INNER <- max(1L, floor(N_CORES / N_CORES_OUTER))
} else {
  # Default: 25 outer workers, each with ~5 inner cores (for 128-core machine)
  N_CORES_OUTER <- min(25L, N_CORES)
  N_CORES_INNER <- max(1L, floor(N_CORES / N_CORES_OUTER))
}
# Cap outer workers at 60 for R socket limits (128 connections)
N_CORES_OUTER <- min(N_CORES_OUTER, 60L)
N_CORES_INNER <- max(1L, floor(N_CORES / N_CORES_OUTER))

cat("Core allocation:", N_CORES_OUTER, "outer x", N_CORES_INNER, "inner =",
    N_CORES_OUTER * N_CORES_INNER, "total (of", N_CORES, "available)\n")

SEED <- 1267

# parscale: match param magnitudes so Nelder-Mead simplex steps are proportionate
# Defined here (not inside SIMULATE block) so consistency study can also use it.
# IMPORTANT: Names must match flat_par names exactly.
p_scale <- c(mu = 1, beta_overall = 0.1, K = 0.1, beta_edges = 0.1, m = 0.1)

make_cluster <- function(n_workers) {
  cl <- makeCluster(n_workers, outfile = "")
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(hawkesNet)
      library(ernm)
      library(network)
      library(sna)
      library(data.table)
      library(hash)
    })
    plogis <- stats::plogis
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    }
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  })
  clusterExport(cl, c("params", "TIME", "TRUNCATION", "SEED", "MAX_ITER"))
  return(cl)
}

if(SIMULATE){
  # 1. Simulation Step: Each sim is single-threaded (~150s at T=50).
  # Use N_CORES workers (capped at 60 for R socket safety).
  t <- proc.time()
  N_SIM_WORKERS <- as.numeric(Sys.getenv("SIM_WORKERS", min(N_CORES, 60L)))
  N_SIM_WORKERS <- max(1L, min(N_SIM_WORKERS, 60L))
  cat("Commencing simulation using", N_SIM_WORKERS, "parallel workers...\n")
  cl_sim <- makeCluster(N_SIM_WORKERS, outfile = "")
  clusterEvalQ(cl_sim, {
    suppressPackageStartupMessages({
      library(hawkesNet)
      library(ernm)
      library(network)
      library(data.table)
    })
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    }
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  })
  clusterExport(cl_sim, c("params", "TIME", "TRUNCATION", "SEED"))
  
  set.seed(SEED)
  sims <- parLapply(cl=cl_sim, 1:N_SIMS, function(x){
    t_start <- proc.time()
    res <- tryCatch({
      sim_hawkesNet(params = params,
                    time_window = c(0, TIME),
                    PMF_mark = PMF_mark_BA,
                    cond_intensity = cond_intensity,
                    hashed_edges = TRUE,
                    verbose = FALSE,
                    mu_multiplier = 3,
                    joint_accept = FALSE,
                    truncation = TRUNCATION)
    }, error = function(e) {
      message("Error in sim_hawkesNet: ", e$message)
      return(NULL)
    })
    t_end <- proc.time()
    if (!is.null(res)) {
      message(sprintf("  [Sim Worker %d] Rep %d complete in %.1f s (events=%d)", 
                      Sys.getpid(), x, (t_end - t_start)[3], length(res$events$t)))
    }
    return(res)
  })
  stopCluster(cl_sim)
  cat("Simulation took:", round((proc.time() - t)[3], 1), "s\n")
  
  # only keep non null sims:
  sims <- sims[sapply(sims, length) != 0]
  gc()
  
  # 2. Fitting Step: Nested parallelism (PSOCK outer x fork inner).
  # Each outer worker runs one fit; inner cores used for intensity cache (mclapply fork).
  # BLAS threads are pre-set to 1 in make_cluster() so forked grandchildren are safe.
  cat("\nCommencing Fitting\n")
  cat("Using nested parallelism:", N_CORES_OUTER, "outer workers x", N_CORES_INNER, "inner cores each\n")
  
  cl_fit <- make_cluster(N_CORES_OUTER)
  clusterExport(cl_fit, c("N_CORES_INNER"))
  
  # Init values
  params_init <- list(mu = 10,
                      beta_overall = 1,
                      K = 0.5,
                      beta_edges = 1,
                      m = 1
  )
  clusterExport(cl_fit, c("params_init", "p_scale", "TIME", "MAX_ITER", "TRUNCATION"))
  
  t1 <- proc.time()
  fits <- parLapply(cl=cl_fit, sims, function(x){
    worker_id <- Sys.getpid()
    n_events <- length(x$events$t)
    message(sprintf("  [Outer Worker %d] Starting fit for sim with %d events...", worker_id, n_events))
    t_start <- proc.time()
    
    fit <- tryCatch({
      fit_hawkesNet(
        params_init = params_init,
        time_window = c(0, TIME),
        mark_filtration = x$net,
        PMF_mark = PMF_mark_BA,
        trace = 0,
        maxit = MAX_ITER,
        truncation = TRUNCATION,
        method = "Nelder-Mead",
        parscale = p_scale,
        cores = N_CORES_INNER,
        cache_intensity = TRUE,
        combine_intensity = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      message(sprintf("  [Outer Worker %d] ERROR: %s", worker_id, e$message))
      return(e$message)
    })
    
    t_end <- proc.time()
    elapsed <- (t_end - t_start)[3]
    
    if (is.list(fit)) {
      message(sprintf("  [Outer Worker %d] Fit complete in %.1f s (events=%d, conv=%d, iters=%d, ll=%.2f)", 
                      worker_id, elapsed, n_events, fit$fit$convergence, fit$fit$counts[1], -fit$fit$value))
      if (!is.null(fit$intens_funcs)) fit$intens_funcs <- NULL
    }
    return(fit)
  })
  
  cat("Fitting took:", round((proc.time() - t1)[3], 1), "s\n")
  
  # 3. Temporal Hawkes: Use all cores again
  cat("\nFitting Temporal Hawkes models...\n")
  temp_hawkes_fits <- parLapply(cl_fit, sims, function(x){
    fit_temporal_hawkes(params_init = list(mu = 0.1, beta = 1, K = 0.1),
                        realiz = data.frame(t = x$events$t, n = rep(x$events$n, length(x$events$t))),
                        windowT = c(0, TIME), trace = 0, maxit = 1000)
  })
  
  stopCluster(cl_fit)
  
  saveRDS(list(sims=sims,
               fits = fits,
               temp_hawkes_fits = temp_hawkes_fits,
               params = params,
               params_init = params_init),
          file = file.path(CLUSTER_OUTPUT_DIR, "results_BA.RDS"))
  cat("Total SIMULATE block took:", round((proc.time() - t)[3], 1), "s\n")
  gc()
}

# ==============================================================================
# STUDY 1: Consistency Analysis (Sliding Window / Increasing T) - BA model
# ==============================================================================
# Goal: Show that as TIME increases, the variance of estimates decreases and
# means converge to truth.
# ==============================================================================

if(RUN_CONSISTENCY){

  # 1. Define Time Windows to test (include T=20 for consistency with diagnostic, T=1000 for large-sample)
  time_windows <- c(5, 10, 20, 25, 50, 75, 100, 200, 500, 1000)
  N_SIMS_CONSISTENCY <- 50

  # Parameters (Standard/Stable regime) - BA model
  params_true <- list(mu = 10,
                      beta_overall = 1,
                      K = 0.5,
                      beta_edges = 1,
                      m = 1)

    # Setup Cluster — nested parallelism (PSOCK outer x fork inner).
    # Each worker does sim + fit sequentially. Use same allocation as main study.
    N_CONS_OUTER <- min(N_SIMS_CONSISTENCY, N_CORES_OUTER, 60L)
    N_CONS_INNER <- max(1L, floor(N_CORES / N_CONS_OUTER))
    
    t_consistency_total <- proc.time()
    cat("=== Consistency Study (BA) ===\n")
    cat("  Time windows:", paste(time_windows, collapse = ", "), "\n")
    cat("  N_SIMS per window:", N_SIMS_CONSISTENCY, "\n")
    cat("  Core allocation:", N_CONS_OUTER, "outer workers x", N_CONS_INNER, "inner =",
        N_CONS_OUTER * N_CONS_INNER, "total (of", N_CORES, "available)\n")
    cat("Setting up cluster...\n")
    t_cluster <- proc.time()
    cl <- make_cluster(N_CONS_OUTER)
    clusterExport(cl, c("params_true", "TRUNCATION", "N_CONS_INNER", "MAX_ITER", "p_scale"))
    # (Inside the loop, workers will use N_CONS_INNER for fitting)
    cat("  Cluster setup:", round((proc.time() - t_cluster)[3], 1), "s\n")
    
    # Storage for results
    consistency_results <- data.frame()
    consistency_example_sims <- list()
    
    for(curr_time in time_windows){
      t_window <- proc.time()
      cat("\n--- T =", curr_time, "(", which(time_windows == curr_time), "/",
          length(time_windows), ") ---\n")
      
      # Export current time to cluster
      clusterExport(cl, "curr_time", envir = environment())
      
      # Parallel Simulation & Fitting Loop
      cat("  Running", N_SIMS_CONSISTENCY, "sim+fit pairs:", N_CONS_OUTER, "parallel x",
          N_CONS_INNER, "inner cores...\n")
      t_simfit <- proc.time()
      res_list <- parLapply(cl = cl, X = 1:N_SIMS_CONSISTENCY, fun = function(i){
        worker_id <- Sys.getpid()
        t_start_pair <- proc.time()
        
        # A. Simulate
        sim_res <- tryCatch({
          sim_hawkesNet(params = params_true,
                              time_window = c(0, curr_time),
                              PMF_mark = PMF_mark_BA,
                              cond_intensity = cond_intensity,
                              hashed_edges = TRUE,
                              mu_multiplier = 3,
                              verbose = FALSE,
                              truncation = TRUNCATION)
        }, error = function(e) return(NULL))
        
        if(is.null(sim_res)) return(NULL)
        
        # B. Fit
        # Initialize near true params + small noise for better convergence.
        params_init <- list(
          mu = max(0.1, params_true$mu * exp(rnorm(1, 0, 0.2))),
          beta_overall = max(0.1, params_true$beta_overall * exp(rnorm(1, 0, 0.2))),
          K = min(0.99, max(0.01, params_true$K * exp(rnorm(1, 0, 0.2)))),
          beta_edges = max(0.1, params_true$beta_edges * exp(rnorm(1, 0, 0.2))),
          m = max(0.1, params_true$m * exp(rnorm(1, 0, 0.2)))
        )
        
        fit_res <- tryCatch({
          fit_hawkesNet(params_init = params_init,
                              time_window = c(0, curr_time),
                              mark_filtration = sim_res$net,
                              PMF_mark = PMF_mark_BA,
                              maxit = MAX_ITER,
                              truncation = TRUNCATION,
                              cache_intensity = TRUE,
                              combine_intensity = TRUE,
                              verbose = FALSE,
                              trace = 0,
                              parscale = p_scale,
                              cores = N_CONS_INNER,
                              method = "Nelder-Mead")
      }, error = function(e) return(NULL))
        
        if(is.null(fit_res) || is.null(fit_res$fit)) return(NULL)
        
        t_end_pair <- proc.time()
        elapsed_pair <- (t_end_pair - t_start_pair)[3]
        message(sprintf("  [Consistency Worker %d] T=%.1f Rep %d complete in %.1f s (events=%d, conv=%d)", 
                        worker_id, curr_time, i, elapsed_pair, length(sim_res$events$t), fit_res$fit$convergence))
        
        # Check if fit succeeded
        par_bo_idx <- which(names(fit_res$fit$par) == "beta_overall")
        keep <- (!is.null(fit_res$fit) && 
                 length(fit_res$fit) > 0 && 
                 fit_res$fit$convergence == 0 &&
                 all(is.finite(fit_res$fit$par)) &&
                 !any(fit_res$fit$par > 100) && 
                 length(fit_res$fit$par) >= 2 &&
                 (length(par_bo_idx) == 0 || fit_res$fit$par[par_bo_idx] <= 10))
        
        # Map true values safely
        par_names <- names(fit_res$fit$par)
        n_pars <- length(par_names)
        true_vals <- setNames(numeric(n_pars), par_names)
        
        # Helper to safely map
        safe_map <- function(target_nm, val) {
          if (target_nm %in% par_names) true_vals[target_nm] <<- val
        }
        
        safe_map("mu", params_true$mu)
        safe_map("beta_overall", params_true$beta_overall)
        safe_map("K", params_true$K)
        safe_map("beta_edges", params_true$beta_edges)
        safe_map("m", params_true$m)
        
        # Final safety check on dimensions before returning
        estimates <- as.numeric(fit_res$fit$par)
        t_vals <- as.numeric(true_vals)
        
        if (length(estimates) != n_pars || length(t_vals) != n_pars) {
           message(sprintf("  [Consistency Worker %d] DIMENSION MISMATCH: names=%d, est=%d, true=%d", 
                           worker_id, n_pars, length(estimates), length(t_vals)))
           return(NULL)
        }

        result <- list(
          df = data.frame(
            keep = keep,
            sim_id = i,
            time_window = curr_time,
            event_count = length(sim_res$events$t),
            param = par_names,
            estimate = estimates,
            true_value = t_vals
          ),
          example_sim = if (i == 1) sim_res else NULL
        )
        return(result)
      })
      elapsed_simfit <- (proc.time() - t_simfit)[3]
      
      # Bind results (each element is list(df=..., example_sim=...))
      non_null <- Filter(Negate(is.null), res_list)
      res_df <- do.call(rbind, lapply(non_null, `[[`, "df"))
      
      # Retain the first non-NULL example sim for this time window
      ex_sim <- NULL
      for (r in non_null) {
        if (!is.null(r$example_sim)) { ex_sim <- r$example_sim; break }
      }
      if (!is.null(ex_sim)) {
        consistency_example_sims[[as.character(curr_time)]] <- ex_sim
      }
      
      if (is.null(res_df) || nrow(res_df) == 0) {
        cat("  WARNING: No successful fits in this window!\n")
        next
      }
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
    cl <- NULL
    elapsed_consistency <- (proc.time() - t_consistency_total)[3]
    cat("\n=== Consistency Study (BA) complete ===\n")
    cat("  Total time:", round(elapsed_consistency / 60, 1), "min (",
        round(elapsed_consistency / 3600, 2), "h)\n")
    
    # --- Cleanup: remove consistency temporaries ---
    cleanup_vars <- c("t_consistency_total", "t_cluster", "t_window", "t_simfit",
                       "elapsed_simfit", "elapsed_window", "res_list", "res_df",
                       "n_success", "n_fail", "remaining_windows", "elapsed_so_far")
    rm(list = intersect(cleanup_vars, ls()), envir = environment())
    gc()
    
    # ==========================
    # Visualization
    # ==========================
    if (nrow(consistency_results) > 0) {
      prop_keep <- consistency_results %>%
        group_by(time_window, param) %>%
        summarise(prop_keep = sum(keep) / N_SIMS_CONSISTENCY, .groups = "drop")
      print(prop_keep)
      
      # Filter to converged fits only
      converged_results <- consistency_results %>% filter(keep == TRUE)
      
      if (nrow(converged_results) > 0) {
        # Calculate Bias and RMSE (converged fits only)
        summary_stats <- converged_results %>%
          group_by(time_window, param) %>%
          summarise(
            mean_est = mean(estimate),
            sd_est = sd(estimate),
            rmse = sqrt(mean((estimate - true_value)^2)),
            true_val = mean(true_value),
            mean_events = mean(event_count),
            n_conv = n(),
            .groups = "drop"
          )
        
        print(summary_stats, n = 100)
        
        # Plot 1: Boxplots of convergence (converged fits only)
        p_cons <- ggplot(converged_results, aes(x = factor(time_window), y = estimate)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightblue") +
          geom_jitter(width = 0.2, alpha = 0.3) +
          geom_hline(aes(yintercept = true_value), color = "red", linetype = "dashed", linewidth = 1) +
          facet_wrap(~param, scales = "free_y") +
          labs(title = "Parameter Consistency vs Time Window (T) - BA model",
               subtitle = "Red dashed line indicates true parameter value",
               x = "Time Window Length (T)",
               y = "Parameter Estimate") +
          theme_minimal()
        
        print(p_cons)
        
        # Plot 1b: Boxplots of convergence vs Event Count
        # Bin event counts into deciles for boxplotting
        converged_results$event_bin <- cut(converged_results$event_count, breaks = 10)
        p_cons_ev <- ggplot(converged_results, aes(x = event_bin, y = estimate)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightgreen") +
          geom_jitter(width = 0.2, alpha = 0.3) +
          geom_hline(aes(yintercept = true_value), color = "red", linetype = "dashed", linewidth = 1) +
          facet_wrap(~param, scales = "free_y") +
          labs(title = "Parameter Consistency vs Event Count - BA model",
               subtitle = "Red dashed line indicates true parameter value",
               x = "Event Count (Binned)",
               y = "Parameter Estimate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(p_cons_ev)
        
        # Plot 2: RMSE decay (The "Getting Better" plot)
        p_rmse <- ggplot(summary_stats, aes(x = time_window, y = rmse)) +
          geom_line(linewidth = 1) +
          geom_point(size = 3) +
          facet_wrap(~param, scales = "free_y") +
          labs(title = "RMSE Decay as Data Increases (BA model)",
               x = "Time Window Length (T)",
               y = "Root Mean Squared Error") +
          theme_bw()
        
        print(p_rmse)

        # Plot 2b: RMSE decay vs Mean Event Count
        p_rmse_ev <- ggplot(summary_stats, aes(x = mean_events, y = rmse)) +
          geom_line(linewidth = 1, color = "darkgreen") +
          geom_point(size = 3, color = "darkgreen") +
          facet_wrap(~param, scales = "free_y") +
          labs(title = "RMSE Decay vs Mean Event Count (BA model)",
               x = "Mean Event Count",
               y = "Root Mean Squared Error") +
          theme_bw()
        
        print(p_rmse_ev)
      } else {
        cat("  WARNING: No converged fits. Skipping consistency plots.\n")
        summary_stats <- NULL
        p_cons <- NULL
        p_cons_ev <- NULL
        p_rmse <- NULL
        p_rmse_ev <- NULL
      }
    } else {
      cat("  WARNING: consistency_results is empty. Skipping consistency plots.\n")
      summary_stats <- NULL
      p_cons <- NULL
      p_cons_ev <- NULL
      p_rmse <- NULL
      p_rmse_ev <- NULL
    }
  }

  # Save checkpoint: main + consistency results (before explosive, which may fail)
  save_list_ckpt <- list()
  if (exists("sims") && !is.null(sims)) {
    save_list_ckpt$sims <- sims
    save_list_ckpt$fits <- fits
    save_list_ckpt$temp_hawkes_fits <- temp_hawkes_fits
    save_list_ckpt$params <- params
    save_list_ckpt$params_init <- params_init
    save_list_ckpt$TIME <- TIME
    save_list_ckpt$N_SIMS <- N_SIMS
  }
  if (exists("consistency_results")) {
    save_list_ckpt$consistency_results <- consistency_results
    save_list_ckpt$summary_stats <- if (exists("summary_stats")) summary_stats else NULL
    save_list_ckpt$p_cons <- if (exists("p_cons")) p_cons else NULL
    save_list_ckpt$p_cons_ev <- if (exists("p_cons_ev")) p_cons_ev else NULL
    save_list_ckpt$p_rmse <- if (exists("p_rmse")) p_rmse else NULL
    save_list_ckpt$p_rmse_ev <- if (exists("p_rmse_ev")) p_rmse_ev else NULL
    save_list_ckpt$N_SIMS_CONSISTENCY <- N_SIMS_CONSISTENCY
    save_list_ckpt$time_windows <- time_windows
    if (exists("consistency_example_sims")) save_list_ckpt$consistency_example_sims <- consistency_example_sims
  }
  tryCatch({
    saveRDS(save_list_ckpt, file.path(CLUSTER_OUTPUT_DIR, "results_BA_full.RDS"))
    cat("Checkpoint saved (main + consistency) to results_BA_full.RDS\n")
  }, error = function(e) message("Checkpoint save failed: ", conditionMessage(e)))

# ==============================================================================
# STUDY 2: Explosive Regime Analysis - BA model
# ==============================================================================
# Goal: Set beta_overall and beta_edges -> 0.
# 1. beta_overall -> 0 means the memory of the process never decays.
#    If K > 0, the integral of intensity diverges (Explosive / Super-critical).
# 2. beta_edges -> 0 means the preferential attachment logic considers ALL past
#    nodes equally (no time decay on degree relevance).
# ==============================================================================

# --- Cleanup before explosive study ---
gc()

if(RUN_EXPLOSIVE){
  tryCatch({
  # Define Explosive Parameters - BA model
  params_explosive <- list(
    mu = 10,
    beta_overall = 0.1,
    K = 0.99,
    beta_edges = 0.1,
    m = 1
  )

  # Compare with Stable Parameters - BA model
  params_stable <- list(
    mu = 10,
    beta_overall = 2.0,
    K = 0.5,
    beta_edges = 1.0,
    m = 1
  )

  cat("Simulating Explosive Regime (BA model)...\n")

  T_explode <- 5

  # Simulate Explosive (may fail with "full networks" if params too extreme)
  sim_exp <- sim_hawkesNet(params = params_explosive,
                                 time_window = c(0, T_explode),
                                 PMF_mark = PMF_mark_BA,
                                 cond_intensity = cond_intensity,
                                 hashed_edges = TRUE,
                                 verbose = FALSE,
                                 truncation = TRUNCATION,
                                 mu_multiplier = 50)

  cat("Simulating Stable Regime (BA model)...\n")
  sim_stable <- sim_hawkesNet(params = params_stable,
                                    time_window = c(0, T_explode),
                                    PMF_mark = PMF_mark_BA,
                                    cond_intensity = cond_intensity,
                                    hashed_edges = TRUE,
                                    verbose = FALSE,
                                    truncation = TRUNCATION,
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
    labs(title = "Explosive vs Stable Process Dynamics (BA model)",
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
    labs(title = "Node Growth Over Time: Explosive vs Stable (BA)",
         x = "Time",
         y = "Number of Nodes") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p_node_growth)

  # ==========================
  # Visualization: Network Structure Impact
  # ==========================
  op <- par(mfrow = c(1, 2))

  plot(sim_stable$net, main = "Stable Network (BA)\n(Recent Activity Matters)",
       vertex.cex = 0.5, edge.col = "gray")

  plot(sim_exp$net, main = "Explosive/Memory Network (BA)\n(History Never Dies)",
       vertex.cex = 0.5, edge.col = "gray")

  par(op)

  max_deg_stable <- max(degree(sim_stable$net))
  max_deg_exp <- max(degree(sim_exp$net))

  cat("Max Degree Stable:", max_deg_stable, "\n")
  cat("Max Degree Explosive:", max_deg_exp, "\n")
  }, error = function(e) {
    message("Explosive regime simulation failed: ", conditionMessage(e))
    message("Skipping explosive block; main and consistency results will still be saved.")
  })
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
  save_list$p_cons_ev <- p_cons_ev
  save_list$p_rmse <- p_rmse
  save_list$p_rmse_ev <- p_rmse_ev
  save_list$N_SIMS_CONSISTENCY <- N_SIMS_CONSISTENCY
  save_list$time_windows <- time_windows
  if (exists("consistency_example_sims")) save_list$consistency_example_sims <- consistency_example_sims
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
# Ensure output dir exists (getwd() can become invalid on some clusters)
out_dir <- CLUSTER_OUTPUT_DIR
out_file <- file.path(out_dir, "results_BA_full.RDS")
ok <- tryCatch({
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(save_list, out_file)
  TRUE
}, error = function(e) {
  message("Save to ", out_file, " failed: ", conditionMessage(e))
  FALSE
})
if (ok) {
  print(paste("Saved full state to", out_file))
} else {
  fallback <- "results_BA_full.RDS"
  tryCatch({
    saveRDS(save_list, fallback)
    CLUSTER_OUTPUT_DIR <- "."  # so PAPER_OUTPUT can find the file
    print(paste("Saved full state to fallback:", normalizePath(fallback, mustWork = FALSE)))
  }, error = function(e) message("Fallback save also failed: ", conditionMessage(e)))
}

# ==============================================================================
# PAPER OUTPUT (at end: main study + consistency + explosive)
# Re-hydrate from results_BA_full.RDS and produce all figures/tables.
# ==============================================================================
if(PAPER_OUTPUT){
  dat <- readRDS(file.path(CLUSTER_OUTPUT_DIR, "results_BA_full.RDS"))
  list2env(dat, envir = .GlobalEnv)

  # ---------- Main study output ----------
  if(!is.null(dat$sims)){
    net_stats <- do.call(rbind, lapply(sims, function(s){
      net <- s$net
      degs <- calculateStatistics(net ~ degree(0:20,"in"))
      esps <- calculateStatistics(net ~ esp(0:20))
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

    # --- Main study summary ---
    keep_idx <- which(sapply(fits, function(x) {
      is.list(x) && !is.null(x$fit) && x$fit$convergence == 0 && 
      all(is.finite(x$fit$par)) && !any(x$fit$par > 100)
    }))
    # Additional check for beta_overall (usually index 2, but let's be safe)
    if (length(keep_idx) > 0) {
      keep_idx <- keep_idx[sapply(fits[keep_idx], function(x) {
        idx <- which(names(x$fit$par) == "beta_overall")
        if (length(idx) > 0) x$fit$par[idx] <= 10 else TRUE
      })]
    }
    
    print(paste0("keeping ", length(keep_idx), " of ", length(fits), " fits"))
    
    if (length(keep_idx) > 0) {
      # Extract all fitted parameters into a data frame
      estims <- do.call(rbind, lapply(keep_idx, function(idx) {
        as.data.frame(t(fits[[idx]]$fit$par))
      }))
      
      # Create a named vector for true and init values that matches the fit_table names
      sample_fit <- fits[[keep_idx[1]]]
      par_names <- names(sample_fit$fit$par)
      
      # Map true and init values (only for params present in fit$par)
      true_vec <- setNames(numeric(length(par_names)), par_names)
      init_vec <- setNames(numeric(length(par_names)), par_names)
      safe_set <- function(vec, nm, val) { if (nm %in% par_names) vec[nm] <- val; vec }
      true_vec <- safe_set(true_vec, "mu", params$mu)
      true_vec <- safe_set(true_vec, "beta_overall", params$beta_overall)
      true_vec <- safe_set(true_vec, "K", params$K)
      true_vec <- safe_set(true_vec, "beta_edges", params$beta_edges)
      true_vec <- safe_set(true_vec, "m", params$m)
      
      # Map init values for scalar params
      init_vec <- safe_set(init_vec, "mu", params_init$mu)
      init_vec <- safe_set(init_vec, "beta_overall", params_init$beta_overall)
      init_vec <- safe_set(init_vec, "K", params_init$K)
      init_vec <- safe_set(init_vec, "beta_edges", params_init$beta_edges)
      init_vec <- safe_set(init_vec, "m", params_init$m)
      
      results <- data.frame(
        mean = colMeans(estims),
        sd   = apply(estims, 2, sd),
        true = true_vec[colnames(estims)],
        init = init_vec[colnames(estims)]
      )
      print(results)
      
      # Plotting
      estim_long <- pivot_longer(estims, cols = everything(), names_to = "param", values_to = "estimate")
      ref_lines <- data.frame(
        param = par_names,
        true  = as.numeric(true_vec),
        init  = as.numeric(init_vec)
      )
      
      p_est_dist <- ggplot(estim_long, aes(x = estimate)) +
        geom_histogram(bins = 20, fill = "lightblue", alpha = 0.7) +
        geom_vline(data = ref_lines, aes(xintercept = true), color = "red", linetype = "dashed", linewidth = 1) +
        geom_vline(data = ref_lines, aes(xintercept = init), color = "darkgreen", linetype = "dotted", linewidth = 1) +
        facet_wrap(~param, scales = "free") +
        labs(title = "Distribution of parameter estimates (BA)", 
             subtitle = "Red dashed = true; green dotted = init", x = "Estimate") + 
        theme_minimal()
      print(p_est_dist)
    }

    comps <- lapply(sims, function(x) compensators_hawkesNet(params = params, mark_filtration = x$net, time_window = c(0, TIME)))
    if (length(keep_idx) > 0) {
      marked_p_vals <- sapply(keep_idx, function(sim_idx){
        sim <- sims[[sim_idx]]
        fit_obj <- fits[[sim_idx]]
        # Use full params (includes fixed mu, K) — fit$par alone lacks them
        pfull <- if (!is.null(fit_obj$params)) fit_obj$params else merge_fit_params(fit_obj$fit$par, params_init, fit_obj$fixed_params)
        ks_test_pval_hawkesNet(params = pfull, time_window = c(0, TIME), mark_filtration = sim$net)
      })
      cat("Mean marked KS p-value:", mean(marked_p_vals), "\n")
    }
    temp_p_vals <- sapply(seq_along(sims), function(i){
      ks_test_pval_temporal(realiz = data.frame(t = sims[[i]]$events$t, n = rep(sims[[i]]$events$n, length(sims[[i]]$events$t))), windowT = c(0, TIME), hawkes_par = temp_hawkes_fits[[i]]$par)
    })
    cat("Mean temporal KS p-value:", mean(temp_p_vals), "\n")
  }

  # ---------- Consistency study output ----------
  if(!is.null(dat$consistency_results) && nrow(dat$consistency_results) > 0){
    prop_keep <- consistency_results %>% group_by(time_window, param) %>% summarise(prop_keep = sum(keep) / N_SIMS_CONSISTENCY, .groups = "drop")
    print(prop_keep)
    if (!is.null(summary_stats)) print(summary_stats)
    if (!is.null(p_cons)) print(p_cons)
    if (!is.null(p_cons_ev)) print(p_cons_ev)
    if (!is.null(p_rmse)) print(p_rmse)
    if (!is.null(p_rmse_ev)) print(p_rmse_ev)
  }

  # ---------- Explosive study output ----------
  if(!is.null(dat$sim_exp)){
    print(p_expl)
    if(!is.null(dat$p_node_growth)) print(p_node_growth)
    op <- par(mfrow = c(1, 2))
    plot(sim_stable$net, main = "Stable Network (BA)\n(Recent Activity Matters)", vertex.cex = 0.5, edge.col = "gray")
    plot(sim_exp$net, main = "Explosive/Memory Network (BA)\n(History Never Dies)", vertex.cex = 0.5, edge.col = "gray")
    par(op)
    print(paste("Max Degree Stable:", max_deg_stable))
    print(paste("Max Degree Explosive:", max_deg_exp))
  }

  # ---------- Re-save with PAPER_OUTPUT objects (tables + plots) ----------
  # Add main study results table and estimate distribution plot
  if (exists("results") && !is.null(results))    save_list$results_table   <- results
  if (exists("estims") && !is.null(estims))      save_list$estims          <- estims
  if (exists("p_est_dist") && !is.null(p_est_dist)) save_list$p_est_dist   <- p_est_dist
  if (exists("deg_plot") && !is.null(deg_plot))   save_list$deg_plot        <- deg_plot
  if (exists("esp_plot") && !is.null(esp_plot))   save_list$esp_plot        <- esp_plot
  if (exists("keep_idx"))                         save_list$keep_idx        <- keep_idx
  if (exists("marked_p_vals"))                    save_list$marked_p_vals   <- marked_p_vals
  if (exists("temp_p_vals"))                      save_list$temp_p_vals     <- temp_p_vals
  # Consistency plots (already in save_list from earlier, but re-add for safety)
  if (exists("p_cons") && !is.null(p_cons))       save_list$p_cons          <- p_cons
  if (exists("p_cons_ev") && !is.null(p_cons_ev)) save_list$p_cons_ev       <- p_cons_ev
  if (exists("p_rmse") && !is.null(p_rmse))       save_list$p_rmse          <- p_rmse
  if (exists("p_rmse_ev") && !is.null(p_rmse_ev)) save_list$p_rmse_ev       <- p_rmse_ev
  # Re-save
  tryCatch({
    saveRDS(save_list, file.path(CLUSTER_OUTPUT_DIR, "results_BA_full.RDS"))
    cat("Re-saved full state (with tables + plots) to results_BA_full.RDS\n")
  }, error = function(e) message("Re-save failed: ", conditionMessage(e)))
}
