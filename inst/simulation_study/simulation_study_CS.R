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
CLUSTER_OUTPUT_DIR <- file.path(PKG_ROOT, "cluster_output")
dir.create(CLUSTER_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===================================================
# Change Statistic Mark Generation
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
# takes ~ 6 minutes to simulate at T=50, gives ~600 events
params <- list(mu = 10,
               beta_overall = 2,
               K = 0.5,
               beta_edges = 1,
               node_lambda = 1,
               m = 1,
               CS_params = c(-7, 3, 0.1, -0.1)
               )
TRUNCATION  <- 500
SIMULATE <- TRUE
PAPER_OUTPUT <- TRUE
RUN_EXPLOSIVE <- ON_SLURM          # skip explosive in interactive mode
RUN_CONSISTENCY <- ON_SLURM  # Run on cluster batch jobs only (slow)
MAX_ITER <- if (ON_SLURM) 5000 else 2000
cat(sprintf("Mode: %s | TIME=%d | N_SIMS=%d | N_CORES=%d | MAX_ITER=%d\n",
    if (ON_SLURM) "SLURM batch" else "Interactive (RStudio)", TIME, N_SIMS, N_CORES, MAX_ITER))

if (nzchar(Sys.getenv("CORES_OVERRIDE"))) {
  N_CORES <- as.numeric(Sys.getenv("CORES_OVERRIDE"))
} else if (N_CORES == 128L && nzchar(Sys.getenv("USE_256_WHEN_128"))) {
  # Cluster gave 256 when you requested 128 (squeue shows 256); use them.
  N_CORES <- 256L
  cat("Request was 128; using 256 (USE_256_WHEN_128 set; typical when squeue shows 256)\n")
}
# Core allocation: 128 -> 16 outer x 8 inner; 256 -> 16 outer x 16 inner.
# We prioritize inner cores now because the intensity cache is the bottleneck.
# R < 4.4.0 has a socket limit of 128. We cap outer workers at 64 to be safe.
CORES_OUTER_ENV <- Sys.getenv("CORES_OUTER", "")
if (nzchar(CORES_OUTER_ENV)) {
  N_CORES_OUTER <- as.numeric(CORES_OUTER_ENV)
  N_CORES_INNER <- max(1L, floor(N_CORES / N_CORES_OUTER))
} else if (N_CORES >= 128L) {
  # Use fewer outer workers but many more inner cores for intensity cache speedup
  N_CORES_OUTER <- 16L 
  N_CORES_INNER <- max(1L, floor(N_CORES / N_CORES_OUTER))
  cat("Using high-core mode:", N_CORES_OUTER, "outer workers x", N_CORES_INNER, "inner =", N_CORES_OUTER * N_CORES_INNER, "total\n")
} else {
  N_CORES_INNER <- as.numeric(Sys.getenv("CORES_INNER", 8)) # Default to 8 inner
  N_CORES_OUTER <- max(1L, floor(N_CORES / N_CORES_INNER))
}
# Final safety cap for R socket limits (128 total)
# Each PSOCK worker with outfile="" uses 2 connections.
N_CORES_OUTER <- min(N_CORES_OUTER, 60L)
N_CORES_INNER <- max(1L, floor(N_CORES / N_CORES_OUTER))

cat("Core allocation:", N_CORES_OUTER, "outer x", N_CORES_INNER, "inner =",
    N_CORES_OUTER * N_CORES_INNER, "total (of", N_CORES, "available)\n")

SEED <- 1267

  # parscale: match param magnitudes so Nelder-Mead simplex steps are proportionate
  # Defined here (not inside SIMULATE block) so consistency study can also use it.
  # IMPORTANT: Names must match flat_par names exactly.
  # fixed_params = c("K") -> flat_par has:
  #   mu, beta_overall, beta_edges, node_lambda, m, CS_params1, CS_params2, CS_params3, CS_params4
  # K is fixed as in the BA model (where K=1 usually).
  # mu, beta_overall, beta_edges, node_lambda, and m are now all free.
  p_scale <- c(mu = 1, beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.5, m = 0.1,
               CS_params1 = 1, CS_params2 = 0.1, CS_params3 = 0.1, CS_params4 = 0.1)

make_cluster <- function(n_workers) {
  # PSOCK cluster: each worker runs one sim or one fit at a time.
  # outfile="" ensures worker stdout/stderr are forwarded to the master (visible in .out/.err).
  cl <- makeCluster(n_workers, outfile = "")
  registerDoParallel(cl)
  # export libraries to cluster:
  clusterEvalQ(cl, {
    library(spatstat)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(data.table)
    library(pbapply)
    library(parallel)
    library(doParallel)
    library(R.utils)
    library(ernm)
    library(network)
    library(sna)
    library(hash)
    library(hawkesNet)
    # Ensure plogis is available (needed by PMF_mark_CS closures when serialized to inner workers)
    plogis <- stats::plogis
    # CRITICAL: Pre-set BLAS/OpenMP threads to 1 in each PSOCK worker.
    # When fit_hawkesNet uses inner parallelism (mclapply fork), forked
    # grandchildren inherit the worker's thread state. If BLAS is multi-threaded,
    # the fork deadlocks. Setting threads=1 here prevents that.
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    }
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
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
  # 1. Simulation Step: Use ALL available cores for outer workers (capped at 120 for R socket limit)
  # Simulation is single-threaded, so nested parallelism is not needed here.
  t <- proc.time()
  # R < 4.4.0 has a limit of 128 total connections. PSOCK workers use 1 each.
  # We cap at 120 to leave room for files/stdout/etc.
  # NeSI scheduling/memory: too many PSOCK workers can bloat memory.
  # Default to a smaller simulation worker pool; override with SIM_WORKERS.
  N_SIM_WORKERS <- as.numeric(Sys.getenv("SIM_WORKERS", min(N_CORES, 16L)))
  N_SIM_WORKERS <- max(1L, min(N_SIM_WORKERS, N_CORES))
  cat("Commencing simulation using", N_SIM_WORKERS, "parallel workers...\n")
  cl_sim <- makeCluster(N_SIM_WORKERS)
  registerDoParallel(cl_sim)
  clusterEvalQ(cl_sim, {
    library(hawkesNet)
    library(ernm)
    library(network)
    library(data.table)
  })
  clusterExport(cl_sim, c("params", "TIME", "TRUNCATION", "SEED"))
  
  set.seed(SEED)
  sims <- parLapply(cl=cl_sim, 1:N_SIMS, function(x){
    tryCatch({
      sim_hawkesNet(params = params,
                    time_window = c(0, TIME),
                    PMF_mark = PMF_mark_CS,
                    cond_intensity = cond_intensity,
                    hashed_edges = TRUE,
                    verbose = FALSE,
                    mu_multiplier = 3,
                    joint_accept = FALSE,
                    truncation = TRUNCATION,
                    formula_RHS = "edges + triangles + star(c(2,3))",
                    mark_decay = "node_entrance",
                    growth_only = FALSE)
    }, error = function(e) {
      message("Error in sim_hawkesNet: ", e$message)
      return(NULL)
    })
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
  
  # Init values: node_lambda fixed to true value; mu and K free.
  # node_lambda fixed to prevent co-variation with edges term.
  # mu and K are identifiable from event times (ground intensity only).
  params_init <- list(mu = 10,
                      beta_overall = 1,
                      K = 0.5,
                      beta_edges = 1,
                      node_lambda = 1,
                      m = 1,
                      CS_params = c(-9, 2, 0, 0)
  )
  clusterExport(cl_fit, c("params_init", "p_scale", "TIME", "MAX_ITER", "TRUNCATION"))
  
  t1 <- proc.time()
  fits <- parLapply(cl=cl_fit, sims, function(x){
    worker_id <- Sys.getpid()
    message(sprintf("  [Outer Worker %d] Starting fit for sim with %d events...", worker_id, length(x$events$t)))
    
      fit <- tryCatch({
      fit_hawkesNet(
        params_init = params_init,
        time_window = c(0, TIME),
        mark_filtration = x$net,
        PMF_mark = PMF_mark_CS,
        formula_RHS = "edges + triangles + star(c(2,3))",
        trace = 1,
        maxit = MAX_ITER,
        truncation = TRUNCATION,
        mark_decay = "node_entrance",
        growth_only = FALSE,
        fixed_params = c("K"),
        method = "Nelder-Mead",
        parscale = p_scale,
        cores = N_CORES_INNER,
        cache_intensity = TRUE,
        combine_intensity = TRUE,
        verbose = TRUE
      )
    }, error = function(e) {
      message(sprintf("  [Outer Worker %d] ERROR: %s", worker_id, e$message))
      return(e$message)
    })
    
    if (is.list(fit)) {
      message(sprintf("  [Outer Worker %d] Fit complete (convergence=%d, iterations=%d)", 
                      worker_id, fit$fit$convergence, fit$fit$counts[1]))
      if (!is.null(fit$intens_funcs)) fit$intens_funcs <- NULL
    }
    return(fit)
  })
  
  cat("Fitting took:", round((proc.time() - t1)[3], 1), "s\n")
  
  # 3. Temporal Hawkes: Use all cores again
  cat("\nFitting Temporal Hawkes models...\n")
  # (Reuse cl_fit logic or just use parLapply on all cores)
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
          file = file.path(CLUSTER_OUTPUT_DIR, "results_CS.RDS"))
  cat("Total SIMULATE block took:", round((proc.time() - t)[3], 1), "s\n")
  gc()
}

# ==============================================================================
# STUDY 1: Consistency Analysis (Sliding Window / Increasing T) - CS model
# ==============================================================================
# Goal: Show that as TIME increases, the variance of estimates decreases and
# means converge to truth.
# ==============================================================================

if(RUN_CONSISTENCY){

  # 1. Define Time Windows to test (include T=20 for consistency with diagnostic)
  time_windows <- c(5, 10, 20, 25, 50, 75, 100)
  N_SIMS_CONSISTENCY <- 25

  # Parameters (Standard/Stable regime) - CS model
  params_true <- list(mu = 10,
                      beta_overall = 2,
                      K = 0.5,
                      beta_edges = 1,
                      node_lambda = 1,
                      CS_params = c(-7, 3, 0.1, -0.1))

    # Setup Cluster — nested parallelism (PSOCK outer x fork inner).
    # BLAS threads pre-set to 1 in make_cluster() so forked grandchildren are safe.
    t_consistency_total <- proc.time()
    cat("=== Consistency Study (CS) ===\n")
    cat("  Time windows:", paste(time_windows, collapse = ", "), "\n")
    cat("  N_SIMS per window:", N_SIMS_CONSISTENCY, "\n")
    cat("  Core allocation:", N_CORES_OUTER, "outer x", N_CORES_INNER, "inner =",
        N_CORES_OUTER * N_CORES_INNER, "total (of", N_CORES, "available)\n")
    cat("Setting up cluster...\n")
    t_cluster <- proc.time()
    cl <- make_cluster(N_CORES_OUTER)
    clusterExport(cl, c("params_true", "TRUNCATION", "N_CORES_INNER", "MAX_ITER", "p_scale"))
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
          sim_hawkesNet(params = params_true,
                              time_window = c(0, curr_time),
                              PMF_mark = PMF_mark_CS,
                              cond_intensity = cond_intensity,
                              hashed_edges = TRUE,
                              mu_multiplier = 3,
                              verbose = FALSE,
                              truncation = TRUNCATION,
                              formula_RHS = "edges + triangles + star(c(2,3))",
                              mark_decay = "node_entrance",
                              growth_only = FALSE)
        }, error = function(e) return(NULL))
        
        if(is.null(sim_res)) return(NULL)
        
        # B. Fit
        # Initialize near true params + small noise for better convergence.
        # node_lambda fixed to true value; mu and K free with noise.
        params_init <- list(
          mu = max(0.1, params_true$mu * exp(rnorm(1, 0, 0.2))),
          beta_overall = max(0.1, params_true$beta_overall * exp(rnorm(1, 0, 0.2))),
          K = min(0.99, max(0.01, params_true$K * exp(rnorm(1, 0, 0.2)))),
          beta_edges = max(0.1, params_true$beta_edges * exp(rnorm(1, 0, 0.2))),
          node_lambda = params_true$node_lambda,
          m = max(0.1, params_true$m * exp(rnorm(1, 0, 0.2))),
          CS_params = params_true$CS_params + rnorm(length(params_true$CS_params), 0, 0.5)
        )
        # Ensure CS_params are finite
        params_init$CS_params[!is.finite(params_init$CS_params)] <- params_true$CS_params[!is.finite(params_init$CS_params)]
        
        fit_res <- tryCatch({
          fit_hawkesNet(params_init = params_init,
                              time_window = c(0, curr_time),
                              mark_filtration = sim_res$net,
                              PMF_mark = PMF_mark_CS,
                              formula_RHS = "edges + triangles + star(c(2,3))",
                              maxit = MAX_ITER,
                              truncation = TRUNCATION,
                              mark_decay = "node_entrance",
                              growth_only = FALSE,
                              cache_intensity = TRUE,
                              combine_intensity = TRUE,
                              verbose = FALSE,
                              fixed_params = c("K"),
                              parscale = p_scale,
                              cores = N_CORES_INNER,
                              method = "Nelder-Mead")
      }, error = function(e) return(NULL))
        
        if(is.null(fit_res) || is.null(fit_res$fit)) return(NULL)
        
        # Check if fit succeeded
        par_bo_idx <- which(names(fit_res$fit$par) == "beta_overall")
        keep <- (!is.null(fit_res$fit) && 
                 length(fit_res$fit) > 0 && 
                 fit_res$fit$convergence == 0 &&
                 all(is.finite(fit_res$fit$par)) &&
                 !any(fit_res$fit$par > 100) && 
                 length(fit_res$fit$par) >= 2 &&
                 (length(par_bo_idx) == 0 || fit_res$fit$par[par_bo_idx] <= 10))
        
        # Map true values: unlist gives CS_params1,2,...; fit$par may use those or formula names
        par_names <- names(fit_res$fit$par)
        true_vals <- setNames(numeric(length(par_names)), par_names)
        if ("mu" %in% par_names) true_vals["mu"] <- params_true$mu
        true_vals["beta_overall"] <- params_true$beta_overall
        if ("K" %in% par_names) true_vals["K"] <- params_true$K
        true_vals["beta_edges"] <- params_true$beta_edges
        true_vals["node_lambda"] <- params_true$node_lambda
        if ("m" %in% par_names) true_vals["m"] <- params_true$m
        for (k in seq_along(params_true$CS_params)) {
          nm <- paste0("CS_params", k)
          if (nm %in% par_names) true_vals[nm] <- params_true$CS_params[k]
        }
        # Fallback: formula names from expected_params
        exp_cs <- tryCatch(expected_params_PMF_mark_CS(sim_res$net, "edges + triangles + star(c(2,3))"), error = function(e) NULL)
        if (!is.null(exp_cs$CS_params_names) && length(params_true$CS_params) >= length(exp_cs$CS_params_names)) {
          for (k in seq_along(exp_cs$CS_params_names)) {
            if (exp_cs$CS_params_names[k] %in% par_names) true_vals[exp_cs$CS_params_names[k]] <- params_true$CS_params[k]
          }
        }
        # Return row
        return(data.frame(
          keep = keep,
          sim_id = i,
          time_window = curr_time,
          param = par_names,
          estimate = as.numeric(fit_res$fit$par),
          true_value = as.numeric(true_vals)
        ))
      })
      elapsed_simfit <- (proc.time() - t_simfit)[3]
      
      # Bind results
      res_df <- do.call(rbind, res_list)
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
    cat("\n=== Consistency Study (CS) complete ===\n")
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
        # Calculate Bias and RMSE (converged fits only - non-converged inflate RMSE and break monotonicity)
        summary_stats <- converged_results %>%
          group_by(time_window, param) %>%
          summarise(
            mean_est = mean(estimate),
            sd_est = sd(estimate),
            rmse = sqrt(mean((estimate - true_value)^2)),
            true_val = mean(true_value),
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
          labs(title = "Parameter Consistency vs Time Window (T) - CS model",
               subtitle = "Red dashed line indicates true parameter value",
               x = "Time Window Length (T)",
               y = "Parameter Estimate") +
          theme_minimal()
        
        print(p_cons)
        
        # Plot 2: RMSE decay (The "Getting Better" plot)
        p_rmse <- ggplot(summary_stats, aes(x = time_window, y = rmse)) +
          geom_line(linewidth = 1) +
          geom_point(size = 3) +
          facet_wrap(~param, scales = "free_y") +
          labs(title = "RMSE Decay as Data Increases (CS model)",
               x = "Time Window Length (T)",
               y = "Root Mean Squared Error") +
          theme_bw()
        
        print(p_rmse)
      } else {
        cat("  WARNING: No converged fits. Skipping consistency plots.\n")
        summary_stats <- NULL
        p_cons <- NULL
        p_rmse <- NULL
      }
    } else {
      cat("  WARNING: consistency_results is empty. Skipping consistency plots.\n")
      summary_stats <- NULL
      p_cons <- NULL
      p_rmse <- NULL
    }
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

# --- Cleanup before explosive study ---
gc()

if(RUN_EXPLOSIVE){

  # Define Explosive Parameters - CS model (include node_lambda, CS_params)
  params_explosive <- list(
    mu = 10,
    beta_overall = 0.1,
    K = 0.99,
    beta_edges = 0.1,
    node_lambda = 1,
    CS_params = c(-7, 3, 0.1, -0.1)
  )

  # Compare with Stable Parameters - CS model
  params_stable <- list(
    mu = 10,
    beta_overall = 2.0,
    K = 0.5,
    beta_edges = 1.0,
    node_lambda = 1,
    CS_params = c(-7, 3, 0.1, -0.1)
  )

  cat("Simulating Explosive Regime (CS model)...\n")

  T_explode <- 5

  # Simulate Explosive
  sim_exp <- sim_hawkesNet(params = params_explosive,
                                 time_window = c(0, T_explode),
                                 PMF_mark = PMF_mark_CS,
                                 cond_intensity = cond_intensity,
                                 hashed_edges = TRUE,
                                 verbose = FALSE,
                                 truncation = TRUNCATION,
                                 formula_RHS = "edges + triangles + star(c(2,3))",
                                 mark_decay = "node_entrance",
                                 growth_only = FALSE,
                                 mu_multiplier = 50)

  cat("Simulating Stable Regime (CS model)...\n")
  sim_stable <- sim_hawkesNet(params = params_stable,
                                    time_window = c(0, T_explode),
                                    PMF_mark = PMF_mark_CS,
                                    cond_intensity = cond_intensity,
                                    hashed_edges = TRUE,
                                    verbose = FALSE,
                                    truncation = TRUNCATION,
                                    formula_RHS = "edges + triangles + star(c(2,3))",
                                    mark_decay = "node_entrance",
                                    growth_only = FALSE,
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

  cat("Max Degree Stable:", max_deg_stable, "\n")
  cat("Max Degree Explosive:", max_deg_exp, "\n")
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
# Ensure output dir exists (getwd() can become invalid on some clusters)
out_dir <- file.path(PKG_ROOT, "cluster_output")
out_file <- file.path(out_dir, "results_CS_full.RDS")
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
  fallback <- "results_CS_full.RDS"
  tryCatch({
    saveRDS(save_list, fallback)
    CLUSTER_OUTPUT_DIR <- "."  # so PAPER_OUTPUT can find the file
    print(paste("Saved full state to fallback:", normalizePath(fallback, mustWork = FALSE)))
  }, error = function(e) message("Fallback save also failed: ", conditionMessage(e)))
}

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
      # (e.g. "edges" instead of "CS_params1")
      # We use the first successful fit to get the mapping
      sample_fit <- fits[[keep_idx[1]]]
      par_names <- names(sample_fit$fit$par)
      
      # Map true and init values
      true_vec <- setNames(numeric(length(par_names)), par_names)
      init_vec <- setNames(numeric(length(par_names)), par_names)
      if ("mu" %in% par_names) true_vec["mu"] <- params$mu
      true_vec["beta_overall"] <- params$beta_overall
      if ("K" %in% par_names) true_vec["K"] <- params$K
      true_vec["beta_edges"] <- params$beta_edges
      true_vec["node_lambda"] <- params$node_lambda
      if ("m" %in% par_names) true_vec["m"] <- params$m
      # Map CS params: fit may use CS_params1,2,... or formula names (edges, triangles, star.2, star.3)
      exp_cs <- expected_params_PMF_mark_CS(sims[[1]]$net, "edges + triangles + star(c(2,3))")
      for (i in seq_along(params$CS_params)) {
        nm <- paste0("CS_params", i)
        if (nm %in% par_names) {
          true_vec[nm] <- params$CS_params[i]
          init_vec[nm] <- params_init$CS_params[i]
        }
      }
      if (!is.null(exp_cs$CS_params_names)) {
        for (i in seq_along(exp_cs$CS_params_names)) {
          name <- exp_cs$CS_params_names[i]
          if (name %in% par_names) {
            true_vec[name] <- params$CS_params[i]
            init_vec[name] <- params_init$CS_params[i]
          }
        }
      }
      
      # Map init values for scalar params (CS done above)
      if ("mu" %in% par_names) init_vec["mu"] <- params_init$mu
      init_vec["beta_overall"] <- params_init$beta_overall
      if ("K" %in% par_names) init_vec["K"] <- params_init$K
      init_vec["beta_edges"] <- params_init$beta_edges
      init_vec["node_lambda"] <- params_init$node_lambda
      if ("m" %in% par_names) init_vec["m"] <- params_init$m
      
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
        labs(title = "Distribution of parameter estimates (CS)", 
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
    prop_keep <- consistency_results %>% group_by(time_window, param) %>% summarise(prop_keep = sum(keep) / N_SIMS_CONSISTENCY)
    print(prop_keep)
    if (!is.null(summary_stats)) print(summary_stats)
    if (!is.null(p_cons)) print(p_cons)
    if (!is.null(p_rmse)) print(p_rmse)
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
  if (exists("p_rmse") && !is.null(p_rmse))       save_list$p_rmse          <- p_rmse
  # Re-save
  tryCatch({
    saveRDS(save_list, file.path(CLUSTER_OUTPUT_DIR, "results_CS_full.RDS"))
    cat("Re-saved full state (with tables + plots) to results_CS_full.RDS\n")
  }, error = function(e) message("Re-save failed: ", conditionMessage(e)))
}
