# This script fits a LOLOG style model and a BA model to messaging data:

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

# ===================================================
# Change Statistic Mark Generation
# ===================================================
INVESTIGATE <- FALSE
SIMULATE <- TRUE
PAPER_OUTPUT = TRUE
RUN_EXPLOSIVE <- TRUE
RUN_CONSISTENCY <- TRUE

TIME <- 10
params <- list(mu = 10,
               beta_overall = 1,
               K = 0.5,
               beta_edges = 1
)
TRUNCATION  = 100


DEBUG = FALSE
MAX_ITER = 2000

N_SIMS = 21
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 7))

SEED <- 01267

make_cluster <- function(N_CORES){
  # setup the cluster:
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  estimatedparams<-list()
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
  # make the cluster:
  t <- proc.time()
  cl <- make_cluster(N_CORES)
  set.seed(SEED)
  
  # Ensure the cluster will be stopped no matter what.
  on.exit({
    if (!is.null(cl)) {
      stopCluster(cl)
    }
  }, add = TRUE)
  
  sims <- parLapply(cl=cl,1:N_SIMS,function(x){
    results <- tryCatch({
      sim_hawkesGrowthNet(params =  params,
                          time_window = c(0,TIME),
                          PMF_mark = PMF_mark_BA,
                          cond_intensity = cond_intensity,
                          hashed_edges = T,
                          verbose = F,
                          mu_multiplier = 3,
                          joint_accept = F,
                          truncation = TRUNCATION)},
      error = function(e) {
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
  params_init <- list(mu = 0.1,
                      beta_overall = 0.1,
                      beta_edges = 0.1,
                      K = params$K
                      )
  clusterExport(cl, c("params_init"))
  t1 <- proc.time()
  fits <- parLapply(cl=cl,sims,function(x){
    fit <- tryCatch({
      fit_hawkesGrowthNet(
        params_init = params_init,
        time_window = c(0, TIME),
        mark_filtration = x$net,
        PMF_mark = PMF_mark_BA,
        grad = FALSE,
        trace = 0,
        maxit = MAX_ITER,
        truncation = TRUNCATION,
        get_hessian = TRUE#,
        # fixed_params = c("K")
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
  saveRDS(list(sims=sims,
               fits = fits,
               temp_hawkes_fits = temp_hawkes_fits,
               params = params,
               params_init = params_init
               ),
          file = "results_BA.RDS")
  stopCluster(cl)
  print("Simulating and fitting took:")
  print((t - proc.time())[3])
}

if(PAPER_OUTPUT){
  
  results_BA <- readRDS("results_BA.RDS")
  sims <- results_BA$sims
  fits <- results_BA$fits
  temp_hawkes_fits <- results_BA$temp_hawkes_fits
  
  # ==========================
  # Network Descriptive Stats
  # ==========================
  # make the ESP graphs
  net_stats <- do.call(rbind,lapply(sims,function(s){
    net <- s$net
    degs <- ernm::calculateStatistics(net ~ degree(0:20,"in"))
    esps <- ernm::calculateStatistics(net ~ esp(0:20))
    
    tmp <- as.data.frame(cbind(c(degs,esps),
                               rep(0:20,times = 2),
                               c(rep("degree",21),rep("esp",21))
    ))
    rownames(tmp) <- NULL
    
    return(tmp)
  }))
  names(net_stats) <- c("value","var","type")
  net_stats$value <- as.numeric(net_stats$value)
  net_stats$var <- as.numeric(net_stats$var)
  net_stats <- net_stats[net_stats$var <=15,]
  # Ensure 'var' is ordered as a factor
  net_stats <- net_stats %>%
    mutate(var = factor(var, levels = sort(unique(var))))
  
  # make the boxplots
  deg_plot <- ggplot(net_stats[net_stats$type == "degree",],aes(x = var,y = value)) +
    geom_boxplot() +
    labs(title = "Degree Distribution",
         x = "Degree",
         y = "Value")+
    theme_minimal()
  deg_plot
  
  # esp plot:
  esp_plot <- ggplot(net_stats[net_stats$type == "esp",],aes(x = var,y = value)) +
    geom_boxplot() +
    labs(title = "ESP Distribution",
         x = "ESP",
         y = "Value")+
    theme_minimal()
  esp_plot
  
  # mean degree :
  mean_degs <- sapply(sims,function(s){
    mean(degree(s$net,gmode = "graph"))
  })
  mean_deg_df <- data.frame(mean_deg = mean_degs)
  hist(mean_deg_df$mean_deg,
       main = "Histogram of Mean Degrees",
       xlab = "Mean Degree",
       breaks = 10)
  vlines <- mean(mean_deg_df$mean_deg)
  abline(v = vlines, col = "red", lwd = 2)
  
  # ==========================
  # RESULTS TABLE
  # ==========================
  # get mean and sd of parameters from fits:
  keep <- which(sapply(fits,function(x){length(x)!=0 & x$fit$convergence==0 & !any(x$fit$par > 100) & !any(x$fit$par[2] >10)}))
  paste0("keeping ",length(keep), " of ", N_SIMS," fits")
  
  estims <- do.call(rbind,lapply(fits[keep],function(x){
    return(as.data.frame(t(x$fit$par),names = names(x$fit$par)))
  }))
  
  params_vec <- unlist(params)
  params_vec <- params_vec[names(params_vec) %in% colnames(estims)]
  params_init_vec <- unlist(params_init)
  params_init_vec <- params_init_vec[names(params_init_vec) %in% colnames(estims)]
  
  
  results <- data.frame(mean = colMeans(estims),
                        sd = apply(estims,2,sd),
                        true = params_vec,
                        init = params_init_vec)
  print(results)
  estims
  
  # ==========================
  # KS TEST Table
  # ==========================
  # With true params:
  comps <- lapply(sims,function(x){
    compensators_hawkesGrowthNet(params = params,
                                 mark_filtration = x$net,
                                 time_window = c(0,TIME)
    )
  })
  marked_p_vals <- mapply(sims[keep],lapply(1:length(fits[keep]),function(x){fits[keep][[x]]$fit$par}),FUN = function(x,y){
    times <- get_times(x$net)$times
    ks_test_pval_temporal(realiz = data.frame(t = times,
                                              n = rep(length(times),length(times))),
                          windowT = c(0,TIME),
                          hawkes_par = y
    )
  })
  
  temp_p_vals <- mapply(sims,temp_hawkes_fits,FUN = function(x,y){
    ks_test_pval_temporal(realiz = data.frame(t = x$events$t,
                                              n = rep(x$events$n,length(x$events$t))),
                          windowT = c(0,TIME),
                          hawkes_par = y$par
                          
    )
  })
  
  # need to investigate this! - expect higher pvals for true model
  # LOOK INTO PARAMETIZATION OF K !
  mean(marked_p_vals)
  mean(temp_p_vals)
}

if(INVESTIGATE){
  t <- proc.time()
  results <- sim_hawkesGrowthNet(params =  params,
                                 time_window = c(0,TIME),
                                 PMF_mark = PMF_mark_BA,
                                 cond_intensity = cond_intensity,
                                 hashed_edges = T,
                                 verbose = T,
                                 mu_multiplier = 5,
                                 joint_accept = F,
                                 truncation = TRUNCATION
  )
  print("Simulation took:")
  print(proc.time()-t)
  
  ernm::calculateStatistics(results$net ~ degree(0:15))
  ernm::calculateStatistics(results$net ~ esp(0:15))
  
  degs <- degree(results$net)
  mean(degs)
  
  plot(results$net)
  
  if(DEBUG){
    # ==================================
    # Verify Simulation is reasonable
    # ==================================
    
    # check accept probabilites:
    length(results$accept_probs)
    summary(results$accept_probs)
    length(results$events$t)
    plot(results$accept_probs)
    
    # Should be "spikey" due to hawkesian arrival times
    times <- results$net %v% 'time'
    plot(results$net,
         vertex.cex = times/TIME,
         main = '')
    
    # Set up an empty plot with appropriate x-limits and no y-axis ticks
    # plot the times on a number line
    plot(c(0,TIME), c(-1, 1), type = "n", yaxt = "n",
         xlab = "Value", ylab = "", main = "Vector on a Number Line")
    abline(h = 0, col = "gray", lwd = 2)
    points(results$events$t, rep(0, length(results$events$t)), pch = 19, col = "blue", cex = 1.5)
    
    
    # ==================================
    # Plot degrees and ESP distributions
    # ==================================
    degs <- ernm::calculateStatistics(results$net ~ degree(0:15,"in"))
    print(degs)
    plot(degs, col = 'red')
    
    plot(y=degs[2:10]/(results$net %n% 'n'),x=2:10,col = 'red')
    points((2:10)**-3,x = 2:10,col = 'blue')
    
    esps <- ernm::calculateStatistics(results$net ~ esp(0:10))
    print(esps)
    plot(esps)
    
    
    # K
    l_k <- sapply(seq(0,1,length.out=10),function(x){
      print(x)
      tmp <- params
      tmp$K <- x
      loglik_hawkesGrowthNet(params = tmp,
                             time_window = c(0,TIME),
                             mark_filtration = results$net,
                             truncation = TRUNCATION,
                             PMF_mark = PMF_mark_BA,
                             verbose = F
      )$loglik
    })
    plot(y = l_k,x = seq(0,1,length.out = 10))
    
    # beta overall
    l_b <- sapply(seq(0,10,length.out =20),function(x){
      print(x)
      tmp <- params
      tmp$beta_overall <- x
      loglik_hawkesGrowthNet(params = tmp,
                             time_window = c(0,TIME),
                             mark_filtration = results$net,
                             truncation = TRUNCATION,
                             PMF_mark = PMF_mark_BA,
                             verbose = F
      )$loglik
    })
    plot(y = l_b,x = seq(0,10,length.out = 20))
    
    # beta edges
    l_e <- sapply(seq(0,2,length.out =40),function(x){
      print(x)
      tmp <- params
      tmp$CS_params[2] <- x
      loglik_hawkesGrowthNet(params = tmp,
                             time_window = c(0,TIME),
                             mark_filtration = results$net,
                             truncation = TRUNCATION,
                             PMF_mark = PMF_mark_BA,
                             verbose = F
      )$loglik
    })
    plot(y = l_e,x = seq(0,2,length.out = 40))
    
    # K-beta plane :
    K_values <- seq(0, 1, length.out = 20)
    beta_values <- seq(0, 5, length.out = 10)
    
    # Initialize matrix to store log-likelihood values
    loglik_matrix <- matrix(NA, nrow = length(K_values), ncol = length(beta_values))
    # Compute log-likelihood over the grid
    for (i in seq_along(K_values)) {
      for (j in seq_along(beta_values)) {
        tmp <- params
        tmp$K <- K_values[i]
        tmp$beta_overall <- beta_values[j]
        
        loglik_matrix[i, j] <-   loglik_hawkesGrowthNet(params = tmp,
                                                        time_window = c(0,TIME),
                                                        mark_filtration = results$net,
                                                        truncation = TRUNCATION,
                                                        PMF_mark = PMF_mark_BA,
                                                        verbose = F)$loglik
      }
    }
    # Plot using persp (Base R 3D Plot)
    persp(K_values, beta_values, loglik_matrix,
          theta = 30, phi = 30,
          col = "lightblue", shade = 0.5,
          xlab = "K", ylab = "Beta", zlab = "Log-Likelihood",
          main = "Log-Likelihood Surface")
    
    # I think there are large regions where the likelihood is very flat:
  }
  
  # ======================================================
  # Fit the model
  # =====================================================
  
  # start at mu over time
  params_init <- list(mu = 10,
                      beta_overall = 0.1,
                      K = 0.1,
                      beta_edges = 0.1
  )
  
  t <- proc.time()
  l <- loglik_hawkesGrowthNet(params = params,
                              time_window = c(0,TIME),
                              mark_filtration = results$net,
                              PMF_mark = PMF_mark_BA,
                              edge_hash_list = NULL,
                              truncation = TRUNCATION,
                              verbose = TRUE
  )
  l$loglik
  print("1 iteration of log likelihoods took:")
  print(proc.time()-t)
  
  t <- proc.time()
  fit <- fit_hawkesGrowthNet(params_init = params_init,
                             time_window = c(0,TIME),
                             mark_filtration = results$net,
                             PMF_mark = PMF_mark_BA,
                             grad = T,
                             trace = 1,
                             truncation = TRUNCATION,
                             maxit = 1000,
                             verbose = TRUE,
                             get_hessian = TRUE
  )
  print("results summary")
  data.frame(fit = fit$fit$par,
             true = unlist(params),
             init = unlist(params_init)
  )
  print("model fit took:")
  print(proc.time()-t)
  
  info <- fit$fit$hessian
  std_err <- sqrt(diag(solve(-info)))
  std_err
  print("results summary")
  data.frame(fit = fit$fit$par,
             sd = std_err,
             true = unlist(params),
             init = unlist(params_init)
  )
  
  
  # Fit ERGM to latest network
  ergm_1 <- ergm(results$net ~ edges + gwesp(0.5,fixed = T) + gwdegree(0.5,fixed =T))
  print("ergm summary")
  summary(ergm_1)
  
  # temporal hawkes Fit:
  times <- get_times(results$net)$times
  plot(c(0,TIME), c(-1, 1), type = "n", yaxt = "n",
       xlab = "Value", ylab = "", main = "Vector on a Number Line")
  abline(h = 0, col = "gray", lwd = 2)
  points(times, rep(0, length(times)), pch = 19, col = "blue", cex = 1.5)
  
  # temporal hawkes fit suggests its not hawkesian ! yes !
  fit_temp <- fit_temporal_hawkes(params_init = list(mu = 0.1,
                                                     beta = 1,
                                                     K = 0.1),
                                  realiz = data.frame(t = times,
                                                      n=length(times)),
                                  windowT = c(0,TIME),
                                  trace = 0,
                                  maxit = 1000
  )
  fit$par
  data.frame(fitted = fit$par,
             se = diag(solve(-fit$hessian))
  )
  
  # goodness of fit:
  # suggests that data could have come from this hawkes process
  KS_test_temp <- ks_test_pval_temporal(realiz = data.frame(t = times,
                                                            n = rep(length(times),length(times))),
                                        windowT = c(0,TIME),
                                        hawkes_par = fit$par
  )
  
  KS_test_net <- ks_test_pval_hawkesGrowthNet(params = params,
                                              mark_filtration = results$net,
                                              time_window = c(0,TIME))
  
  
}

# ==============================================================================
# STUDY 1: Consistency Analysis (Sliding Window / Increasing T)
# ==============================================================================
# Goal: Show that as TIME increases, the variance of estimates decreases and 
# means converge to truth.
# ==============================================================================

if(RUN_CONSISTENCY){
  
  # 1. Define Time Windows to test
  # We will simulate independent realizations of length T = 10, 30, 50, 100
  time_windows <- c(5, 10, 20, 50,100) 
  N_SIMS_CONSISTENCY <- 20 # Keep small for demonstration, increase for paper
  
  # Parameters (Standard/Stable regime)
  params_true <- list(mu = 10,
                      beta_overall = 1,
                      K = 0.5,
                      beta_edges = 1
                      )
  
  # Setup Cluster
  cl <- make_cluster(N_CORES)
  clusterExport(cl, c("params_true", "PMF_mark_BA", "cond_intensity", 
                      "sim_hawkesGrowthNet", "fit_hawkesGrowthNet"))
  
  # Storage for results
  consistency_results <- data.frame()
  
  print("Starting Consistency Study...")
  
  for(curr_time in time_windows){
    print(paste0("Simulating and Fitting for Time Window T = ", curr_time))
    
    # Export current time to cluster
    clusterExport(cl, "curr_time", envir = environment())
    
    # Parallel Simulation & Fitting Loop
    res_list <- parLapply(cl, 1:N_SIMS_CONSISTENCY, function(i){
      
      # A. Simulate
      sim_res <- tryCatch({
        sim_hawkesGrowthNet(params = params_true,
                            time_window = c(0, curr_time),
                            PMF_mark = PMF_mark_BA,
                            cond_intensity = cond_intensity,
                            hashed_edges = TRUE,
                            mu_multiplier = 3,
                            verbose = FALSE)
      }, error = function(e) return(NULL))
      
      if(is.null(sim_res)) return(NULL)
      
      # B. Fit
      # Randomized init to test robustness
      params_init <- list(mu = runif(1, 1, 10), 
                          beta_overall = runif(1, 0.5, 2),
                          K = runif(1, 0.1, 0.9), 
                          beta_edges = runif(1, 0.5, 2))
      
      fit_res <- tryCatch({
        fit_hawkesGrowthNet(params_init = params_init,
                            time_window = c(0, curr_time),
                            mark_filtration = sim_res$net,
                            PMF_mark = PMF_mark_BA,
                            maxit = 1000,
                            grad = FALSE, 
                            cache_intensity = FALSE, # Disable cache for BA safety
                            verbose = FALSE)
      }, error = function(e) return(NULL))
      
      if(is.null(fit_res)) return(NULL)
      
      # Return row
      return(data.frame(
        sim_id = i,
        time_window = curr_time,
        param = names(fit_res$fit$par),
        estimate = as.numeric(fit_res$fit$par),
        true_value = as.numeric(unlist(params_true)[names(fit_res$fit$par)])
      ))
    })
    
    # Bind results
    res_df <- do.call(rbind, res_list)
    consistency_results <- rbind(consistency_results, res_df)
  }
  
  stopCluster(cl)
  
  # ==========================
  # Visualization
  # ==========================
  # Calculate Bias and RMSE
  summary_stats <- consistency_results %>%
    group_by(time_window, param) %>%
    summarise(
      mean_est = mean(estimate),
      sd_est = sd(estimate),
      rmse = sqrt(mean((estimate - true_value)^2)),
      true_val = mean(true_value)
    )
  
  print(summary_stats)
  
  # Plot 1: Boxplots of convergence
  p_cons <- ggplot(consistency_results, aes(x = factor(time_window), y = estimate)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, fill="lightblue") +
    geom_jitter(width=0.2, alpha=0.3) +
    geom_hline(aes(yintercept = true_value), color = "red", linetype = "dashed", size=1) +
    facet_wrap(~param, scales = "free_y") +
    labs(title = "Parameter Consistency vs Time Window (T)",
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
    labs(title = "RMSE Decay as Data Increases",
         x = "Time Window Length (T)",
         y = "Root Mean Squared Error") +
    theme_bw()
  
  print(p_rmse)
}

# ==============================================================================
# STUDY 2: Explosive Regime Analysis
# ==============================================================================
# Goal: Set beta_overall and beta_edges -> 0.
# 1. beta_overall -> 0 means the memory of the process never decays. 
#    If K > 0, the integral of intensity diverges (Explosive / Super-critical).
# 2. beta_edges -> 0 means the preferential attachment logic considers ALL past
#    nodes equally (no time decay on degree relevance).
# ==============================================================================

if(RUN_EXPLOSIVE){
  
  # Define Explosive Parameters
  # Low beta with K close to beta (or K > beta) causes criticality/explosion
  params_explosive <- list(
    mu = 2,
    beta_overall = 0.05, # Very slow decay (Long memory)
    K = 0.1,             # Branching ratio n* = K/beta = 2 (Super-critical > 1)
    beta_edges = 0.01    # Degrees from ancient history define attachment just as much as recent
  )
  
  # Compare with Stable Parameters
  params_stable <- list(
    mu = 2,
    beta_overall = 2.0,
    K = 0.5,             # Branching ratio n* = 0.25 (Sub-critical < 1)
    beta_edges = 1.0
  )
  
  print("Simulating Explosive Regime...")
  
  # Simulate Explosive
  # Note: simulation might get very slow as N grows, use small window
  sim_exp <- sim_hawkesGrowthNet(params = params_explosive,
                                 time_window = c(0, 50), # Longer window to show curve
                                 PMF_mark = PMF_mark_BA,
                                 cond_intensity = cond_intensity,
                                 hashed_edges = TRUE,
                                 verbose = TRUE, # Watch it grow
                                 mu_multiplier = 10) # Need high bound for explosive
  
  print("Simulating Stable Regime...")
  sim_stable <- sim_hawkesGrowthNet(params = params_stable,
                                    time_window = c(0, 50),
                                    PMF_mark = PMF_mark_BA,
                                    cond_intensity = cond_intensity,
                                    hashed_edges = TRUE,
                                    verbose = FALSE,
                                    mu_multiplier = 5)
  
  # ==========================
  # Visualization: Cumulative Events
  # ==========================
  df_exp <- data.frame(t = sim_exp$events$t, 
                       N = 1:length(sim_exp$events$t), 
                       Type = "Explosive (Low Beta)")
  
  df_stable <- data.frame(t = sim_stable$events$t, 
                          N = 1:length(sim_stable$events$t), 
                          Type = "Stable (High Beta)")
  
  df_compare <- rbind(df_exp, df_stable)
  
  p_expl <- ggplot(df_compare, aes(x = t, y = N, color = Type)) +
    geom_line(size = 1.2) +
    labs(title = "Explosive vs Stable Process Dynamics",
         subtitle = "Explosive: Beta -> 0 (Infinite Memory) | Stable: Beta >> 0",
         x = "Time",
         y = "Cumulative Number of Events (N)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p_expl)
  
  # ==========================
  # Visualization: Network Structure Impact
  # ==========================
  # When beta_edges is low, ancient nodes (the first ones) accumulate massive degree
  # because their 'weight' never decays. This creates "Super Hubs" (Star-like).
  
  op <- par(mfrow=c(1,2))
  
  # Plot Stable Network
  plot(sim_stable$net, main="Stable Network\n(Recent Activity Matters)", 
       vertex.cex = 0.5, edge.col="gray")
  
  # Plot Explosive Network
  # We expect the oldest nodes (ID 1, 2, 3) to have disproportionately high degree
  plot(sim_exp$net, main="Explosive/Memory Network\n(History Never Dies)", 
       vertex.cex = 0.5, edge.col="gray")
  
  par(op)
  
  # Check max degree
  max_deg_stable <- max(degree(sim_stable$net))
  max_deg_exp <- max(degree(sim_exp$net))
  
  print(paste("Max Degree Stable:", max_deg_stable))
  print(paste("Max Degree Explosive:", max_deg_exp))
}
