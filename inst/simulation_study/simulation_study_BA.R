# This script fits a BA model:

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
library(ndtv)
library(hash)
library(hawkesGrowthNet)

# ===================================================
# Change Statistic Mark Generation
# ===================================================
TIME <- 50
INVESTIGATE = FALSE
SIMULATE = TRUE
PAPER_OUTPUT = FALSE

N_SIMS = 21
N_CORES <- detectCores() - 1
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 16))

SEED <- 01267

params  <- list(mu = 10,
                beta_overall = 2,
                K = 0.5,
                beta_edges = 0.5
                )

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
    library(ndtv)
    library(hash)
    library(hawkesGrowthNet)
  })
  clusterExport(cl, c("params",
                      "TIME",
                      "TRUNCATION",
                      "SEED"
  ))
  return(cl)
}


if(SIMULATE){
  # make the cluster:
  t <- proc.time()
  cl <- make_cluster(N_CORES)
  set.seed(SEED)
  
  sims <- parLapply(cl=cl,1:N_SIMS,function(x){
    results <- sim_hawkesGrowthNet(params =  params,
                                   time_window = c(0,TIME),
                                   PMF_mark = PMF_mark_BA,
                                   cond_intensity = cond_intensity,
                                   hashed_edges = T,
                                   verbose = F,
                                   mu_multiplier = 3,
                                   joint_accept = F
    )
    return(results)
  })
  fits <- NULL
  params_init <- list(mu = 10,
                      beta_overall = 0.1,
                      K = 0.1,
                      beta_edges = 0.1)
  clusterExport(cl, c("params_init"))
  fits <- parLapply(cl=cl,sims,function(x){
    fit <- tryCatch({
      fit_hawkesGrowthNet(
        params_init = params_init,
        time_window = c(0, TIME),
        mark_filtration = x$net,
        PMF_mark = PMF_mark_BA,
        grad = FALSE,
        trace = 0,
        maxit = 300
      )
    }, error = function(e) {
      message("Error in fit_hawkesGrowthNet: ", e$message)
      return(NULL)  # Return NULL or an alternative default if an error occurs
    })
    return(fit)
  })
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
               temp_hawkes_fits = temp_hawkes_fits),
          file = "results_BA.RDS")
  stopCluster(cl)
  print("Simulating and fitting took:")
  print((t - proc.time())[3])
}

if(PAPER_OUTPUT){
  
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
  
  # ==========================
  # RESULTS TABLE
  # ==========================
  # get mean and sd of parameters from fits:
  keep <- which(sapply(fits,length)!=0)
  paste0("keeping ",length(keep), " of ", N_SIMS," fits")
  sapply(fits,function(x){x$convergence==1})
  # check if any converged:
  
  estims <- do.call(rbind,lapply(fits[keep],function(x){
    return(as.data.frame(t(x$par),names = names(x$par)))
  }))
  
  results <- data.frame(mean = colMeans(estims),
                        sd = apply(estims,2,sd),
                        true = unlist(params),
                        init = unlist(params_init)
  )
  print(results)
  
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
  marked_p_vals <- mapply(sims[keep],lapply(1:length(fits[keep]),function(x){fits[keep][[x]]$par}),FUN = function(x,y){
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
  results <- sim_hawkesGrowthNet(params =  params,
                                 time_window = c(0,TIME),
                                 PMF_mark = PMF_mark_BA,
                                 cond_intensity = cond_intensity,
                                 hashed_edges = T,
                                 verbose = T,
                                 mu_multiplier = 10,
                                 joint_accept = F)
  
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
       vertex.cex = times/5,
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
  
  # start at mu over time
  tmp <- list(mu = 5,
              beta_overall = 0.1,
              K = 0.1,
              beta_edges = 0.1)
  
  # ============================================================
  # Plot single varialbe logliklihoood functions - for debugging
  # ============================================================
  
  # K
  l_k <- sapply(seq(0,1,length.out=10),function(x){
    print(x)
    tmp <- params
    tmp$K <- x
    loglik_hawkesGrowthNet(params = tmp,
                           time_window = c(0,TIME),
                           mark_filtration = results$net,
                           PMF_mark = PMF_mark_BA,
                           verbose = TRUE
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
                           PMF_mark = PMF_mark_BA,
                           verbose = TRUE
    )$loglik
  })
  plot(y = l_b,x = seq(0,10,length.out = 20))
  
  # beta edges
  l_e <- sapply(seq(0,2,length.out =20),function(x){
    print(x)
    tmp <- params
    tmp$beta_edges <- x
    loglik_hawkesGrowthNet(params = tmp,
                           time_window = c(0,TIME),
                           mark_filtration = results$net,
                           PMF_mark = PMF_mark_BA,
                           verbose = TRUE
    )$loglik
  })
  plot(y = l_e,x = seq(0,2,length.out = 20))
  
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
      
      loglik_matrix[i, j] <- loglik_hawkesGrowthNet(params = tmp,
                                                    time_window = c(0,TIME),
                                                    mark_filtration = results$net,
                                                    PMF_mark = PMF_mark_BA,
                                                    verbose = TRUE)$loglik
    }
  }
  # Plot using persp (Base R 3D Plot)
  persp(K_values, beta_values, loglik_matrix,
        theta = 30, phi = 30,
        col = "lightblue", shade = 0.5,
        xlab = "K", ylab = "Beta", zlab = "Log-Likelihood",
        main = "Log-Likelihood Surface")
  # I think there are large regions where the likelihood is very flat:
  
  # ======================================================
  # Fit the model
  # =====================================================
  
  l <- loglik_hawkesGrowthNet(params = params,
                              time_window = c(0,TIME),
                              mark_filtration = results$net,
                              PMF_mark = PMF_mark_BA,
                              edge_hash_list = NULL,
                              verbose = TRUE
  )
  l$loglik
  print(proc.time()-t)
  Rprof(NULL)
  
  fit <- fit_hawkesGrowthNet(params_init = params_init,
                             time_window = c(0,TIME),
                             mark_filtration = results$net,
                             PMF_mark = PMF_mark_BA,
                             grad = F,
                             trace = 1,
                             maxit = 500)
  data.frame(fit = fit$fit$par,
             true = unlist(params),
             init = unlist(params_init)
  )
  
  
  fit$value
  length(results$events$t)
  
  # Fit ERGM to latest network
  ergm_1 <- ergm(results$net ~ edges + gwesp(0.5,fixed = T) + gwdegree(0.5,fixed =T))
  summary(ergm_1)
}
  
  