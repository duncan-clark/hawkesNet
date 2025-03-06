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
library(ndtv)
library(hash)
library(hawkesGrowthNet)

# ===================================================
# 2) LOLOG Transitivity Style Arrival
# ===================================================
TIME = 50
params <- list(mu = 10,
               beta_overall =2,
               K = 0.5,
               beta_edges = 1,
               CS_params = c(-1,-0.2,0.1)
)
results_2 <- simThin_hawkesGrowthNet(params =  params,
                                     time_window = c(0,TIME),
                                     PMF_mark = PMF_mark_CS,
                                     cond_intensity = cond_intensity,
                                     hashed_edges = T,
                                     verbose = T,
                                     mu_multiplier = 10,
                                     joint_accept = F,
                                     formula_RHS = "edges  + star(c(2,3))"
)

# check accept_probs:
length(results_2$accept_probs)
summary(results_2$accept_probs)
length(results_2$events$t)

plot(results_2$accept_probs)

times <- results_2$net %v% 'time'
plot(results_2$net,
     vertex.cex = times/2,
     main = 'Plot of network, larger nodes happened later')

# Set up an empty plot with appropriate x-limits and no y-axis ticks
plot(c(0,TIME), c(-1, 1), type = "n", yaxt = "n",
     xlab = "Value", ylab = "", main = "Vector on a Number Line")
abline(h = 0, col = "gray", lwd = 2)
points(results_2$events$t, rep(0, length(results_2$events$t)), pch = 19, col = "blue", cex = 1.5)
# plot the times on a number line

degs <- ernm::calculateStatistics(results_2$net ~ degree(0:10,"in"))
print(degs)
plot(degs, col = 'red')

plot(y=degs[2:10]/(results_2$net %n% 'n'),x=2:10,col = 'red')
points((2:10)**-3,x = 2:10,col = 'blue')

esps <- ernm::calculateStatistics(results_2$net ~ esp(0:10))
print(esps)
plot(esps)

# start at mu over time
tmp <- list(mu = 10,
            beta_edges = 0.1,
            beta_overall = 0.1,
            K = 0.1,
            CS_params = c(0,0,0)
)

params <- list(mu = 10,
               beta_overall = 2,
               K = 0.5,
               beta_edges = 1,
               CS_params = c(-1,-0.2,0.1)
)

# check how the logliklihoood reacts to K:
l_k <- sapply(seq(0,1,length.out =10),function(x){
  print(x)
  tmp <- params
  tmp$K <- x
  loglik_hawkesGrowthNet(params = tmp,
                         time_window = c(0,TIME),
                         events = results_2$events,
                         PMF_mark = PMF_mark_CS,
                         use_hashing = TRUE,
                         formula_RHS = "edges  + star(c(2,3))"
  )$loglik

}
)
plot(y = l_k,x = seq(0,1,length.out = 10))

# check how the logliklihoood reacts to beta overall:
l_b <- sapply(seq(0,10,length.out =10),function(x){
  print(x)
  tmp <- params
  tmp$beta_overall <- x
  loglik_hawkesGrowthNet(params = tmp,
                         time_window = c(0,TIME),
                         events = results_2$events,
                         PMF_mark = PMF_mark_CS,
                         use_hashing = TRUE,
                         formula_RHS = "edges  + star(c(2,3))"
  )$loglik
}
)
plot(y = l_b,x = seq(0,10,length.out = 10))

# check how the logliklihoood reacts to beta edges:
l_e <- sapply(seq(-2,2,length.out =10),function(x){
  print(x)
  tmp <- params
  tmp$CS_params[2] <- x
  loglik_hawkesGrowthNet(params = tmp,
                         time_window = c(0,TIME),
                         events = results_2$events,
                         PMF_mark = PMF_mark_CS,
                         use_hashing = TRUE,
                         formula_RHS = "edges  + star(c(2,3))"
  )$loglik

}
)
plot(y = l_e,x = seq(-2,2,length.out = 10))

# I think there are large regions where the likelihood is very flat:

loglik_hawkesGrowthNet(params = params,
                       time_window = c(0,TIME),
                       events = results_2$events,
                       PMF_mark = PMF_mark_CS,
                       use_hashing = TRUE,
                       formula_RHS = "edges  + star(c(2,3))"
)$loglik

fit <- fit_hawkesGrowthNet(params_init = tmp,
                           time_window = c(0,TIME),
                           events = results_2$events,
                           PMF_mark = PMF_mark_CS,
                           formula_RHS = "edges + star(c(2,3))",
                           grad = TRUE,
                           trace = 1,
                           maxit = 1000
)

fit$par
fit$value
length(results_2$events$t)
params

# look at the K beta plane:
# Define grid of K and beta values
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

    loglik_matrix[i, j] <- loglik_hawkesGrowthNet(
      params = tmp,
      time_window = c(0, TIME),
      events = results_2$events,
      PMF_mark = PMF_mark_CS,
      use_hashing = TRUE,
      formula_RHS = "edges + star(c(2))"
    )$loglik
  }
}

# Plot using persp (Base R 3D Plot)
persp(K_values, beta_values, loglik_matrix,
      theta = 30, phi = 30,
      col = "lightblue", shade = 0.5,
      xlab = "K", ylab = "Beta", zlab = "Log-Likelihood",
      main = "Log-Likelihood Surface")

# fit temporal hawkes to the times:
fit_t <- fit_temporal_hawkes(params_init = list(mu = 0.1,
                                                beta = 1,
                                                K = 0.1),
                             realiz = data.frame(t = results_1$events$t,
                                                 n = rep(results_1$events$n,length(results_1$events$t))),
                             windowT = c(0,TIME),
                             trace = 1,
                             maxit = 1000
)
fit_t$par
fit_t$value
