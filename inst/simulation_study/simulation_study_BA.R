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
# 1) Preferential Attachment Style Arrival
# ===================================================
TIME = 100
params <- list(mu = 10,
               beta_BA_edges = 1,
               beta_overall = 2,
               K = 0.5)
results_1 <- simThin_hawkesGrowthNet(params =  params,
                                     time_window = c(0,TIME),
                                     PMF_mark = PMF_mark_BA,
                                     cond_intensity = cond_intensity,
                                     hashed_edges = T,
                                     verbose = T,
                                     mu_multiplier = 1,
                                     joint_accept = F
                                     )

# check accept_probs:
length(results_1$accept_probs)
summary(results_1$accept_probs)
length(results_1$events$t)

plot(results_1$accept_probs)

times <- results_1$net %v% 'time'
plot(results_1$net,
     vertex.cex = times/100,
     main = 'Plot of network, larger nodes happened later')

# Set up an empty plot with appropriate x-limits and no y-axis ticks
plot(c(0,TIME), c(-1, 1), type = "n", yaxt = "n",
     xlab = "Value", ylab = "", main = "Vector on a Number Line")
abline(h = 0, col = "gray", lwd = 2)
points(results_1$events$t, rep(0, length(results_1$events$t)), pch = 19, col = "blue", cex = 1.5)
# plot the times on a number line

degs <- ernm::calculateStatistics(results_1$net ~ degree(0:10,"in"))
print(degs)
plot(degs, col = 'red')

plot(y=degs[2:10]/(results_1$net %n% 'n'),x=2:10,col = 'red')
points((2:10)**-3,x = 2:10,col = 'blue')

esps <- ernm::calculateStatistics(results_1$net ~ esp(0:10))
print(esps)
plot(esps)

# start at mu over time
tmp <- list(mu = 0.1,
            beta_BA_edges = 0.1,
            beta_overall = 0.1,
            K = 0.01)

fit <- fit_hawkesGrowthNet(params_init = tmp,
                    time_window = c(0,TIME),
                    events = results_1$events,
                    PMF_mark = PMF_mark_BA,
                    trace = 1,
                    maxit = 10)

fit$par
fit$value
length(results_1$events$t)
params

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

# power law predicted:






# ============================
# SUMMARY OF CURRENT PROBLEMS
# ============================
# - When fitting log likelihood, the optim routine sends K -> 0 and mu to the poisson process mu
# - Not sure why this is happening? It essentially wipes out all the trigerring.
# - Guess that it can't get the intensity sum and the integral to matchup - but they seem fine in code


# I think this is working correctly, but also think that the data to not have enough info in them to fit:





# NOT CURRENTLY WORKING ....
# make a nice animation:
# ani <- make_network_growth_animation(results_1$events$marks,
#                               results_1$events$t,
#                               adjust = 10,
#                               file = 'test.mp4')
# ani.replay()


