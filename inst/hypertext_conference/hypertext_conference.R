library(hawkesGrowthNet)
library(network)
library(sna)
library(ernm)
library(dplyr)
library(parallel)
# library(ergm)

N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 7))
MAX_ITER <- 2000

# explorator function for processing data from socioNet
process_dat <- function(dat,title = NULL,windowT = NULL){
  # make one directional
  swap <- dat$from > dat$to
  dat[swap, c("from", "to")] <- dat[swap, c("to", "from")]
  dat <- dat %>% distinct()
  # do hawkes fit on pure time occurences samples due to memory
  if(length(dat$time) >1000){
    t <- sample(dat$time,size = 1000,replace = F)
  }else{
    t <- dat$time
  }
  t <- (t - min(t))/(max(t))
  print(paste0("raw times are length ",length(dat$time)))
  fit_pure_temp <- fit_temporal_hawkes(params_init = list(mu = 0.1,
                                                          beta = 0,
                                                          K = 0.5),
                                       realiz = data.frame(t = t,
                                                           n=length(t)),
                                       windowT = c(min(t),max(t)),
                                       trace = 0,
                                       maxit = 1000)
  print("params for pure time hawkes on 1000 points are:")
  print(fit_pure_temp$par)
  
  int_arrivals <- diff(t)
  print("sense check")
  print("Mean inter arrival:")
  print(mean(int_arrivals))
  print(" exp(-beta(mean interarrival))")
  print( exp(-fit_pure_temp$par[2] * mean(int_arrivals)))
  
  
  # make undirected and take minimum time:
  edges <- dat %>%
    group_by(from, to) %>%
    mutate(time = min(time)) %>%
    distinct()
  dim(edges)
  
  # scale time to be between 0 and 1 
  edges$time <- (edges$time - min(edges$time))
  edges$time <- edges$time / max(edges$time)
  
  # get node arrival times:
  from_times <- edges %>%
    group_by(from) %>%
    mutate(time = min(time)) %>%
    ungroup()
  to_times <- edges %>%
    group_by(to) %>%
    mutate(time = min(time)) %>%
    ungroup()
  # get all times:
  times <- from_times %>%
    full_join(to_times, by = c("from" = "to")) %>%
    mutate(time = pmin(time.x, time.y, na.rm = T)) %>%
    select(from, time) %>%
    distinct() %>%
    arrange(time)
  times$from_id <- 1:nrow(times)
  times
  
  tmp_edges <- data.frame(from = times$from_id[match(edges$from, times$from)],
                          to = times$from_id[match(edges$to, times$from)],
                          time = edges$time
  )
  tmp_edges
  
  # convert to events and mark filtration:
  net <- as.network(tmp_edges %>% select(from,to) %>% as.matrix(), matrix.type = "edgelist",directed =F)
  set.edge.attribute(net, "time", edges$time)
  set.vertex.attribute(net, "time", times$time)
  
  
  # plot times:
  t <- get_times(net)$times
  print(paste0("network times are length ",length(t)))
  
  
  # check if hawkesian
  if(is.null(windowT)){
    windowT <- c(0,1)
  }
  t <- t[t > windowT[1] & t < windowT[2]]
  
  plot(c(min(t),max(t)), c(-1, 1), type = "n", yaxt = "n",
       xlab = "Value", ylab = "", main = paste0("Vector on a Number Line: ",title))
  abline(h = 0, col = "gray", lwd = 2)
  points(t, rep(0, length(t)), pch = 19, col = "blue", cex = 1.5)
  
  fit_temp <- fit_temporal_hawkes(params_init = list(mu = 0.1,
                                                     beta = 0,
                                                     K = 0.5),
                                  realiz = data.frame(t = t,
                                                      n=length(t)),
                                  windowT = windowT,
                                  trace = 0,
                                  maxit = 1000
  )
  print(fit_temp$par)
  # check average inter arrival * exp(-beta(t-t))
  int_arrivals <- diff(t)
  print("sense check")
  print("Mean inter arrival:")
  print(mean(int_arrivals))
  print(" exp(-beta(mean interarrival))")
  print( exp(-fit_temp$par[2] * mean(int_arrivals)))
  
  
  return(list(net = net,
              times = times,
              edges = edges,
              fit_temp = fit_temp))
}

dat <- read.table('hawkesGrowthNet/data/ht09_contact_list.dat')
dat <- data.frame(
  from = dat$V2,
  to = dat$V3,
  time = dat$V1/20
)
# make not simple
dat$time <- dat$time + rnorm(nrow(dat),0,0.001)
# first day only 
dat <- dat[dat$time < (3600*24)/20,]
hist(dat$time,breaks = 100)
# remove part - presumabley after conference where there are few interactions
dat <- dat[dat$time < 2200,]
hist(dat$time,breaks = 100)
dat_clean <- process_dat(dat,title = 'hypertext 2009 conference',windowT = c(0,1))
net <- dat_clean$net
# chhck the network times :
hist(dat_clean$times$time,breaks = 100)
# looks good!

TRUNCATION <- net %n% 'n'

# fit the model to the network:
params_init <- list(mu = length(get_times(net)$times)/max(get_times(net)$times),
                   beta_overall = 0.1,
                   K = 0.1,
                   beta_edges = 0.1,
                   node_lambda = 1,
                   CS_params = c(-10,0,0,0)
)

if(FALSE){
  init_lik <- loglik_hawkesGrowthNet(params = params_init,
                                     time_window = c(0,1),
                                     mark_filtration = net,
                                     PMF_mark = PMF_mark_CS,
                                     formula_RHS = "edges + triangles + star(c(2,3))",
                                     truncation = TRUNCATION,
                                     verbose = TRUE,
                                     cores = N_CORES
                                     )
  init_lik$loglik
}

fit <- fit_hawkesGrowthNet(params_init = params_init,
                           time_window = c(0,max(get_times(net)$times)),
                           mark_filtration = net,
                           PMF_mark = PMF_mark_CS,
                           formula_RHS = "edges + triangles + star(c(2,3))",
                           truncation = TRUNCATION,
                           grad = FALSE,
                           trace = 1,
                           verbose = TRUE,
                           get_hessian = T,
                           maxit = MAX_ITER,
                           cores = N_CORES
)
fit

# info <- fit$fit$hessian
# std_err <- sqrt(diag(solve(-info)))
# std_err
# print("results summary")
# data.frame(fit = fit$fit$par,
#            sd = std_err,
#            init = unlist(params_init)
#            )
# 
# # investigate singular hessian:
# m <- info
# qr_m <- qr(m, LAPACK = TRUE)  
# qr_m$rank
# qr_m$pivo

times <- get_times(net)$times
temp_fit <- fit_temporal_hawkes(params_init = list(mu = 0.1,
                                                   beta = 1,
                                                   K = 0.1),
                                realiz = data.frame(t = times,
                                                    n = length(times)),
                                windowT = c(0,max(times)),
                                trace = 1,
                                maxit = 1000
)

# fit an ergm to the networks
times <- get_times(net)$node_times
diff_mat <- outer(times, times, FUN = function(a, b) abs(a - b))
ergm_fit <- ergm(net ~ edges + gwesp(0.5,fixed = T) + gwdegree(0.5,fixed =T) + edgecov(diff_mat))
print("ergm summary")
summary(ergm_fit)

# Save all the results from the fitting:
saveRDS(list(fit=fit,
             temp_fit = temp_fit,
             ergm_fit = ergm_fit,
             net = net
             ),
        file = "hypertext_conference_results.rds"
        )




