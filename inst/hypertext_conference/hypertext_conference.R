library(hawkesNet)
library(network)
library(sna)
library(ernm)
library(dplyr)
library(parallel)

N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 7))
MAX_ITER <- 5000

# function to plot pp on line:
pp_line_plot <- function(t,title= NULL){
  plot(c(min(t),max(t)), c(-1, 1), type = "n", yaxt = "n",
       xlab = "Value", ylab = "", main = paste0("Vector on a Number Line: ",title))
  abline(h = 0, col = "gray", lwd = 2)
  points(t, rep(0, length(t)), pch = 19, col = "blue", cex = 1.5)
}

# explorator function for processing data from socioNet
process_dat <- function(dat,
                        title = NULL,
                        windowT = NULL,
                        head_drop = 0,
                        tail_drop = 0
                        ){
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
  # make between 0 and 1
  print(paste0("raw tiems summary"))
  print(summary(t))
  t <- (t - min(t))/(max(t))
  print(paste0("scaled times summary"))
  print(summary(t))
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
  
  # restrict to window:
  if(is.null(windowT)){
    windowT <- c(0,1)
  }
  edges <- edges[edges$time >= windowT[1] & edges$time <= windowT[2],]
  
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
  # convert to events and mark filtration:
  tmp_edges <- tmp_edges[(head_drop+1):nrow(tmp_edges)-tail_drop,]
  # rescale:
  tmp_edges <- tmp_edges %>%
    mutate(time = (time - min(time)) / (max(time) - min(time)))
  
  net <- as.network(tmp_edges %>% select(from,to) %>% as.matrix(), matrix.type = "edgelist",directed =F)
  set.edge.attribute(net, "time", tmp_edges$time)
  set.vertex.attribute(net, "time", tmp_edges$time)
  
  
  # plot times:
  t <- get_times(net)$times
  print(paste0("network times are length ",length(t)))
  
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
              times = get_times(net)$times,
              edges = tmp_edges,
              fit_temp = fit_temp))
}

net_dat <- read.table('hawkesNet/data/ht09_contact_list.dat')
net_dat <- data.frame(
  from = net_dat$V2,
  to = net_dat$V3,
  time = net_dat$V1/20
)
# make not simple
net_dat$time <- net_dat$time + rnorm(nrow(net_dat),0,0.001)
pp_line_plot(net_dat$time)
# first day only 
net_dat <- net_dat[net_dat$time < (3600*24)/20,]
pp_line_plot(net_dat$time)
# remove part - presumabley after conference where there are few interactions
net_dat <- net_dat[net_dat$time < 2200,]
pp_line_plot(net_dat$time)
# use the window since isolated points can make the hawkes optimization fail
# drop first 2 and last isolated point:
dat_clean <- process_dat(net_dat,
                         title = 'hypertext 2009 conference',
                         windowT = c(0,1),
                         head_drop = 3,
                         tail_drop = 1
                         )
dat_clean$fit_temp$par
# drop the first two and last isolated points:
dat_clean$times
# net times 
pp_line_plot(get_times(dat_clean$net)$times)
# edge times:
pp_line_plot(dat_clean$edges$time)
# node times:
pp_line_plot(dat_clean$times)
# all times:
pp_line_plot(c(dat_clean$edges$time,dat_clean$times))

# plot_kde_intensity(get_times(dat_clean$net)$times)

net <- dat_clean$net
# check network times
hist(dat_clean$times,breaks = 100)
TRUNCATION <- net %n% 'n'

# turns out initial param matter quite a bit:
# if nodes form in [0,0.1] mu = 100 we expect 10 events -> need node_lambda = 10
# but this will never work almost all nodes are added as singletons 
# need mu ~1000, node_lambda ~ 1
params_init <- list(mu = 1000,
                   beta_overall = 0.1,
                   K = 1,
                   beta_edges = 0.1,
                   node_lambda = 1,
                   CS_params = c(-10,1,1,-1)
)

# note in this network the nodes do not keep arriving!
times <- get_times(net)
max(times$node_times)
max(times$edge_times)

if(FALSE){
  init_lik <- loglik_hawkesNet(params = params_init,
                                     time_window = c(0,1),
                                     mark_filtration = net,
                                     PMF_mark = PMF_mark_CS,
                                     formula_RHS = "edges + triangles + star(c(2,3))",
                                     truncation = TRUNCATION,
                                     verbose = T,
                                     mark_decay ='activity',
                                     max_node_time = max(times$node_times),
                                     cores = 7
                                     )
  init_lik$loglik
}

fit <- fit_hawkesNet(params_init = params_init,
                           time_window = c(0,max(get_times(net)$times)),
                           mark_filtration = net,
                           PMF_mark = PMF_mark_CS,
                           formula_RHS = "edges + triangles + star(c(2,3))",
                           truncation = TRUNCATION,
                           mark_decay ='activity',
                           max_node_time = max(times$node_times),
                           trace = 1,
                           reltol = 1e-8,
                           verbose = FALSE,
                           maxit = MAX_ITER,
                           cores = N_CORES,
                           fixed_params = c("K")
                           )

results <- data.frame(
  param = names(fit$fit$par),
  estimate = fit$fit$par,
  se = sqrt(diag(solve(-fit$fit$hessian)))
)
results

# 1. Get the negative Hessian
H <- -fit$fit$hessian
eigen_H <- eigen(H)
print(eigen_H$values)

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

ergm_fit <- NULL

# in case doesn't work befreo 
saveRDS(list(fit=fit,
             temp_fit = temp_fit,
             ergm_fit = ergm_fit,
             net = net
),
file = "hypertext_conference_results.rds"
)

cl <- makeForkCluster(N_CORES)
net_list <- parLapply(1:100,
                      cl=cl,
                      function(x){
                        sim_hawkesNet(params = relist(fit$fit$par, skeleton = params_init),
                                            time_window = c(0,1),
                                            PMF_mark = PMF_mark_CS,
                                            cond_intensity = cond_intensity,
                                            hashed_edges = T,
                                            verbose = F,
                                            mu_multiplier = 2,
                                            joint_accept = F,
                                            truncation = results$net %n% 'n',
                                            formula_RHS = "edges + triangles + star(c(2,3))",
                                            mark_decay = 'activity'
                                            )
                      }
)
stopCluster(cl)

# Save all the results from the fitting:
saveRDS(list(fit=fit,
             temp_fit = temp_fit,
             ergm_fit = ergm_fit,
             net = net
             ),
        file = "hypertext_conference_results.rds"
        )

# ====================
# Goodness of Fit
# ====================

# 1. KS TEST 

# 2. Network Fit Test

# 3. Structural Waiting Times
  - edges
  - triangles
  - stars

# 4. Network Distributions Over Time:

  - Surface of good fit

# Compare against purely separable specification i.e. Poisson Process


