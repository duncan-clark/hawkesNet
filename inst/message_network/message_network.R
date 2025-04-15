library(hawkesGrowthNet)
library(network)
library(sna)
library(ernm)
library(dplyr)
library(parallel)
# library(ergm)

N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 7))


# read data in:
# data collected from : https://snap.stanford.edu/data/CollegeMsg.html
dat <- read.table('hawkesGrowthNet/data/CollegeMsg.txt')
names(dat) <- c('from', 'to', 'time')
# convert to days:
dat$time <- (dat$time - min(dat$tim))/(24*60*60)
# make one directional
swap <- dat$from > dat$to
dat[swap, c("from", "to")] <- dat[swap, c("to", "from")]
dat <- dat %>% distinct()
dim(dat)

# make undirected and take minimum time:
edges <- dat %>%
  group_by(from, to) %>%
  mutate(time = min(time)) %>%
  distinct()
dim(edges)


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

# convert to events and mark filtration:
net <- as.network(edges, matrix.type = "edgelist",directed =F)
set.edge.attribute(net, "time", edges$time)
set.vertex.attribute(net, "time", times$time)

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
                                     time_window = c(0,max(edges$time)),
                                     mark_filtration = filtration_to_net(net,10) ,
                                     PMF_mark = PMF_mark_CS,
                                     formula_RHS = "edges + triangles + star(c(2,3))",
                                     truncation = 10,
                                     verbose = TRUE
  )
  init_lik$loglik
  
  # investigate truncating the network:
  trunc_net <- filtration_to_net(net,21)
  summary(trunc_net,print.adj = F)
  plot(net)
  plot(trunc_net)
  
  tmp <- trunc_net
  delete.vertices(tmp,isolates(tmp))
  summary(tmp,print.adj = F)
  plot(tmp)
  
  init_lik <- loglik_hawkesGrowthNet(params = params_init,
                                     time_window = c(0,max(edges$time)),
                                     mark_filtration = tmp,
                                     PMF_mark = PMF_mark_CS,
                                     formula_RHS = "edges + triangles + star(c(2,3))",
                                     truncation = 100,
                                     verbose = TRUE,
                                     cores = N_CORES
  )
  init_lik$loglik
}


# investigate truncating the network:
net_7 <- filtration_to_net(net,7)
net_14 <- filtration_to_net(net,14)
net_21 <- filtration_to_net(net,21)

# remove isolates:
delete.vertices(net_7,isolates(net_7))
delete.vertices(net_14,isolates(net_14))
delete.vertices(net_21,isolates(net_21))

# net_list <- list(net_7,net_14,net_21)
net_list <- list(net_7,net_14)

fits <- lapply(net_list,function(net){
  fit <- fit_hawkesGrowthNet(params_init = params_init,
                             time_window = c(0,max(get_times(net)$times)),
                             mark_filtration = net,
                             PMF_mark = PMF_mark_CS,
                             formula_RHS = "edges + triangles + star(c(2,3))",
                             truncation = 100,
                             grad = FALSE,
                             trace = 1,
                             verbose = TRUE,
                             get_hessian = T,
                             maxit = 1000,
                             cores = N_CORES
  )
  fit
})

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

temp_fits <- lapply(net_list,function(net){
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
  return(temp_fit)
}
)

# fit an ergm to the networks
# ergm_fits <- lapply(list(net_7,net_14,net_21),function(net){
#   ergm_fit <- tryCatch({ergm(net ~ edges + gwesp(0.5,fixed = T) + gwdegree(0.5,fixed =T))},
#                        error = function(e) {
#                          message("Error in ergm fit: ", e)
#                          return(NA)
#                        })
#   return(ergm_fit)
# })
ergm_fits <- NULL

# Save all the results from the fitting:
saveRDS(list(fits=fits,
             temp_fits = temp_fits,
             ergm_fits = ergm_fits),
        file = "message_network_results.rds"
        )




