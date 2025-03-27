# read data in:
dat <- read.table('data/CollegeMsg.txt')
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
trunc_net <- filtration_to_net(net,7)
summary(trunc_net,print.adj = F)
plot(net)
plot(trunc_net)

init_lik <- loglik_hawkesGrowthNet(params = params_init,
                                   time_window = c(0,max(edges$time)),
                                   mark_filtration = trunc_net,
                                   PMF_mark = PMF_mark_CS,
                                   formula_RHS = "edges + triangles + star(c(2,3))",
                                   truncation = 100,
                                   verbose = TRUE)
init_lik$loglik

# profile:
fit <- fit_hawkesGrowthNet(params_init = params_init,
                           time_window = c(0,max(get_times(trunc_net)$times)),
                           mark_filtration = trunc_net,
                           PMF_mark = PMF_mark_CS,
                           formula_RHS = "edges + triangles + star(c(2,3))",
                           truncation = 200,
                           grad = FALSE,
                           trace = 1,
                           verbose = TRUE,
                           get_hessian = T,
                           maxit = 1000)

info <- fit$fit$hessian
std_err <- sqrt(diag(solve(-info)))
std_err
print("results summary")
data.frame(fit = fit$fit$par,
           sd = std_err,
           init = unlist(params_init)
           )

# investigate sinular hessian:
m <- info
qr_m <- qr(m, LAPACK = TRUE)  
qr_m$rank
qr_m$pivo

# looks like the arrival times are NOT hawkesian..... 
times <- get_times(trunc_net)$times
times <- get_times(net)$times
temp_fit <- fit_temporal_hawkes(params_init = list(mu = 0.1,
                                                  beta = 1,
                                                  K = 0.1),
                               realiz = data.frame(t = times,
                                                   n = length(times)),
                               windowT = c(0,max(times)),
                               trace = 1,
                               maxit = 100
                               )
temp_fit$par



