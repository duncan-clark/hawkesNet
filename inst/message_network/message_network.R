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
                      truncation = 200,
                      verbose = TRUE
                      )
init_lik$loglik

trunc_net <- filtration_to_net(net,7)
summary(trunc_net,print.adj = F)
plot(trunc_net)

fit <- fit_hawkesGrowthNet(params_init = params_init,
                           time_window = c(0,max(get_times(trunc_net)$times)),
                           mark_filtration = trunc_net,
                           PMF_mark = PMF_mark_CS,
                           formula_RHS = "edges + triangles + star(c(2,3))",
                           truncation = 200,
                           grad = FALSE,
                           trace = 1,
                           maxit = 100)
fit$par
cbind(fit$fit$par,
      unlist(params_init)
      )

# interpretation of parameters:



