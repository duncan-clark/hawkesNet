devtools::install_github("https://github.com/duncan-clark/hawkesNet",
                         auth_token = "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")

load("~/MEGA/papers/crime/knife/all_data.RData")

library(hawkesGrowthNet)
library(dplyr)
## pkg example
## wrangle as Duncan
data(CollegeMsg, package = "hawkesGrowthNet")
dat <- CollegeMsg
names(dat) <- c('from', 'to', 'time')
# convert to days:
dat$time <- (dat$time - min(dat$tim))/(24*60*60)
swap <- dat$from > dat$to
dat[swap, c("from", "to")] <- dat[swap, c("to", "from")]
dat <- dat |> dplyr::distinct()
## test set
#dat <- dat [1:1000,]
######################
edges <- dat %>%
  group_by(from, to) %>%
  mutate(time = min(time)) %>%
    distinct()
## get node arrival times:
from_times <- edges %>%
  group_by(from) %>%
  mutate(time = min(time)) %>%
  ungroup()
to_times <- edges %>%
  group_by(to) %>%
  mutate(time = min(time)) %>%
  ungroup()
## get all times:
times <- from_times %>%
  full_join(to_times, by = c("from" = "to"), relationship = "many-to-many") %>%
  mutate(time = pmin(time.x, time.y, na.rm = T)) %>%
  select(from, time) %>%
  distinct() %>%
  arrange(time)

## convert to events and mark filtration:
net <- network::as.network(edges, matrix.type = "edgelist",directed = FALSE)
network::set.edge.attribute(net, "time", edges$time)
network::set.vertex.attribute(net, "time", times$time)

require(network)
params_init <- list(mu = length(get_times(net)$times)/max(get_times(net)$times),
                   beta_overall = 0.1,
                   K = 0.1,
                   beta_edges = 0.1,
                   node_lambda = 1,
                   CS_params = c(-10,0,0,0))
require(ernm)
require(sna)
loglik_hawkesGrowthNet(params = params_init,
                                     time_window = c(0,max(edges$time)),
                                     mark_filtration = filtration_to_net(net,10) ,
                                     PMF_mark = PMF_mark_CS,
                                     formula_RHS = "edges + triangles + star(c(2,3))",
                                     truncation = 10,
                                     verbose = TRUE
