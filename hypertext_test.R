dat <- results_hypertext_full

t0 <- get_times(dat$net)
t1 <- get_times(dat$gof$nets_sim[[1]])
t2 <- get_times(dat$gof_decay05$nets_sim[[1]])
t3 <- get_times(dat$gof_star_esp$nets_sim[[1]])

hist(t0$times)
hist(t1$times)
hist(t2$times)
hist(t3$times)
