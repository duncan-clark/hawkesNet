test_that("gof works with BA model", {
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(1)
  sim <- sim_hawkesNet(params, c(0, 1), PMF_mark_BA, cond_intensity,
                       verbose = FALSE, mu_multiplier = 5, truncation = 10)
  
  # Mock fit object
  fit <- list(par = unlist(params), value = -100)
  
  # Run GOF with 1 simulation for speed
  g <- gof(fit, sim$net, params, PMF_mark_BA, cond_intensity,
           time_window = c(0, 1), n_sim = 1, truncation = 10, verbose = FALSE)
  
  expect_type(g, "list")
  expect_named(g, c("degree_obs", "degree_sim", "esp_obs", "esp_sim",
                    "geodist_obs", "geodist_sim", "wait_obs", "wait_sim",
                    "nodemix_obs", "nodemix_sim", "nets_sim", "plots"))
  expect_s3_class(g$plots$degree_plot, "ggplot")
})

test_that("gof handles NULL fit", {
  g <- gof(NULL, network::network(0), list(), PMF_mark_BA, cond_intensity, verbose = FALSE)
  expect_null(g$degree_obs)
})
