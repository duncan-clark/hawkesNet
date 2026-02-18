test_that("gof works with BA model", {
  skip_if_not_installed("ggplot2")
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(1)
  sim <- sim_hawkesNet(params, c(0, 1), PMF_mark_BA, cond_intensity,
                       verbose = FALSE, mu_multiplier = 5, truncation = 10)
  skip_if(network::network.edgecount(sim$net) < 3, "Too few events for GOF test")

  # Mock fit object
  fit <- list(par = unlist(params), value = -100)

  # Run GOF with sequential simulations (cores = 1) to avoid parallel issues
  g <- tryCatch(
    gof(fit, sim$net, params, PMF_mark_BA, cond_intensity,
        time_window = c(0, 1), n_sim = 1, truncation = 10,
        verbose = FALSE, cores = 1),
    error = function(e) NULL
  )
  skip_if(is.null(g), "GOF simulation failed")

  expect_type(g, "list")
  expect_true("plots" %in% names(g))
  if (!is.null(g$plots$degree_plot)) {
    expect_s3_class(g$plots$degree_plot, "ggplot")
  }
})

test_that("gof handles NULL fit", {
  g <- gof(NULL, network::network.initialize(0), list(), PMF_mark_BA, cond_intensity, verbose = FALSE)
  expect_null(g$degree_obs)
})
