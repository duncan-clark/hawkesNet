test_that("KDE background estimation works", {
  set.seed(1)
  t <- sort(runif(100, 0, 10))

  mu_fit <- estimate_mu_kde(t, windowT = c(0, 10))
  expect_type(mu_fit, "list")
  expect_type(mu_fit$mu_fun, "closure")
  expect_equal(mu_fit$windowT, c(0, 10))

  ch <- make_cumhaz_fun(mu_fit)
  expect_type(ch$Lambda_fun, "closure")
  expect_equal(ch$Lambda_fun(0), 0, tolerance = 1e-5)
  expect_gt(ch$Lambda_fun(10), 90)

  rescaled <- time_rescale_by_baseline(t, mu_fit)
  expect_length(rescaled$tau, 100)
  expect_true(all(diff(rescaled$tau) >= 0))
})

test_that("network_rescale_times_by_kde works", {
  net <- network::network.initialize(10, directed = FALSE)
  set.seed(1)
  t <- sort(runif(10, 0, 100))
  network::set.vertex.attribute(net, "time", t)

  for (i in 1:9) network::add.edge(net, i, i + 1)
  network::set.edge.attribute(net, "time", sort(runif(9, 0, 100)))

  res <- network_rescale_times_by_kde(net)
  expect_type(res, "list")
  expect_s3_class(res$net_rescaled, "network")

  new_t <- network::get.vertex.attribute(res$net_rescaled, "time")
  expect_equal(min(new_t), 0, tolerance = 1e-5)
})

test_that("prepare_inhomogeneous_background works", {
  net <- network::network.initialize(10, directed = FALSE)
  set.seed(1)
  t <- sort(runif(10, 0, 100))
  network::set.vertex.attribute(net, "time", t)

  inhom <- prepare_inhomogeneous_background(net)
  expect_type(inhom, "list")
  expect_length(inhom$mu_vec, 10)
  expect_type(inhom$integral_bg, "double")
})

test_that("plot_kde_background runs without error", {
  set.seed(1)
  t <- sort(runif(50, 0, 10))
  mu_fit <- estimate_mu_kde(t, windowT = c(0, 10))
  expect_no_error(plot_kde_background(mu_fit))
})
