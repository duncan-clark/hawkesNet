test_that("KDE background estimation works", {
  set.seed(1)
  t <- sort(runif(100, 0, 10))
  
  mu_fit <- estimate_mu_kde(t, windowT = c(0, 10))
  expect_type(mu_fit, "list")
  expect_s3_class(mu_fit$mu_fun, "function")
  expect_equal(mu_fit$windowT, c(0, 10))
  
  ch <- make_cumhaz_fun(mu_fit)
  expect_s3_class(ch$Lambda_fun, "function")
  expect_equal(ch$Lambda_fun(0), 0, tolerance = 1e-5)
  expect_gt(ch$Lambda_fun(10), 90) # Should be close to n=100
  
  rescaled <- time_rescale_by_baseline(t, mu_fit)
  expect_length(rescaled$tau, 100)
  expect_true(all(diff(rescaled$tau) >= 0))
})

test_that("network_rescale_times_by_kde works", {
  net <- network::network(10, directed = FALSE)
  set.seed(1)
  t <- sort(runif(10, 0, 100))
  network::set.vertex.attribute(net, "time", t)
  
  # Add some edges with times
  network::add.edge(net, 1, 2)
  network::set.edge.attribute(net, "time", 50, e = 1)
  
  res <- network_rescale_times_by_kde(net)
  expect_type(res, "list")
  expect_s3_class(res$net_rescaled, "network")
  
  new_t <- network::get.vertex.attribute(res$net_rescaled, "time")
  expect_equal(min(new_t), 0, tolerance = 1e-5)
})

test_that("prepare_inhomogeneous_background works", {
  net <- network::network(10, directed = FALSE)
  set.seed(1)
  t <- sort(runif(10, 0, 100))
  network::set.vertex.attribute(net, "time", t)
  
  inhom <- prepare_inhomogeneous_background(net)
  expect_type(inhom, "list")
  expect_length(inhom$mu_vec, 10)
  expect_type(inhom$integral_bg, "double")
})
