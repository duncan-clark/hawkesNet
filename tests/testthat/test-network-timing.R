timing_triangle_history <- function() {
  net <- network::network.initialize(3, directed = FALSE)
  network::set.vertex.attribute(net, "time", c(1, 1, 1))
  network::add.edges(net, c(1, 2, 1), c(2, 3, 3))
  network::set.edge.attribute(net, "time", c(2, 2, 3))
  net
}

test_that("timing constructor makes the two alpha specifications explicit", {
  expect_equal(hawkes_timing()$model, "M0")
  expect_equal(hawkes_timing("M1", "degree")$link, "exp")
  expect_equal(hawkes_timing("M2", "triangles")$link, "logistic")
  expect_equal(hawkes_timing("M2", "degree", "logistic", 5)$scale, 5)
  expect_error(hawkes_timing("M3"), "arg")
  expect_error(hawkes_timing(scale = 0), "positive")
  expect_error(hawkes_timing(scale = Inf), "finite")
  expect_error(hawkes_timing(link = "linear"), "arg")
  expect_equal(.network_timing_bound(list(feedback_gamma = log(3)), hawkes_timing("M1")), 3)
  expect_equal(.network_timing_bound(list(feedback_gamma = -2), hawkes_timing("M2")), 1)
  expect_equal(.network_timing_bound(list(feedback_gamma = 1e4),
                                      hawkes_timing("M2", "triangles")), 2)
  expect_error(.network_timing_bound(list(feedback_gamma = 1e4), hawkes_timing("M1")),
               "non-finite alpha bound")
})

test_that("one timestamp is one complete update and M2 reweights earlier marks", {
  net <- timing_triangle_history()
  fixed <- .prepare_network_timing(net, hawkes_timing("M1", "triangles"))
  evolving <- .prepare_network_timing(net, hawkes_timing("M2", "triangles"))
  expect_equal(fixed$times, c(1, 2, 3))
  expect_equal(fixed$birth_features, c(0, 0, 1))
  expect_equal(evolving$features_after, list(0, c(0, 0), c(0, 1, 1)))
  expect_equal(.network_timing_features(evolving, 3), c(0, 0))
  params <- list(mu = .2, K = .7, beta_overall = .5, feedback_gamma = log(3))
  t <- 3.5
  alpha_closed <- 1.5
  expected_m1 <- .2 + .7 * (exp(-.5 * (t - 1)) + exp(-.5 * (t - 2)) +
                             alpha_closed * exp(-.5 * (t - 3)))
  expected_m2 <- .2 + .7 * (exp(-.5 * (t - 1)) +
                             alpha_closed * sum(exp(-.5 * (t - c(2, 3)))))
  expect_equal(hawkes_ground_intensity(t, params, net, hawkes_timing("M1", "triangles")),
               expected_m1)
  expect_equal(hawkes_ground_intensity(t, params, net, hawkes_timing("M2", "triangles")),
               expected_m2)
})

test_that("degree features use original IDs and the complete post-update graph", {
  net <- network::network.initialize(4, directed = FALSE)
  network::set.vertex.attribute(net, "time", c(1, 3, 1, 2))
  # Deliberately insert edges out of endpoint and temporal order.
  network::add.edges(net, c(2, 1, 1, 3), c(4, 3, 4, 4))
  network::set.edge.attribute(net, "time", c(3, 1, 2, 2))
  timing <- hawkes_timing("M2", "degree")
  cache <- .prepare_network_timing(net, timing)
  expect_equal(cache$birth_features, c(1 / 2, 2 / 3, 2 / 3))
  expect_equal(cache$features_after, list(1 / 2, c(2 / 3, 2 / 3),
                                         c(2 / 3, 7 / 10, 2 / 3)))
  scaled <- .prepare_network_timing(net, hawkes_timing("M1", "degree", scale = 2))
  expect_equal(scaled$birth_features, c(1 / 3, 1 / 2, 1 / 2))
  params <- list(mu = .3, K = .4, beta_overall = .7, feedback_gamma = .8)
  expected <- .3 + .4 * sum(exp(.8 * c(2 / 3, 7 / 10, 2 / 3)) *
                              exp(-.7 * (3.2 - c(1, 2, 3))))
  expect_equal(hawkes_ground_intensity(3.2, params, net, timing), expected)
  expect_equal(network::get.vertex.attribute(net, "time"), c(1, 3, 1, 2))
})

test_that("left limits exclude both the new event and its topology changes", {
  net <- timing_triangle_history()
  truncated <- network::network.copy(net)
  network::delete.edges(truncated, 3)
  params <- list(mu = .2, K = .7, beta_overall = .5, feedback_gamma = log(3))
  for (feature in c("degree", "triangles")) {
    timing <- hawkes_timing("M2", feature)
    query <- c(0, 1, 2, 2.5, 3)
    expect_equal(hawkes_ground_intensity(query, params, net, timing),
                 hawkes_ground_intensity(query, params, truncated, timing))
    expect_equal(hawkes_ground_compensator(params, c(0, 3), net, timing),
                 hawkes_ground_compensator(params, c(0, 3), truncated, timing))
  }
})

test_that("zero feedback recovers M0 for both links and both features", {
  net <- timing_triangle_history()
  params <- list(mu = .2, K = .7, beta_overall = .5)
  query <- c(.5, 1, 1.7, 2, 2.8, 3, 3.2, 5)
  expected <- vapply(query, function(t)
    .2 + .7 * sum(exp(-.5 * (t - c(1, 2, 3)[c(1, 2, 3) < t]))), numeric(1))
  base_integral <- hawkes_ground_compensator(params, c(1.8, 5), net)
  for (model in c("M0", "M1", "M2")) {
    for (feature in c("degree", "triangles")) {
      for (link in c("exp", "logistic")) {
        timing <- hawkes_timing(model, feature, link)
        expect_equal(hawkes_ground_intensity(query, params, net, timing), expected)
        expect_equal(hawkes_ground_compensator(params, c(1.8, 5), net, timing), base_integral)
      }
    }
  }
})

test_that("M2 compensator includes the reweighting between graph updates", {
  net <- timing_triangle_history()
  params <- list(mu = .2, K = .7, beta_overall = .5, feedback_gamma = log(3))
  timing <- hawkes_timing("M2", "triangles")
  expected <- .2 * 5 + .7 / .5 * (
    (1 - exp(-.5 * 4)) +
      (1 - exp(-.5)) + 1.5 * (exp(-.5) - exp(-.5 * 3)) +
      1.5 * (1 - exp(-.5 * 2)))
  expect_equal(hawkes_ground_compensator(params, c(0, 5), net, timing), expected)
  numerical <- sum(vapply(seq_len(4), function(i) {
    cuts <- c(0, 1, 2, 3, 5)
    stats::integrate(function(t) hawkes_ground_intensity(t, params, net, timing),
                     cuts[i], cuts[i + 1], rel.tol = 1e-10)$value
  }, numeric(1)))
  expect_equal(hawkes_ground_compensator(params, c(0, 5), net, timing), numerical)
  fixed <- hawkes_ground_compensator(params, c(0, 5), net, hawkes_timing("M1", "triangles"))
  expect_equal(expected - fixed, .7 * .5 / .5 * (exp(-.5) - exp(-1.5)))
})

test_that("integration handles prehistory, future events and time shifts", {
  net <- timing_triangle_history()
  shifted <- network::network.copy(net)
  network::set.vertex.attribute(shifted, "time", c(11, 11, 11))
  network::set.edge.attribute(shifted, "time", c(12, 12, 13))
  params <- list(mu = .2, K = .7, beta_overall = .5, feedback_gamma = -.6)
  for (model in c("M0", "M1", "M2")) {
    for (feature in c("degree", "triangles")) {
      timing <- hawkes_timing(model, feature)
      window <- c(2.5, 3.8)
      numerical <- sum(vapply(list(c(2.5, 3), c(3, 3.8)), function(w)
        stats::integrate(function(t) hawkes_ground_intensity(t, params, net, timing),
                         w[1], w[2], rel.tol = 1e-10)$value, numeric(1)))
      actual <- hawkes_ground_compensator(params, window, net, timing)
      expect_equal(actual, numerical)
      expect_equal(actual, hawkes_ground_compensator(params, window + 10, shifted, timing))
      expect_equal(hawkes_ground_intensity(c(2.5, 3.8), params, net, timing),
                   hawkes_ground_intensity(c(12.5, 13.8), params, shifted, timing))
      expect_equal(hawkes_ground_compensator(params, c(-1, 0), net, timing), .2)
      expect_equal(hawkes_ground_compensator(params, c(2, 2), net, timing), 0)
      expect_equal(hawkes_ground_compensator(params, window, net, timing, integral_bg = 4),
                   actual - .2 * diff(window) + 4)
    }
  }
})

test_that("empty histories and static initial graphs have no artificial events", {
  params <- list(mu = .2, K = .7, beta_overall = .5, feedback_gamma = 1)
  empty <- network::network.initialize(0, directed = FALSE)
  for (model in c("M0", "M1", "M2")) {
    timing <- hawkes_timing(model)
    expect_equal(.prepare_network_timing(empty, timing)$times, numeric(0))
    expect_equal(hawkes_ground_intensity(c(0, 1), params, empty, timing), c(.2, .2))
    expect_equal(hawkes_ground_compensator(params, c(0, 5), empty, timing), 1)
  }
  initial <- network::network.initialize(2, directed = FALSE)
  network::set.vertex.attribute(initial, "time", c(NA, NA))
  expect_equal(.prepare_network_timing(initial, hawkes_timing("M2"))$times, numeric(0))
  expect_equal(hawkes_ground_intensity(1, params, initial, hawkes_timing("M2")), .2)
  net <- network::network.initialize(4, directed = FALSE)
  network::set.vertex.attribute(net, "time", c(NA, NA, 1, 2))
  network::add.edges(net, c(1, 1, 2), c(2, 3, 3))
  network::set.edge.attribute(net, "time", c(NA, 1, 1))
  triangles <- .prepare_network_timing(net, hawkes_timing("M2", "triangles"))
  expect_equal(triangles$times, c(1, 2))
  expect_equal(triangles$birth_features, c(1, 0))
  degree <- .prepare_network_timing(net, hawkes_timing("M2", "degree"))
  expect_equal(degree$birth_features, c(2 / 3, 0))
})

test_that("invalid feedback inputs fail explicitly and caches check configuration", {
  net <- timing_triangle_history()
  params <- list(mu = .2, K = .7, beta_overall = .5)
  directed <- network::network.initialize(0, directed = TRUE)
  multiple <- network::network.initialize(0, directed = FALSE, multiple = TRUE)
  expect_error(.prepare_network_timing(directed, hawkes_timing("M1")), "undirected simple")
  expect_error(.prepare_network_timing(multiple, hawkes_timing("M2")), "undirected simple")
  data <- data.frame(t = c(1, 1, 2), i = c(1, 2, 1), j = c(2, 3, 3))
  expect_equal(.prepare_network_timing(data)$times, c(1, 2))
  expect_error(.prepare_network_timing(data, hawkes_timing("M1")), "data.frame")
  invalid <- network::network.copy(net)
  network::set.vertex.attribute(invalid, "time", c(1, Inf, 1))
  expect_error(.prepare_network_timing(invalid, hawkes_timing("M2")), "finite")
  network::set.vertex.attribute(invalid, "time", c(1, 4, 1))
  expect_error(.prepare_network_timing(invalid, hawkes_timing("M2")), "precede")
  cache <- .prepare_network_timing(net, hawkes_timing("M1"))
  expect_error(hawkes_ground_intensity(3.5, params, net, hawkes_timing("M2"), cache),
               "same timing")
  expect_equal(hawkes_ground_intensity(3.5, params, net, hawkes_timing("M1"), cache),
               hawkes_ground_intensity(3.5, params, net, hawkes_timing("M1")))
  expect_error(hawkes_ground_intensity(Inf, params, net), "finite numeric query")
  params$feedback_gamma <- NA_real_
  expect_error(hawkes_ground_intensity(1, params, net), "feedback_gamma")
  params$feedback_gamma <- 0
  params$beta_overall <- 0
  expect_error(hawkes_ground_intensity(1, params, net), "beta_overall")
})
