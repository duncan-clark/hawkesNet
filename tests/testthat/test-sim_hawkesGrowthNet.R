# Unit tests for sim_hawkesGrowthNet

test_that("sim_hawkesGrowthNet errors when params are invalid", {
  bad_params <- list(mu = 0, beta_overall = 1, K = 0.5, beta_edges = 1, m = 1)
  expect_error(
    sim_hawkesGrowthNet(
      params = bad_params,
      time_window = c(0, 1),
      PMF_mark = PMF_mark_BA,
      cond_intensity = cond_intensity,
      verbose = FALSE
    ),
    "Invalid point process parameters"
  )
})

test_that("sim_hawkesGrowthNet returns list with events, net, accept_probs for valid params", {
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(42)
  out <- sim_hawkesGrowthNet(
    params = params,
    time_window = c(0, 2),
    PMF_mark = PMF_mark_BA,
    cond_intensity = cond_intensity,
    verbose = FALSE,
    mu_multiplier = 5,
    truncation = 50
  )
  expect_type(out, "list")
  expect_true("events" %in% names(out))
  expect_true("net" %in% names(out))
  expect_true("accept_probs" %in% names(out))
  expect_true(is.list(out$events))
  expect_true("t" %in% names(out$events))
  expect_true("n" %in% names(out$events))
  expect_s3_class(out$net, "network")
  expect_type(out$accept_probs, "double")
})

test_that("sim_hawkesGrowthNet event times lie in time_window", {
  params <- list(mu = 0.3, beta_overall = 1, K = 0.4, beta_edges = 0.5, m = 1)
  set.seed(123)
  out <- sim_hawkesGrowthNet(
    params = params,
    time_window = c(0, 3),
    PMF_mark = PMF_mark_BA,
    cond_intensity = cond_intensity,
    verbose = FALSE,
    mu_multiplier = 5,
    truncation = 50
  )
  if (length(out$events$t) > 0) {
    expect_true(all(out$events$t >= 0 & out$events$t <= 3))
  }
})

test_that("sim_hawkesGrowthNet net has time attribute on vertices when non-empty", {
  params <- list(mu = 1, beta_overall = 1, K = 0.5, beta_edges = 0.5, m = 1)
  set.seed(999)
  out <- sim_hawkesGrowthNet(
    params = params,
    time_window = c(0, 5),
    PMF_mark = PMF_mark_BA,
    cond_intensity = cond_intensity,
    verbose = FALSE,
    mu_multiplier = 5,
    truncation = 100
  )
  n_verts <- network::network.size(out$net)
  if (n_verts > 0) {
    expect_true("time" %in% network::list.vertex.attributes(out$net))
  }
})
