# Unit tests for fit_hawkesNet and loglik_hawkesNet

# Helper: minimal filtration (small network with times) for BA
make_minimal_filtration_BA <- function() {
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(1)
  sim <- sim_hawkesNet(
    params = params,
    time_window = c(0, 3),
    PMF_mark = PMF_mark_BA,
    cond_intensity = cond_intensity,
    verbose = FALSE,
    mu_multiplier = 5,
    truncation = 30
  )
  sim$net
}

test_that("loglik_hawkesNet returns list with loglik and intens_funcs", {
  net <- make_minimal_filtration_BA()
  times <- get_times(net)$times
  skip_if(length(times) < 2, "Need at least 2 event times for loglik")
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  out <- loglik_hawkesNet(
    params = params,
    time_window = range(times),
    mark_filtration = net,
    PMF_mark = PMF_mark_BA,
    verbose = FALSE,
    truncation = 50
  )
  expect_type(out, "list")
  expect_true("loglik" %in% names(out))
  expect_true("intens_funcs" %in% names(out))
  expect_type(out$loglik, "double")
  expect_length(out$loglik, 1)
  expect_true(is.finite(out$loglik))
})

test_that("fit_hawkesNet returns expected structure", {
  net <- make_minimal_filtration_BA()
  times <- get_times(net)$times
  skip_if(length(times) < 2, "Need at least 2 event times for fit")
  params_init <- list(mu = 0.3, beta_overall = 0.8, beta_edges = 0.4, K = 0.4, m = 0.8)
  suppressMessages({
    fit <- fit_hawkesNet(
      params_init = params_init,
      time_window = range(times),
      mark_filtration = net,
      PMF_mark = PMF_mark_BA,
      maxit = 50,
      trace = 0,
      verbose = FALSE,
      cache_intensity = FALSE,
      truncation = 50
    )
  })
  expect_type(fit, "list")
  expect_true("fit" %in% names(fit))
  expect_true("fit_table" %in% names(fit))
  expect_true("hessian" %in% names(fit))
  expect_true("par" %in% names(fit$fit))
  expect_true("value" %in% names(fit$fit))
  expect_true(is.finite(fit$fit$value))
})

test_that("compensators_hawkesNet returns numeric vector", {
  net <- make_minimal_filtration_BA()
  times <- get_times(net)$times
  skip_if(length(times) < 2, "Need at least 2 event times")
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  comp <- compensators_hawkesNet(
    params = params,
    time_window = range(times),
    mark_filtration = net
  )
  expect_type(comp, "double")
  expect_length(comp, length(times))
  expect_true(all(is.finite(comp)))
})

test_that("cond_intensity returns list with result and func", {
  net <- make_minimal_filtration_BA()
  times <- get_times(net)$times
  skip_if(length(times) < 2, "Need at least 2 event times")
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  out <- cond_intensity(
    new_net = net,
    t = times[length(times)],
    mark_filtration = net,
    PMF_mark = PMF_mark_BA,
    params = params,
    truncation = 50
  )
  expect_type(out, "list")
  expect_true("result" %in% names(out))
  expect_true("func" %in% names(out))
  expect_true(is.finite(out$result))
  expect_true(out$result > 0)
})
