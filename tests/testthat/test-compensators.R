test_that("merge_fit_params merges fitted + fixed correctly", {
  params_init <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  par_vec <- c(mu = 0.6, beta_overall = 1.1, K = 0.4, beta_edges = 0.6, m = 1.2)

  pfull <- merge_fit_params(par_vec, params_init)
  expect_equal(pfull$mu, 0.6)
  expect_equal(pfull$m, 1.2)
})

test_that("merge_fit_params works with fixed_params", {
  params_init <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  par_vec <- c(mu = 0.6, beta_overall = 1.1, beta_edges = 0.6, m = 1.2)

  pfull <- merge_fit_params(par_vec, params_init, fixed_params = "K")
  expect_equal(pfull$K, 0.3)
  expect_equal(pfull$mu, 0.6)
})

test_that("compensators_hawkesNet returns increasing values", {
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(42)
  sim <- sim_hawkesNet(params, c(0, 2), PMF_mark_BA, cond_intensity,
                       verbose = FALSE, mu_multiplier = 5, truncation = 10)
  comp <- compensators_hawkesNet(params, c(0, 2), sim$net)
  expect_type(comp, "double")
  expect_true(all(diff(comp) >= 0))
})

test_that("ks_test_pval_hawkesNet returns valid p-value", {
  params <- list(mu = 2, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(42)
  sim <- sim_hawkesNet(params, c(0, 5), PMF_mark_BA, cond_intensity,
                       verbose = FALSE, mu_multiplier = 5, truncation = 10)
  n_events <- length(get_times(sim$net)$times)
  skip_if(n_events < 3, "Simulation produced too few events for KS test")
  pval <- ks_test_pval_hawkesNet(params, c(0, 5), sim$net)
  expect_type(pval, "double")
  expect_gte(pval, 0)
  expect_lte(pval, 1)
})

test_that("cond_intensity returns a list with result and func", {
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(42)
  sim <- sim_hawkesNet(params, c(0, 2), PMF_mark_BA, cond_intensity,
                       verbose = FALSE, mu_multiplier = 5, truncation = 10)
  intens <- cond_intensity(sim$net, 2.5, sim$net, PMF_mark_BA, params, truncation = 10)
  expect_type(intens, "list")
  expect_true("result" %in% names(intens))
  expect_true("func" %in% names(intens))
  expect_gt(intens$result, 0)
})

test_that("cond_intensity_inhom works", {
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(42)
  sim <- sim_hawkesNet(params, c(0, 2), PMF_mark_BA, cond_intensity,
                       verbose = FALSE, mu_multiplier = 5, truncation = 10)
  intens <- cond_intensity_inhom(sim$net, 2.5, sim$net, PMF_mark_BA, params,
                                  mu_at_t = 0.8, truncation = 10)
  expect_type(intens, "list")
  expect_true("result" %in% names(intens))
  expect_gt(intens$result, 0)
})
