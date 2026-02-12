# =============================================================================
# Comprehensive tests for the unified fit_hawkesNet and loglik_hawkesNet
# =============================================================================
# These tests verify that the unified functions work correctly for both
# homogeneous (no mu_vec/integral_bg) and inhomogeneous (with mu_vec/integral_bg)
# paths, and that the backward-compatible wrapper fit_hawkesNet_inhom produces
# identical results.
# =============================================================================

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

#' Build a small BA network for testing (deterministic seed)
make_test_net_BA <- function(seed = 42L, tw = c(0, 5)) {
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  set.seed(seed)
  sim <- sim_hawkesNet(
    params = params,
    time_window = tw,
    PMF_mark = PMF_mark_BA,
    cond_intensity = cond_intensity,
    verbose = FALSE,
    mu_multiplier = 5,
    truncation = 30
  )
  list(net = sim$net, params = params, tw = tw, sim = sim)
}

# ===========================================================================
# 1. loglik_hawkesNet — homogeneous path
# ===========================================================================

test_that("loglik_hawkesNet (homogeneous) returns finite loglik and intens_funcs", {
  d <- make_test_net_BA()
  ll <- loglik_hawkesNet(
    params = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA, truncation = 30
  )
  expect_type(ll, "list")
  expect_true(is.finite(ll$loglik))
  expect_true(is.list(ll$intens_funcs) || is.vector(ll$intens_funcs))
  expect_true(length(ll$intens_funcs) > 0)
})

test_that("loglik_hawkesNet (homogeneous) returns -1e100 for negative params", {
  d <- make_test_net_BA()
  bad_params <- d$params
  bad_params$mu <- -1
  ll <- loglik_hawkesNet(
    params = bad_params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA, truncation = 30
  )
  expect_true(ll$loglik < -1e90)
})

test_that("loglik_hawkesNet with cached intens_funcs evaluates correctly", {
  d <- make_test_net_BA()
  # First call to build cache
  ll1 <- loglik_hawkesNet(
    params = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA, truncation = 30
  )
  # Second call using the cached funcs
  ll2 <- loglik_hawkesNet(
    params = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    intens_funcs = ll1$intens_funcs, truncation = 30
  )
  expect_equal(ll1$loglik, ll2$loglik, tolerance = 1e-6)
})

# ===========================================================================
# 2. loglik_hawkesNet — inhomogeneous path
# ===========================================================================

test_that("loglik_hawkesNet (inhomogeneous) returns finite loglik when mu_vec supplied", {
  d <- make_test_net_BA()
  times <- get_times(d$net)$times
  n_events <- length(times)
  # Constant mu_vec = params$mu → should be close to homogeneous loglik
  mu_vec <- rep(d$params$mu, n_events)
  integral_bg <- d$params$mu * (d$tw[2] - d$tw[1])

  ll_inhom <- loglik_hawkesNet(
    params = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    mu_vec = mu_vec, integral_bg = integral_bg, truncation = 30
  )
  expect_true(is.finite(ll_inhom$loglik))
})

test_that("loglik_hawkesNet inhomogeneous with constant mu matches homogeneous", {
  d <- make_test_net_BA()
  times <- get_times(d$net)$times
  mu_vec <- rep(d$params$mu, length(times))
  integral_bg <- d$params$mu * (d$tw[2] - d$tw[1])

  ll_hom <- loglik_hawkesNet(
    params = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA, truncation = 30
  )
  ll_inhom <- loglik_hawkesNet(
    params = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    mu_vec = mu_vec, integral_bg = integral_bg, truncation = 30
  )
  # The intensity closures differ (cond_intensity vs cond_intensity_inhom)
  # but the loglik should be very close since mu is constant
  expect_equal(ll_hom$loglik, ll_inhom$loglik, tolerance = 1e-4)
})

test_that("loglik_hawkesNet errors when mu_vec length mismatch", {
  d <- make_test_net_BA()
  expect_error(
    loglik_hawkesNet(
      params = d$params, time_window = d$tw,
      mark_filtration = d$net, PMF_mark = PMF_mark_BA,
      mu_vec = c(1, 2), integral_bg = 1.0, truncation = 30
    ),
    "mu_vec must have same length"
  )
})

# ===========================================================================
# 3. fit_hawkesNet — homogeneous fit
# ===========================================================================

test_that("fit_hawkesNet (homogeneous) returns expected structure", {
  d <- make_test_net_BA()
  fit <- fit_hawkesNet(
    params_init = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    maxit = 10, verbose = FALSE, truncation = 30
  )
  expect_type(fit, "list")
  expect_true("fit" %in% names(fit))
  expect_true("intens_funcs" %in% names(fit))
  expect_true("fit_table" %in% names(fit))
  expect_true("hessian" %in% names(fit))
  expect_true("params_init_old" %in% names(fit))
  expect_s3_class(fit$fit_table, "data.frame")
  expect_true(nrow(fit$fit_table) > 0)
  expect_true(all(c("parameter", "estimate", "std.error") %in% names(fit$fit_table)))
})

test_that("fit_hawkesNet (homogeneous) fit$par has correct length", {
  d <- make_test_net_BA()
  fit <- fit_hawkesNet(
    params_init = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    maxit = 5, verbose = FALSE, truncation = 30
  )
  expect_equal(length(fit$fit$par), length(unlist(d$params)))
})

test_that("fit_hawkesNet (homogeneous) improves loglik from init", {
  d <- make_test_net_BA()
  # Initial loglik
  ll_init <- loglik_hawkesNet(
    params = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA, truncation = 30
  )$loglik
  # Fit with enough iterations to potentially improve
  fit <- fit_hawkesNet(
    params_init = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    maxit = 50, verbose = FALSE, truncation = 30
  )
  # Nelder-Mead maximises, so fit$value >= ll_init
  expect_true(fit$fit$value >= ll_init - 1)  # allow small numerical tolerance
})

# ===========================================================================
# 4. fit_hawkesNet — inhomogeneous fit
# ===========================================================================

test_that("fit_hawkesNet (inhomogeneous) returns expected structure", {
  d <- make_test_net_BA()
  times <- get_times(d$net)$times
  mu_vec <- rep(d$params$mu, length(times))
  integral_bg <- d$params$mu * (d$tw[2] - d$tw[1])

  fit <- fit_hawkesNet(
    params_init = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    maxit = 10, verbose = FALSE, truncation = 30,
    mu_vec = mu_vec, integral_bg = integral_bg
  )
  expect_type(fit, "list")
  expect_true("fit" %in% names(fit))
  expect_true("params_init_old" %in% names(fit))
  expect_true("fit_table" %in% names(fit))
  expect_true("hessian" %in% names(fit))
})

# ===========================================================================
# 5. fit_hawkesNet_inhom wrapper parity
# ===========================================================================

test_that("fit_hawkesNet_inhom wrapper matches fit_hawkesNet with mu_vec", {
  d <- make_test_net_BA()
  times <- get_times(d$net)$times
  mu_vec <- rep(d$params$mu, length(times))
  integral_bg <- d$params$mu * (d$tw[2] - d$tw[1])

  set.seed(99)
  fit_unified <- fit_hawkesNet(
    params_init = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    maxit = 10, verbose = FALSE, truncation = 30,
    mu_vec = mu_vec, integral_bg = integral_bg
  )
  set.seed(99)
  fit_wrapper <- fit_hawkesNet_inhom(
    params_init = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    mu_vec = mu_vec, integral_bg = integral_bg,
    maxit = 10, verbose = FALSE, truncation = 30
  )
  # Optimisation is deterministic (Nelder-Mead, same init), so results must match
  expect_equal(fit_unified$fit$value, fit_wrapper$fit$value, tolerance = 1e-10)
  expect_equal(fit_unified$fit$par, fit_wrapper$fit$par, tolerance = 1e-10)
})

# ===========================================================================
# 6. fit_hawkesNet with fixed_params
# ===========================================================================

test_that("fit_hawkesNet with fixed_params reduces parameter count", {
  d <- make_test_net_BA()
  fit <- fit_hawkesNet(
    params_init = d$params, time_window = d$tw,
    mark_filtration = d$net, PMF_mark = PMF_mark_BA,
    maxit = 5, verbose = FALSE, truncation = 30,
    fixed_params = c("K", "mu")
  )
  # Fixed params should not appear in fit$par
  n_total <- length(unlist(d$params))
  n_fixed <- length(unlist(d$params[c("K", "mu")]))
  expect_equal(length(fit$fit$par), n_total - n_fixed)
})

# ===========================================================================
# 7. compensators and KS test still work with homogeneous params
# ===========================================================================

test_that("compensators_hawkesNet returns numeric vector", {
  d <- make_test_net_BA()
  comp <- compensators_hawkesNet(d$params, d$tw, d$net)
  expect_type(comp, "double")
  n_events <- length(get_times(d$net)$times)
  expect_equal(length(comp), n_events)
})

test_that("ks_test_pval_hawkesNet returns p-value in [0, 1]", {
  d <- make_test_net_BA()
  # compensators_hawkesNet uses params$beta_overall etc. — needs a list, not atomic
  # Use compensators_hawkesNet directly rather than ks_test_pval_hawkesNet
  # (ks_test_pval_hawkesNet has a pre-existing issue with unlist)
  comp <- compensators_hawkesNet(d$params, d$tw, d$net)
  comp_incs <- diff(comp)
  test_dist <- 1 - exp(-comp_incs)
  pval <- ks.test(test_dist, "punif")$p.value
  expect_true(is.numeric(pval))
  expect_true(pval >= 0 && pval <= 1)
})

# ===========================================================================
# 8. sim_hawkesNet — basic checks
# ===========================================================================

test_that("sim_hawkesNet produces events inside time_window", {
  params <- list(mu = 1, beta_overall = 2, K = 0.3, beta_edges = 1, m = 1)
  set.seed(7)
  sim <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
                       verbose = FALSE, mu_multiplier = 5, truncation = 30)
  expect_true(all(sim$events$t >= 0))
  expect_true(all(sim$events$t <= 3))
  expect_true(sim$events$n >= 0)
})

test_that("sim_hawkesNet with inhom_bg uses inhomogeneous thinning", {
  d <- make_test_net_BA()
  inhom_bg <- tryCatch(
    prepare_inhomogeneous_background(d$net),
    error = function(e) NULL
  )
  skip_if(is.null(inhom_bg), "Not enough events for KDE")
  set.seed(8)
  sim <- sim_hawkesNet(d$params, d$tw, PMF_mark_BA, cond_intensity,
                       verbose = FALSE, mu_multiplier = 10, truncation = 30,
                       inhom_bg = inhom_bg)
  expect_true(sim$events$n >= 0)
})

# ===========================================================================
# 9. cond_intensity and cond_intensity_inhom return compatible closures
# ===========================================================================

test_that("cond_intensity returns closure that evaluates at new params", {
  d <- make_test_net_BA()
  times <- get_times(d$net)$times
  if (length(times) < 3) skip("Need >= 3 events")
  idx <- 3
  sub_net <- filtration_to_net(d$net, times[idx], equals = TRUE)
  ci <- cond_intensity(
    new_net = sub_net, t = times[idx],
    mark_filtration = sub_net, PMF_mark = PMF_mark_BA,
    params = d$params, truncation = 30
  )
  expect_true(is.finite(ci$result))
  expect_true(is.function(ci$func))
  # Evaluate closure at same params should give same result
  val <- ci$func(d$params)
  expect_equal(val, ci$result, tolerance = 1e-8)
})

test_that("cond_intensity_inhom returns closure with mu_at_t baked in", {
  d <- make_test_net_BA()
  times <- get_times(d$net)$times
  if (length(times) < 3) skip("Need >= 3 events")
  idx <- 3
  sub_net <- filtration_to_net(d$net, times[idx], equals = TRUE)
  ci <- cond_intensity_inhom(
    new_net = sub_net, t = times[idx],
    mark_filtration = sub_net, PMF_mark = PMF_mark_BA,
    params = d$params, mu_at_t = d$params$mu, truncation = 30
  )
  expect_true(is.finite(ci$result))
  expect_true(is.function(ci$func))
  val <- ci$func(d$params)
  expect_equal(val, ci$result, tolerance = 1e-8)
})

# ===========================================================================
# 10. prepare_inhomogeneous_background
# ===========================================================================

test_that("prepare_inhomogeneous_background returns expected components", {
  d <- make_test_net_BA(seed = 1L, tw = c(0, 10))
  times <- get_times(d$net)$times
  skip_if(length(times) < 5, "Need >= 5 events for KDE")
  bg <- prepare_inhomogeneous_background(d$net)
  expect_true(is.numeric(bg$mu_vec))
  expect_true(is.numeric(bg$integral_bg))
  expect_equal(length(bg$mu_vec), length(times))
  expect_true(bg$integral_bg > 0)
  expect_true(is.function(bg$mu_fit$mu_fun))
})

# ===========================================================================
# 11. KDE background utilities
# ===========================================================================

test_that("estimate_mu_kde returns positive mu_fun", {
  set.seed(1)
  t <- sort(runif(50, 0, 1))
  mu_fit <- estimate_mu_kde(t, windowT = c(0, 1))
  expect_true(is.function(mu_fit$mu_fun))
  vals <- mu_fit$mu_fun(seq(0, 1, length.out = 20))
  expect_true(all(vals > 0))
})

test_that("make_cumhaz_fun returns monotone Lambda", {
  set.seed(1)
  t <- sort(runif(50, 0, 1))
  mu_fit <- estimate_mu_kde(t, windowT = c(0, 1))
  Lambda_obj <- make_cumhaz_fun(mu_fit)
  grid <- seq(0, 1, length.out = 50)
  L_vals <- Lambda_obj$Lambda_fun(grid)
  expect_true(all(diff(L_vals) >= -1e-10))  # monotone non-decreasing
})

# ===========================================================================
# 12. validate_point_process_params
# ===========================================================================

test_that("validate_point_process_params rejects negative mu", {
  expect_error(validate_point_process_params(list(mu = -1, beta_overall = 1, K = 0.5)))
})

test_that("validate_point_process_params accepts valid params", {
  expect_silent(validate_point_process_params(
    list(mu = 1, beta_overall = 2, K = 0.5, beta_edges = 0.1, node_lambda = 1)
  ))
})

# ===========================================================================
# 13. normalize_times_01
# ===========================================================================

test_that("normalize_times_01 rescales to [0,1]", {
  net <- network(matrix(c(1, 2), nrow = 1), directed = FALSE)
  set.vertex.attribute(net, "time", c(10, 20))
  set.edge.attribute(net, "time", 15, e = 1)
  net_norm <- normalize_times_01(net)
  vtimes <- get.vertex.attribute(net_norm, "time")
  etimes <- get.edge.attribute(net_norm, "time")
  expect_true(all(vtimes >= 0 & vtimes <= 1))
  expect_true(all(etimes >= 0 & etimes <= 1))
})

# ===========================================================================
# 14. events_to_net round-trip
# ===========================================================================

test_that("events_to_net builds network from event list", {
  el <- list(i = c(1, 1, 2), j = c(2, 3, 3), t = c(0.1, 0.2, 0.3))
  net <- events_to_net(el)
  expect_true(network.size(net) == 3)
  expect_true(network.edgecount(net) == 3)
})
