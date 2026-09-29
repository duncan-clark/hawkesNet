m0_underflow_fixture <- function() {
  net <- network::network.initialize(3, directed = FALSE)
  network::set.vertex.attribute(net, "time", c(.2, .4, .8))
  list(net = net, params = list(mu = 2, K = 0, beta_overall = 1,
                               beta_edges = 0, m = 1000))
}

test_that("M0 scores rare whole marks without probability floors", {
  x <- m0_underflow_fixture()
  expected <- 3 * log(2) - 2 * x$params$m - 2
  for (combine in c(FALSE, TRUE)) {
    direct <- suppressMessages(loglik_hawkesNet(x$params, c(0, 1), x$net,
      PMF_mark_BA, combine_intensity = combine))
    expect_equal(direct$loglik, expected, tolerance = 1e-10)
    expect_true(is.function(attr(direct$intens_funcs, "log_intensities")))
    changed <- x$params; changed$m <- 1200; changed$K <- .3
    cached <- loglik_hawkesNet(changed, c(0, 1), x$net, PMF_mark_BA,
                              intens_funcs = direct$intens_funcs)
    times <- c(.2, .4, .8)
    ground <- vapply(times, function(t) 2 + .3 * sum(exp(-(t-times[times<t]))), 0)
    reference <- sum(log(ground)) - 2400 - 2 - .3 * sum(-expm1(-(1-times)))
    expect_equal(cached$loglik, reference, tolerance = 1e-10)
    inh <- suppressMessages(loglik_hawkesNet(x$params, c(0, 1), x$net,
      PMF_mark_BA, mu_vec = c(1, 2, 3), integral_bg = 2,
      combine_intensity = combine))
    expect_equal(inh$loglik, log(6) - 2000 - 2, tolerance = 1e-10)
  }
})

test_that("M0 reports impossible observed marks as minus infinity", {
  x <- m0_underflow_fixture()
  network::add.edges(x$net, 1, 2, names.eval = "time", vals.eval = .4)
  x$params$m <- 0
  direct <- suppressMessages(loglik_hawkesNet(x$params, c(0, 1), x$net, PMF_mark_BA))
  expect_identical(direct$loglik, -Inf)
  cached <- loglik_hawkesNet(x$params, c(0, 1), x$net, PMF_mark_BA,
                            intens_funcs = direct$intens_funcs)
  expect_identical(cached$loglik, -Inf)
})

test_that("cached fitting escapes underflow and agrees with direct likelihood", {
  x <- m0_underflow_fixture()
  fit <- suppressMessages(fit_hawkesNet(x$params, c(0, 1), x$net, PMF_mark_BA,
    fixed_params = c("mu", "K", "beta_overall", "beta_edges"),
    method = "L-BFGS-B", maxit = 100, get_hessian = FALSE,
    verbose = FALSE, combine_intensity = TRUE))
  expect_lt(fit$params$m, 1)
  expect_equal(fit$fit$convergence, 0)
  direct <- suppressMessages(loglik_hawkesNet(fit$params, c(0, 1), x$net, PMF_mark_BA))
  expect_equal(fit$fit$value, direct$loglik, tolerance = 1e-8)
})

test_that("M0 logarithmic recurrence matches direct history sums", {
  times <- c(.2, .4, 1.1, 3, 8)
  marks <- lapply(seq_along(times), function(i) { force(i); function(p) -i * p$m })
  cache <- .m0_log_intensity_cache(marks, times)
  for (beta in c(.01, 1, 100)) {
    p <- list(mu = .7, K = 1.2, beta_overall = beta, m = 250)
    expected <- vapply(seq_along(times), function(i) {
      -i * p$m + log(p$mu + p$K * sum(exp(-beta * (times[i]-times[times<times[i]]))))
    }, 0)
    expect_equal(cache(p), expected, tolerance = 1e-12)
  }
})
