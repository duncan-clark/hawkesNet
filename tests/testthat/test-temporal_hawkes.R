# Tests for temporal Hawkes process: simulation and fitting.
# Logic moved from the former if(FALSE) block in R/temporal_hawkes.R.

# Helper: general Hawkes simulator for both exponential and power-law kernels
# (used only in tests; package exports simulate_hawkes_branching for exp only).
simulate_hawkes_general <- function(kernel = c("exp", "powerlaw"),
                                   mu, K,
                                   beta = NULL,
                                   p = NULL, c = NULL,
                                   T_max, seed = NULL) {
  kernel <- match.arg(kernel)
  if (!is.null(seed)) set.seed(seed)
  N0 <- rpois(1, lambda = mu * T_max)
  if (N0 > 0) {
    events <- sort(runif(N0, min = 0, max = T_max))
  } else {
    events <- numeric(0)
  }
  i <- 1
  while (i <= length(events)) {
    parent_time <- events[i]
    num_children <- rpois(1, K)
    if (num_children > 0) {
      if (kernel == "exp") {
        offsets <- rexp(num_children, rate = beta)
      } else {
        u <- runif(num_children)
        offsets <- c * ((1 - u)^(-1 / (p - 1)) - 1)
      }
      child_times <- parent_time + offsets
      child_times <- child_times[child_times <= T_max]
      if (length(child_times) > 0) {
        events <- c(events, child_times)
      }
    }
    i <- i + 1
  }
  sort(events)
}

test_that("simulate_hawkes_branching returns valid event times", {
  set.seed(1)
  times <- simulate_hawkes_branching(mu = 0.2, K = 0.5, beta = 1, T = 10)
  expect_type(times, "double")
  expect_true(length(times) >= 0)
  expect_true(all(times >= 0 & times <= 10))
  expect_true(all(diff(sort(times)) >= 0))
})

test_that("fit_temporal_hawkes runs and returns expected structure for exp kernel", {
  set.seed(2)
  times <- simulate_hawkes_branching(mu = 2, K = 0.5, beta = 1.5, T = 20)
  expect_gte(length(times), 10)
  realiz <- data.frame(t = times)
  windowT <- c(0, 20)
  fit <- fit_temporal_hawkes(
    params_init = list(gamma = length(times) / 2, beta = 1, K = 0.2),
    realiz = realiz,
    windowT = windowT,
    kernel = "exp",
    trace = 0
  )
  expect_type(fit, "list")
  expect_true("par" %in% names(fit))
  expect_true("value" %in% names(fit))
  expect_length(fit$par, 3)
  expect_true(all(is.finite(fit$par)))
})

test_that("fit_temporal_hawkes runs for powerlaw kernel on simulated powerlaw data", {
  set.seed(3)
  times <- simulate_hawkes_general(
    "powerlaw", mu = 2, K = 0.5, p = 2.5, c = 0.5, T_max = 30
  )
  expect_gte(length(times), 15)
  realiz <- data.frame(t = times)
  windowT <- c(0, 30)
  fit <- fit_temporal_hawkes(
    params_init = list(gamma = length(times) / 2, p = 2, K = 0.2, c = 0.1),
    realiz = realiz,
    windowT = windowT,
    kernel = "powerlaw",
    trace = 0
  )
  expect_type(fit, "list")
  expect_true("par" %in% names(fit))
  expect_length(fit$par, 4)
  expect_true(all(is.finite(fit$par)))
})

test_that("loglik_temporal_hawkes returns finite value for valid inputs", {
  set.seed(4)
  times <- simulate_hawkes_branching(mu = 0.2, K = 0.4, beta = 1, T = 15)
  skip_if(length(times) < 5, "Too few events")
  realiz <- data.frame(t = times)
  windowT <- c(0, 15)
  dens <- density(times, from = windowT[1], to = windowT[2])
  mu_at_events <- approx(dens$x, dens$y, xout = times)$y
  kde_bg <- list(mu_vec = mu_at_events, total_int = 1)
  params <- c(gamma = 3, beta = 1, K = 0.4)
  ll <- loglik_temporal_hawkes(params, realiz, windowT, kde_bg = kde_bg, kernel = "exp")
  expect_length(ll, 1)
  expect_true(is.finite(ll))
})

test_that("temporal likelihood uses strict past times independent of row order", {
  times <- c(0.2, 0.7, 0.7, 1.4)
  bg <- list(mu_vec = rep(0.5, 4), total_int = 1)
  p <- c(3, 2, 0.4)
  expected <- sum(vapply(times, function(t)
    log(1.5 + 0.4 * sum(2 * exp(-2 * (t - times[times < t])))), numeric(1))) -
    3 - 0.4 * sum(1 - exp(-2 * (2 - times)))
  for (ord in list(1:4, c(4, 2, 1, 3))) {
    expect_equal(loglik_temporal_hawkes(p, data.frame(t = times[ord]), c(0, 2),
                                      kde_bg = bg), expected)
  }
  expect_error(loglik_temporal_hawkes(p, data.frame(t = -1), c(0, 2), kde_bg = bg),
               "inside windowT")
  expect_equal(loglik_temporal_hawkes(c(3, -2, .4), data.frame(t = times),
                                    c(0, 2), kde_bg = bg), -Inf)
})

test_that("KDE likelihood and compensator integrate the same window-normalized rate", {
  times <- c(0, .01, .03, .12, .2, .3, .6, .8, .9, 1)
  bg <- compute_kde_background(times, c(0, 1))
  rate <- stats::approxfun(bg$grid, bg$density)
  expect_equal(sum(diff(bg$grid) *
                     (head(bg$density, -1) + tail(bg$density, -1)) / 2), 1)
  expect_equal(bg$integral_at_events[c(1, length(times))], c(0, 1))
  for (i in c(4L, 7L)) {
    numerical <- stats::integrate(rate, 0, times[i], subdivisions = 2000,
                                 rel.tol = 1e-5)$value
    expect_equal(bg$integral_at_events[i], numerical, tolerance = 2e-5)
  }
  comp <- compensator_temporal_hawkes(c(2, 3, 0), data.frame(t = times),
                                      c(0, 1), use_kde = TRUE)
  expect_equal(comp, 2 * bg$integral_at_events)
  expect_equal(loglik_temporal_hawkes(c(2, 3, 0), data.frame(t = times),
                                    c(0, 1), kde_bg = bg),
               sum(log(2 * bg$mu_vec)) - tail(comp, 1))
})

test_that("kernel integrals retain very short-lag mass", {
  expect_equal(hawkes_kernel_cdf(1e-18, "exp", beta = 2), 2e-18,
               tolerance = 1e-30)
  expect_equal(hawkes_kernel_cdf(1e-18, "powerlaw", beta = 2, c = 1), 1e-18,
               tolerance = 1e-30)
  expect_warning(kernel <- hawkes_kernel_density(matrix(c(-10, 0, 1)),
                        "powerlaw", beta = 1.5, c = .2), NA)
  expect_equal(as.vector(kernel)[1:2], c(0, 0))
})

test_that("temporal fits report missing Wald information at a branching-ratio boundary", {
  set.seed(3)
  fit <- fit_temporal_hawkes(c(gamma = 20, beta = 1, K = .1),
                            data.frame(t = sort(runif(50, 0, 10))), c(0, 10))
  expect_equal(unname(fit$par["K"]), 1e-5, tolerance = 1e-8)
  expect_true(is.finite(fit$value))
  expect_true("K" %in% fit$boundary)
  expect_true(all(is.na(fit$hessian)))
  expect_true(all(is.na(fit$standard_errors)))
  expect_match(fit$hessian_status, "boundary")
})
