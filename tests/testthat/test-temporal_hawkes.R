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
  times <- simulate_hawkes_branching(mu = 0.2, K = 0.5, beta = 1.5, T = 20)
  skip_if(length(times) < 10, "Too few events for fit")
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
    "powerlaw", mu = 0.2, K = 0.5, p = 2.5, c = 0.5, T_max = 30
  )
  skip_if(length(times) < 15, "Too few events for powerlaw fit")
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
