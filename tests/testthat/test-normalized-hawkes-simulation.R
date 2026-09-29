reference_hawkes_thinning <- function(mu, K, beta, window, history = numeric(0)) {
  t <- max(c(window[1L], history))
  events <- numeric(0)
  acceptance <- numeric(0)
  while (t < window[2L]) {
    bound <- mu + K * sum(exp(-beta * (t - c(history, events))))
    candidate <- t + rexp(1L, bound)
    if (candidate > window[2L]) break
    intensity <- mu + K * sum(exp(-beta * (candidate - c(history, events))))
    prob <- intensity / bound
    acceptance <- c(acceptance, prob)
    t <- candidate
    if (runif(1L) <= prob) events <- c(events, t)
  }
  list(times = events, acceptance = acceptance)
}

test_that("adaptive homogeneous recurrence agrees with a full-history reference", {
  params <- list(mu = 1.3, K = 0.7, beta_overall = 1.2, beta_edges = 0.3, m = 1.5)
  for (seed in c(2, 39, 113)) {
    set.seed(seed)
    reference <- reference_hawkes_thinning(params$mu, params$K, params$beta_overall, c(0, 8))
    set.seed(seed)
    out <- sim_hawkesNet(params, c(0, 8), PMF_mark_BA, cond_intensity)
    expect_equal(out$events$t, reference$times, tolerance = 1e-12)
    expect_equal(out$accept_probs, reference$acceptance, tolerance = 1e-12)
    expect_true(all(out$accept_probs >= 0 & out$accept_probs <= 1))
    expect_length(out$events$mark_density, out$events$n)
    expect_length(out$events$log_mark_density, out$events$n)
    expect_equal(get_times(out$net)$times, out$events$t)
    expect_equal(out$events$mark_density, exp(out$events$log_mark_density))
    direct <- vapply(out$events$t, function(t) {
      PMF_mark_BA(t, params, out$net)$log_mark_density
    }, numeric(1))
    expect_equal(out$events$log_mark_density, direct)
  }
})

test_that("ground times do not depend on mark draws or obsolete homogeneous envelope", {
  params <- list(mu = 2, K = 0.4, beta_overall = 1, beta_edges = 0.3, m = 1)
  set.seed(819)
  a <- sim_hawkesNet(params, c(0, 5), PMF_mark_BA, cond_intensity, mu_multiplier = 0.01)
  params$m <- 4
  params$beta_edges <- 2
  set.seed(819)
  b <- sim_hawkesNet(params, c(0, 5), PMF_mark_BA, cond_intensity, mu_multiplier = 100)
  expect_equal(a$events$t, b$events$t)
  expect_equal(a$accept_probs, b$accept_probs)
  expect_error(sim_hawkesNet(params, c(0, 1), PMF_mark_BA, cond_intensity,
                            joint_accept = TRUE), "joint_accept = FALSE")
  expect_error(sim_hawkesNet(params, c(0, 1), PMF_mark_BA, cond_intensity,
                            n_mark_sample = 2), "n_mark_sample = NULL")
})

test_that("conditional simulation initializes decayed seed excitation correctly", {
  params <- list(mu = 1, K = 0.7, beta_overall = 1.2, beta_edges = 0.3, m = 1.5)
  history <- c(-2, -0.5, 0)
  net <- network::network.initialize(3, directed = FALSE)
  network::set.vertex.attribute(net, "time", history)
  set.seed(133)
  reference <- reference_hawkes_thinning(params$mu, params$K, params$beta_overall,
                                        c(1, 5), history)
  set.seed(133)
  out <- sim_hawkesNet(params, c(1, 5), PMF_mark_BA, cond_intensity,
                       seed_net = net, seed_times = history)
  expect_equal(out$events$t, reference$times, tolerance = 1e-12)
  expect_equal(out$accept_probs, reference$acceptance, tolerance = 1e-12)
  expect_equal(network::network.size(net), 3)
  set.seed(133)
  derived <- sim_hawkesNet(params, c(1, 5), PMF_mark_BA, cond_intensity, seed_net = net)
  expect_equal(derived$events$t, out$events$t)
})

test_that("empty-history survival has no forced initial event", {
  params <- list(mu = 0.8, K = 0.7, beta_overall = 1.2, beta_edges = 0.3, m = 1.5)
  # A complete-update point mass avoids attachment RNG in this temporal check.
  point_mass_mark <- function(time, params, mark_filtration, ...) {
    out <- network::network.copy(mark_filtration)
    network::add.vertices(out, 1L)
    network::set.vertex.attribute(out, "time",
      c(network::get.vertex.attribute(mark_filtration, "time"), time))
    list(mark_sample = out, mark_sample_density = 1, log_mark_sample_density = 0)
  }
  set.seed(672)
  counts <- replicate(600, sim_hawkesNet(params, c(0, 0.5), point_mass_mark,
                                         cond_intensity)$events$n)
  expect_lt(abs(mean(counts == 0) - exp(-params$mu * 0.5)), 0.06)
  rate <- params$beta_overall - params$K
  expected <- params$mu * params$beta_overall / rate * 0.5 -
    params$mu * params$K / rate^2 * (-expm1(-rate * 0.5))
  expect_lt(abs(mean(counts) - expected), 5 * sd(counts) / sqrt(length(counts)))
})

test_that("Poisson timing and no-attempt marks are valid boundary models", {
  params <- list(mu = 2, K = 0, beta_overall = 1, beta_edges = 0, m = 0)
  set.seed(302)
  reference <- reference_hawkes_thinning(2, 0, 1, c(0, 3))
  set.seed(302)
  out <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity)
  expect_equal(out$events$t, reference$times)
  expect_true(all(out$accept_probs == 1))
  expect_equal(network::network.edgecount(out$net), 0)
  expect_equal(network::network.size(out$net), out$events$n)
  expect_true(all(out$events$mark_density == 1))
  ll <- loglik_hawkesNet(params, c(0, 3), out$net, PMF_mark_BA, verbose = FALSE)
  expect_equal(ll$loglik, out$events$n * log(2) - 2 * 3)
})

test_that("empty observations have the survival likelihood for both mark models", {
  net <- network::network.initialize(0, directed = FALSE)
  network::set.vertex.attribute(net, "time", numeric(0))
  params <- list(mu = 2, K = .4, beta_overall = 1, beta_edges = .3,
                 m = 1, node_lambda = 1, CS_params = c(0, .2))
  for (pmf in list(PMF_mark_BA, PMF_mark_CS)) {
    q <- loglik_hawkesNet(params, c(7, 10), net, pmf,
                          formula_RHS = "edges + triangles")
    expect_equal(q$loglik, -6)
    expect_length(q$intens_funcs, 0)
    qi <- loglik_hawkesNet(params, c(7, 10), net, pmf,
                           mu_vec = numeric(0), integral_bg = 4.5,
                           formula_RHS = "edges + triangles")
    expect_equal(qi$loglik, -4.5)
  }
})

test_that("the likelihood uses absolute observation endpoints consistently", {
  params <- list(mu = 2, K = .4, beta_overall = 1, beta_edges = .3, m = 1)
  net <- network::network.initialize(3, directed = FALSE)
  network::set.vertex.attribute(net, "time", c(.2, .7, 1.1))
  # Isolated one-node arrivals have valid positive probabilities under BA.
  q <- loglik_hawkesNet(params, c(0, 2), net, PMF_mark_BA)
  shifted <- network::network.copy(net)
  network::set.vertex.attribute(shifted, "time", c(.2, .7, 1.1) + 9)
  qs <- loglik_hawkesNet(params, c(9, 11), shifted, PMF_mark_BA)
  expect_equal(qs$loglik, q$loglik, tolerance = 1e-11)
  expect_error(loglik_hawkesNet(params, c(0, 2), shifted, PMF_mark_BA), "outside time window")
})
