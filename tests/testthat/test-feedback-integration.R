# End-to-end checks use a complete-update fixture with a two-edge arrival.
# At 1.1, vertex 3 and both edges 1--3, 2--3 arrive together: one mark.
feedback_integration_fixture <- function(shift = 0) {
  net <- network::network.initialize(4, directed = FALSE)
  network::set.vertex.attribute(net, "time", c(.2, .7, 1.1, 1.8) + shift)
  network::add.edges(net, c(1, 1, 2, 3), c(2, 3, 3, 4),
    names.eval = "time", vals.eval = c(.7, 1.1, 1.1, 1.8) + shift)
  net
}

# Independent feature values obtained directly from the four small graphs.
feedback_integration_features <- function(model, feature) {
  if (feature == "degree") {
    evolving <- rbind(c(0, 0, 0, 0), c(1/2, 1/2, 0, 0),
                      c(2/3, 2/3, 2/3, 0), c(2/3, 2/3, 7/10, 2/3))
    fixed <- c(0, 1/2, 2/3, 2/3)
  } else {
    evolving <- rbind(c(0, 0, 0, 0), c(0, 0, 0, 0),
                      c(0, 1, 1, 0), c(0, 1, 1, 0))
    fixed <- c(0, 0, 1, 0)
  }
  if (model == "M1") matrix(rep(fixed, each = 4), nrow = 4) else evolving
}

feedback_integration_reference <- function(params, model, feature, end,
                                           start = 0, shift = 0) {
  times <- c(.2, .7, 1.1, 1.8) + shift
  features <- feedback_integration_features(model, feature)
  alpha <- if (model == "M0") matrix(1, 4, 4) else if (feature == "degree") {
    exp(params$feedback_gamma * features)
  } else 2 * stats::plogis(params$feedback_gamma * features)
  intensity <- function(t) {
    k <- sum(times < t)
    if (!k) return(params$mu)
    j <- seq_len(k)
    params$mu + params$K * sum(alpha[k, j] * exp(-params$beta_overall * (t - times[j])))
  }
  cuts <- sort(unique(c(start, times[times > start & times < end], end)))
  integral <- params$mu * (end - start)
  for (r in seq_len(length(cuts) - 1L)) {
    left <- cuts[r]
    right <- cuts[r + 1L]
    k <- sum(times <= left)
    if (!k) next
    j <- seq_len(k)
    integral <- integral + params$K / params$beta_overall *
      sum(alpha[k, j] * (exp(-params$beta_overall * (left - times[j])) -
                          exp(-params$beta_overall * (right - times[j]))))
  }
  list(intensity = intensity, integral = integral)
}

test_that("M0 and zero-feedback simulation preserve the existing seeded path", {
  params <- list(mu = 2, K = .3, beta_overall = 1.2, beta_edges = .4,
                 m = 1.5, feedback_gamma = 0)
  set.seed(185)
  old <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
                      verbose = FALSE)
  for (model in c("M0", "M1", "M2")) {
    set.seed(185)
    out <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
                        timing = hawkes_timing(model), verbose = FALSE)
    expect_equal(out$events, old$events)
    expect_equal(out$accept_probs, old$accept_probs)
    expect_equal(network::as.edgelist(out$net), network::as.edgelist(old$net))
  }
})

test_that("complete-update content changes future simulated waiting times", {
  # These point masses consume no random draws. Timing differences therefore
  # arise from network feedback, rather than a shifted mark RNG stream.
  point_mass_update <- function(dense) {
    force(dense)
    function(time, params, mark_filtration, ...) {
      out <- network::network.copy(mark_filtration)
      n <- network::network.size(out)
      network::add.vertices(out, 1L)
      network::set.vertex.attribute(out, "time",
        c(network::get.vertex.attribute(mark_filtration, "time"), time))
      if (dense && n > 0) {
        network::add.edges(out, rep(n + 1L, n), seq_len(n),
                          names.eval = "time", vals.eval = rep(time, n))
      }
      list(mark_sample = out, mark_sample_density = 1, log_mark_sample_density = 0,
           log_mark_density = 0, log_density_func = function(params) 0)
    }
  }
  params <- list(mu = 2, K = .4, beta_overall = 1.2, beta_edges = .4,
                 m = 1, feedback_gamma = .8)
  paths <- list()
  for (model in c("M0", "M1", "M2")) {
    timing <- hawkes_timing(model, "degree")
    set.seed(202)
    isolated <- sim_hawkesNet(params, c(0, 4), point_mass_update(FALSE),
                             cond_intensity, timing = timing, verbose = FALSE)
    set.seed(202)
    dense <- sim_hawkesNet(params, c(0, 4), point_mass_update(TRUE),
                          cond_intensity, timing = timing, verbose = FALSE)
    expect_gte(dense$events$n, 3)
    expect_equal(dense$events$t[1:2], isolated$events$t[1:2])
    if (model == "M0") expect_equal(dense$events$t, isolated$events$t)
    else expect_false(isTRUE(all.equal(dense$events$t, isolated$events$t)))
    paths[[model]] <- dense$events$t
  }
  expect_false(isTRUE(all.equal(paths$M1, paths$M2)))
})

test_that("feedback simulation retains the complete BA update probability", {
  params <- list(mu = 2, K = .2, beta_overall = 1.1, beta_edges = .3,
                 m = 2.5, feedback_gamma = .6)
  for (model in c("M1", "M2")) {
    set.seed(813)
    out <- sim_hawkesNet(params, c(0, 4), PMF_mark_BA, cond_intensity,
                        timing = hawkes_timing(model), verbose = FALSE)
    expect_gte(out$events$n, 3)
    direct <- vapply(out$events$t, function(t) {
      PMF_mark_BA(t, params, out$net)$log_mark_density
    }, numeric(1))
    expect_equal(out$events$log_mark_density, direct, tolerance = 1e-11)
    expect_equal(out$events$mark_density, exp(direct), tolerance = 1e-11)
    expect_equal(get_times(out$net)$times, out$events$t)
    expect_true(all(out$accept_probs >= 0 & out$accept_probs <= 1))
  }
})

test_that("full likelihood and residuals use whole-mark interval feedback", {
  params <- list(mu = .8, K = .4, beta_overall = 1.3, beta_edges = .2,
                 m = 1.6, feedback_gamma = .7)
  net <- feedback_integration_fixture()
  times <- c(.2, .7, 1.1, 1.8)
  direct_logq <- vapply(times, function(t) PMF_mark_BA(t, params, net)$log_mark_density,
                       numeric(1))
  expect_length(get_times(net)$times, 4) # eight node/edge additions at four events
  for (model in c("M0", "M1", "M2")) {
    for (feature in c("degree", "triangles")) {
      timing <- hawkes_timing(model, feature)
      ref <- feedback_integration_reference(params, model, feature, end = 2.6)
      expected <- sum(log(vapply(times, ref$intensity, numeric(1)))) +
        sum(direct_logq) - ref$integral
      ll <- suppressMessages(loglik_hawkesNet(params, c(0, 2.6), net,
                            PMF_mark_BA, timing = timing, combine_intensity = TRUE))
      expect_equal(ll$loglik, expected, tolerance = 1e-10)
      comp <- compensators_hawkesNet(params, c(0, 2.6), net, timing = timing)
      expected_comp <- vapply(times, function(t) {
        feedback_integration_reference(params, model, feature, end = t)$integral
      }, numeric(1))
      expect_equal(comp, expected_comp, tolerance = 1e-10)
      expected_ks <- stats::ks.test(-expm1(-diff(c(0, comp))), "punif")$p.value
      expect_equal(ks_test_pval_hawkesNet(params, c(0, 2.6), net, timing = timing),
                   expected_ks, tolerance = 1e-10)
      shifted <- suppressMessages(loglik_hawkesNet(params, c(9, 11.6),
        feedback_integration_fixture(9), PMF_mark_BA, timing = timing))
      expect_equal(shifted$loglik, expected, tolerance = 1e-10)
      expect_equal(compensators_hawkesNet(params, c(9, 11.6),
        feedback_integration_fixture(9), timing = timing), comp, tolerance = 1e-10)
    }
  }
})

test_that("feedback intensity caching updates gamma and an inhomogeneous baseline", {
  params <- list(mu = .8, K = .4, beta_overall = 1.3, beta_edges = .2,
                 m = 1.6, feedback_gamma = .7)
  net <- feedback_integration_fixture()
  timing <- hawkes_timing("M2", "triangles")
  cached <- suppressMessages(loglik_hawkesNet(params, c(0, 2.6), net,
    PMF_mark_BA, timing = timing, mu_vec = c(.5, .7, .9, 1.1), integral_bg = 2))
  for (gamma in c(-.9, 0, 1.2)) {
    trial <- params
    trial$feedback_gamma <- gamma
    direct <- suppressMessages(loglik_hawkesNet(trial, c(0, 2.6), net,
      PMF_mark_BA, timing = timing, mu_vec = c(.5, .7, .9, 1.1), integral_bg = 2))
    from_cache <- loglik_hawkesNet(trial, c(0, 2.6), net, PMF_mark_BA,
      timing = timing, mu_vec = c(.5, .7, .9, 1.1), integral_bg = 2,
      intens_funcs = cached$intens_funcs)
    expect_equal(from_cache$loglik, direct$loglik, tolerance = 1e-11)
    reference <- feedback_integration_reference(trial, "M2", "triangles", end = 2.6)
    times <- c(.2, .7, 1.1, 1.8)
    ground <- vapply(times, reference$intensity, numeric(1)) - trial$mu + c(.5, .7, .9, 1.1)
    logq <- vapply(times, function(t) PMF_mark_BA(t, trial, net)$log_mark_density, numeric(1))
    expect_equal(direct$loglik, sum(log(ground)) + sum(logq) -
      (reference$integral - trial$mu * 2.6 + 2), tolerance = 1e-10)
  }
})

test_that("cached and uncached fits agree while estimating feedback gamma", {
  params <- list(mu = .2, K = .3, beta_overall = 1.1, beta_edges = .2,
                 m = 1.6, feedback_gamma = .1)
  net <- feedback_integration_fixture(9)
  fixed <- setdiff(names(params), "feedback_gamma")
  for (model in c("M1", "M2")) {
    timing <- hawkes_timing(model, "degree")
    fit_one <- function(cache) suppressMessages(fit_hawkesNet(params, c(9, 13.5),
      net, PMF_mark_BA, timing = timing, fixed_params = fixed,
      method = "L-BFGS-B", maxit = 80, cache_intensity = cache,
      combine_intensity = TRUE, get_hessian = FALSE, verbose = FALSE))
    cached <- fit_one(TRUE)
    direct <- fit_one(FALSE)
    expect_equal(cached$fit$par, direct$fit$par, tolerance = 1e-5)
    expect_equal(cached$fit$value, direct$fit$value, tolerance = 1e-8)
    expect_named(cached$fit$par, "feedback_gamma")
    expect_equal(cached$params$feedback_gamma, unname(cached$fit$par["feedback_gamma"]))
    checked <- suppressMessages(loglik_hawkesNet(cached$params, c(9, 13.5),
                               net, PMF_mark_BA, timing = timing))
    expect_equal(cached$fit$value, checked$loglik, tolerance = 1e-8)
    merged <- merge_fit_params(cached$fit$par, params, fixed)
    expect_equal(merged$feedback_gamma, cached$params$feedback_gamma)
  }
})

test_that("feedback simulation and caches retain both complete-update CS laws", {
  params <- list(mu = 3, K = .2, beta_overall = 1.2, m = 1.3, beta_edges = .4,
                 node_lambda = .7, CS_params = c(0, .6), feedback_gamma = .5)
  for (mode in c("independent", "size_conditional")) {
    timing <- if (mode == "independent") hawkes_timing("M1", "degree") else
      hawkes_timing("M2", "triangles")
    options <- list(cs_mode = mode, formula_RHS = "edges + triangles", truncation = 3,
                    timing = timing, verbose = FALSE)
    set.seed(928)
    sim <- do.call(sim_hawkesNet, c(list(params = params, time_window = c(0, 2),
      PMF_mark = PMF_mark_CS, cond_intensity = cond_intensity), options))
    expect_gte(sim$events$n, 3)
    logq <- vapply(sim$events$t, function(t) PMF_mark_CS(t, params, sim$net,
      cs_mode = mode, formula_RHS = "edges + triangles", truncation = 3)$log_mark_density,
      numeric(1))
    expect_equal(sim$events$log_mark_density, logq, tolerance = 1e-10)
    common <- c(list(time_window = c(0, 2), mark_filtration = sim$net,
                     PMF_mark = PMF_mark_CS, combine_intensity = TRUE), options)
    initial <- suppressMessages(do.call(loglik_hawkesNet, c(list(params = params), common)))
    trial <- params
    trial$feedback_gamma <- -.7
    fresh <- suppressMessages(do.call(loglik_hawkesNet, c(list(params = trial), common)))
    cached <- do.call(loglik_hawkesNet, c(list(params = trial,
      intens_funcs = initial$intens_funcs), common))
    expect_true(is.finite(fresh$loglik))
    expect_equal(cached$loglik, fresh$loglik, tolerance = 1e-10)
  }
})

test_that("inhomogeneous feedback simulation uses an explicit background bound", {
  params <- list(mu = 2, K = .2, beta_overall = 1.2, beta_edges = .3,
                 m = 2, feedback_gamma = .5)
  timing <- hawkes_timing("M2", "triangles")
  background <- list(mu_fit = list(mu_fun = function(t) rep(2, length(t)), mu_bound = 2))
  set.seed(581)
  homogeneous <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
                              timing = timing)
  set.seed(581)
  inhomogeneous <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
                               timing = timing, inhom_bg = background)
  expect_equal(inhomogeneous$events, homogeneous$events)
  expect_equal(inhomogeneous$accept_probs, homogeneous$accept_probs)
  background$mu_fit$mu_bound <- NULL
  expect_error(sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
                            timing = timing, inhom_bg = background), "mu_bound")
})

test_that("GOF simulations retain fitted feedback and explicit whole-update choices", {
  captured <- list()
  testthat::local_mocked_bindings(
    sim_hawkesNet = function(...) {
      captured[[length(captured) + 1L]] <<- list(...)
      # Return no network so this propagation test avoids unrelated statistics.
      list(net = NULL, error = "intentional test stub")
    },
    safe_parallel_lapply = function(X, FUN, ...) {
      # Exercise the same serialization boundary used for PSOCK closures.
      lapply(X, unserialize(serialize(FUN, NULL)))
    },
    .package = "hawkesNet"
  )
  params <- list(mu = 2, K = 1.4, beta_overall = 1.2, beta_edges = 0,
                 m = 2, node_lambda = 0, feedback_gamma = 0)
  fitted_timing <- hawkes_timing("M2", "triangles")
  fit <- list(fit = list(par = c(feedback_gamma = -.8)), timing = fitted_timing)
  net <- feedback_integration_fixture()
  out <- suppressMessages(gof(fit, net, params, PMF_mark_CS, cond_intensity,
    time_window = c(0, 2), n_sim = 1, cores = 1, verbose = FALSE,
    cs_mode = "size_conditional", max_candidates = 8L))
  expect_length(captured, 1)
  expect_identical(captured[[1]]$timing, fitted_timing)
  expect_equal(unname(captured[[1]]$params$feedback_gamma), -.8)
  expect_identical(captured[[1]]$cs_mode, "size_conditional")
  expect_identical(captured[[1]]$max_candidates, 8L)
  expect_identical(out$timing, fitted_timing)
  expect_identical(captured[[1]]$params$K, 1.4)
  expect_identical(captured[[1]]$params$beta_edges, 0)
  expect_identical(captured[[1]]$params$node_lambda, 0)

  override <- hawkes_timing("M1", "degree", scale = 4)
  params$node_lambda <- NULL
  suppressMessages(gof(fit, net, params, PMF_mark_BA, cond_intensity,
    time_window = c(0, 2), n_sim = 1, cores = 1, verbose = FALSE, timing = override))
  expect_identical(captured[[2]]$timing, override)
  expect_null(captured[[2]]$params$node_lambda)
  fit$timing <- NULL
  params$K <- 0
  suppressMessages(gof(fit, net, params, PMF_mark_BA, cond_intensity,
    time_window = c(0, 2), n_sim = 1, cores = 1, verbose = FALSE))
  expect_identical(captured[[3]]$timing, hawkes_timing())
  expect_identical(captured[[3]]$params$K, 0)

  background <- list(mu_fit = list(mu_fun = function(t) rep(2, length(t)), mu_bound = 2))
  params$mu <- NULL
  suppressMessages(gof(fit, net, params, PMF_mark_BA, cond_intensity,
    time_window = c(0, 2), n_sim = 1, cores = 1, verbose = FALSE,
    timing = fitted_timing, inhom_bg = background))
  expect_identical(captured[[4]]$params$mu, 1)
  params$mu <- 2.3
  suppressMessages(gof(fit, net, params, PMF_mark_BA, cond_intensity,
    time_window = c(0, 2), n_sim = 1, cores = 1, verbose = FALSE,
    timing = fitted_timing, inhom_bg = background))
  expect_identical(captured[[5]]$params$mu, 2.3)
  expect_error(gof(fit, net, params, PMF_mark_BA, cond_intensity,
    n_sim = 1, cores = 1, verbose = FALSE, timing = "M2"), "hawkes_timing")
  params$K <- -1
  expect_error(gof(fit, net, params, PMF_mark_BA, cond_intensity,
    n_sim = 1, cores = 1, verbose = FALSE), "Invalid point process parameters")
})

test_that("inhomogeneous feedback intensity and likelihood need no scalar mu", {
  params <- list(K = .4, beta_overall = 1.3, beta_edges = .2,
                 m = 1.6, feedback_gamma = .7)
  with_mu <- params
  with_mu$mu <- .8
  net <- feedback_integration_fixture()
  for (model in c("M0", "M1", "M2")) {
    timing <- hawkes_timing(model, "degree")
    bare <- cond_intensity_inhom(net, 1.8, net, PMF_mark_BA, params,
                                mu_at_t = .9, timing = timing)
    reference <- cond_intensity_inhom(net, 1.8, net, PMF_mark_BA, with_mu,
                                     mu_at_t = .9, timing = timing)
    expect_equal(bare$result, reference$result, tolerance = 1e-11)
    expect_equal(bare$func(params), reference$result, tolerance = 1e-11)
    bare_ll <- suppressMessages(loglik_hawkesNet(params, c(0, 2.6), net, PMF_mark_BA,
      timing = timing, mu_vec = c(.5, .7, .9, 1.1), integral_bg = 2))
    with_mu_ll <- suppressMessages(loglik_hawkesNet(with_mu, c(0, 2.6), net, PMF_mark_BA,
      timing = timing, mu_vec = c(.5, .7, .9, 1.1), integral_bg = 2))
    expect_equal(bare_ll$loglik, with_mu_ll$loglik, tolerance = 1e-11)
    cached <- loglik_hawkesNet(params, c(0, 2.6), net, PMF_mark_BA,
      timing = timing, mu_vec = c(.5, .7, .9, 1.1), integral_bg = 2,
      intens_funcs = with_mu_ll$intens_funcs)
    expect_equal(cached$loglik, with_mu_ll$loglik, tolerance = 1e-11)
  }
})

test_that("feedback seed times must describe exactly the stored complete updates", {
  params <- list(mu = 2, K = .2, beta_overall = 1.2, beta_edges = .3, m = 2)
  history <- feedback_integration_fixture()
  for (model in c("M1", "M2")) {
    for (gamma in c(0, .5)) {
      params$feedback_gamma <- gamma
      expect_error(sim_hawkesNet(params, c(2, 3), PMF_mark_BA, cond_intensity,
        timing = hawkes_timing(model), seed_net = history,
        seed_times = c(.2, .7, 1.1)), "seed_times must match")
    }
  }
})
