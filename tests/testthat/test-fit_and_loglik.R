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

test_that("CS fit with p_scale (formula names) does not error on parscale NA", {
  # Regression: parscale with formula names (edges, triangles, star.2, star.3) but
  # flat_par uses CS_params1,2,3,4 -> parscale[names(flat_par)] gave NA -> optim error.
  params_true <- list(mu = 10, beta_overall = 2, K = 0.5, beta_edges = 1, node_lambda = 1,
                      CS_params = c(-6.7, 2, 0.1, -0.1))
  set.seed(1)
  sim <- sim_hawkesNet(params = params_true, time_window = c(0, 5),
                      PMF_mark = PMF_mark_CS, cond_intensity = cond_intensity,
                      hashed_edges = TRUE, verbose = FALSE, truncation = 500L,
                      formula_RHS = "edges + triangles + star(c(2,3))",
                      mark_decay = "node_entrance", growth_only = FALSE)
  skip_if(network::network.edgecount(sim$net) < 5, "Need at least 5 edges")
  n_nodes <- network::network.size(sim$net)
  params_init <- list(mu = 10, beta_overall = 1, K = 0.5, beta_edges = 1, node_lambda = 1,
                     CS_params = c(-10, 0, 0, 0))
  p_scale <- c(mu = 1, beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.1,
               edges = 1, triangles = 0.1, star.2 = 0.1, star.3 = 0.1)
  suppressMessages({
    fit <- fit_hawkesNet(params_init = params_init,
                         time_window = c(0, 5),
                         mark_filtration = sim$net,
                         PMF_mark = PMF_mark_CS,
                         formula_RHS = "edges + triangles + star(c(2,3))",
                         maxit = 100,
                         trace = 0,
                         truncation = n_nodes,
                         mark_decay = "node_entrance",
                         growth_only = FALSE,
                         fixed_params = c("K"),
                         parscale = p_scale,
                         cores = 1,
                         cache_intensity = TRUE,
                         combine_intensity = TRUE,
                         verbose = FALSE)
  })
  expect_type(fit, "list")
  expect_true("fit" %in% names(fit))
  expect_true(all(is.finite(fit$fit$par)))
})

test_that("CS model sim+fit at T=5 converges (truncation = network size)", {
  # Same setup as simulation_study_CS consistency study - package should easily fit these.
  params_true <- list(mu = 10, beta_overall = 2, K = 0.5, beta_edges = 1, node_lambda = 1,
                      CS_params = c(-6.7, 2, 0.1, -0.1))
  set.seed(1)
  sim <- sim_hawkesNet(
    params = params_true,
    time_window = c(0, 5),
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    hashed_edges = TRUE,
    verbose = FALSE,
    truncation = 500L,
    formula_RHS = "edges + triangles + star(c(2,3))",
    mark_decay = "node_entrance",
    growth_only = FALSE
  )
  skip_if(network::network.edgecount(sim$net) < 5, "Need at least 5 edges for CS fit test")
  n_nodes <- network::network.size(sim$net)
  # Init near true with small perturbation (as in consistency study)
  set.seed(2)
  params_init <- list(
    mu = params_true$mu,
    beta_overall = max(0.1, params_true$beta_overall * exp(rnorm(1, 0, 0.2))),
    K = params_true$K,
    beta_edges = max(0.1, params_true$beta_edges * exp(rnorm(1, 0, 0.2))),
    node_lambda = max(0.1, params_true$node_lambda * exp(rnorm(1, 0, 0.2))),
    CS_params = params_true$CS_params + rnorm(4, 0, 0.3)
  )
  params_init$CS_params[!is.finite(params_init$CS_params)] <- params_true$CS_params[!is.finite(params_init$CS_params)]
  suppressMessages({
    fit <- fit_hawkesNet(
      params_init = params_init,
      time_window = c(0, 5),
      mark_filtration = sim$net,
      PMF_mark = PMF_mark_CS,
      formula_RHS = "edges + triangles + star(c(2,3))",
      maxit = 1500,
      trace = 0,
      truncation = n_nodes,
      mark_decay = "node_entrance",
      growth_only = FALSE,
      fixed_params = c("K"),
      method = "Nelder-Mead",
      cores = 1,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      verbose = FALSE
    )
  })
  expect_type(fit, "list")
  expect_true("fit" %in% names(fit))
  expect_true(all(is.finite(fit$fit$par)), info = "All fitted params should be finite")
  expect_equal(fit$fit$convergence, 0, info = "CS fit at T=5 should converge (truncation = network size)")
})
