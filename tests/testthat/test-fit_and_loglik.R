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

# The former fixtures used the historical no-m CS law and skipped when its
# simulator failed. They did not validate a normalized model. These fixtures
# explicitly exercise CS-1 with one fixed candidate rule in simulation and fit;
# an error now fails the regression rather than silently skipping it.
make_cs1_fit_fixture <- function() {
  params <- list(mu=10,beta_overall=2,K=.5,beta_edges=.5,m=1.4,
                 node_lambda=.8,CS_params=c(0,.6,.1,-.03))
  options <- list(cs_mode="independent",truncation=4L,
    formula_RHS="edges + triangles + star(c(2,3))",
    mark_decay="node_entrance",growth_only=FALSE)
  set.seed(42)
  sim <- do.call(sim_hawkesNet,c(list(params=params,time_window=c(0,5),
    PMF_mark=PMF_mark_CS,cond_intensity=cond_intensity,verbose=FALSE),options))
  list(params=params,options=options,sim=sim)
}

test_that("CS-1 fit with formula-name parscale maps free structural coefficients", {
  # Regression: formula names (edges, triangles, star.2, star.3) must map to
  # CS_params1,2,3,4 after fixed parameters are removed, without introducing NA.
  fixture <- make_cs1_fit_fixture()
  expect_gte(network::network.edgecount(fixture$sim$net),5)
  initial <- fixture$params; initial$m <- .9; initial$CS_params[2:4] <- 0
  p_scale <- c(mu=1,beta_overall=.1,beta_edges=.1,node_lambda=.1,m=1,
               edges=1,triangles=.1,star.2=.1,star.3=.1)
  common <- c(list(time_window=c(0,5),mark_filtration=fixture$sim$net,
                  PMF_mark=PMF_mark_CS,verbose=FALSE),fixture$options)
  initial_ll <- do.call(loglik_hawkesNet,c(list(params=initial),common))$loglik
  fit <- suppressMessages(do.call(fit_hawkesNet,c(list(params_init=initial,
    maxit=100,method="L-BFGS-B",get_hessian=FALSE,parscale=p_scale,
    fixed_params=c("mu","K","beta_overall","beta_edges","node_lambda","CS_params1"),
    cache_intensity=TRUE,combine_intensity=TRUE,cores=1),common)))
  expect_type(fit,"list")
  expect_true(all(is.finite(fit$fit$par)))
  expect_setequal(names(fit$fit$par),c("m","CS_params2","CS_params3","CS_params4"))
  expect_gte(fit$fit$value,initial_ll-1e-8)
  direct <- do.call(loglik_hawkesNet,c(list(params=fit$params),common))$loglik
  expect_equal(fit$fit$value,direct,tolerance=1e-8)
})

test_that("CS-1 simulation and cached fitting converge with a fixed candidate rule", {
  fixture <- make_cs1_fit_fixture()
  expect_gte(length(fixture$sim$events$t),20)
  initial <- fixture$params
  initial$mu <- 8; initial$m <- .9; initial$node_lambda <- 1.1
  initial$CS_params[2] <- .1
  # T=5 is an integration fixture, not a full-parameter consistency experiment.
  # Fit a baseline rate plus count, birth and triangle effects; hold the decay,
  # excitation and other structural coefficients at their specified values.
  fixed <- c("K","beta_overall","beta_edges","CS_params1","CS_params3","CS_params4")
  common <- c(list(time_window=c(0,5),mark_filtration=fixture$sim$net,
                  PMF_mark=PMF_mark_CS,verbose=FALSE),fixture$options)
  initial_ll <- do.call(loglik_hawkesNet,c(list(params=initial),common))$loglik
  fit <- suppressMessages(do.call(fit_hawkesNet,c(list(params_init=initial,
    maxit=300,method="L-BFGS-B",get_hessian=FALSE,fixed_params=fixed,
    cores=1,cache_intensity=TRUE,combine_intensity=TRUE),common)))
  expect_type(fit,"list")
  expect_true(all(is.finite(fit$fit$par)))
  expect_setequal(names(fit$fit$par),c("mu","m","node_lambda","CS_params2"))
  expect_equal(fit$fit$convergence,0)
  expect_gt(fit$fit$value,initial_ll)
  direct <- do.call(loglik_hawkesNet,c(list(params=fit$params),common))$loglik
  expect_equal(fit$fit$value,direct,tolerance=1e-8)
  expect_equal(fit$params$CS_params[c(1,3,4)],initial$CS_params[c(1,3,4)])
})
