# =============================================================================
# Integration checks for the redeveloped exact whole-update CS study.
# The four-node candidate window is an explicit model restriction: it makes
# enumeration exact (at most six candidates), and does not reproduce the old
# 100-node independent-edge study or establish statistical consistency.
# =============================================================================

# ---- Bounded whole-update configuration ----
SIM_STUDY_PARAMS <- list(
  mu = 10,
  beta_overall = 2,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  m = 1,
  CS_params = c(0, 0.4, 0.02, -0.01)
)
SIM_STUDY_FORMULA <- "edges + triangles + star(c(2,3))"
SIM_STUDY_TRUNCATION <- 4

# parscale: must match free flat_par names (CS_params1 edges is fixed)
SIM_STUDY_PSCALE <- c(
  mu = 1, beta_overall = 0.1, K = 0.1, beta_edges = 0.1, node_lambda = 0.5, m = 0.1,
  CS_params2 = 0.1, CS_params3 = 0.1, CS_params4 = 0.1
)

# Hold the edge tilt fixed to separate it from the reference attachment scale m.
# The joint likelihood is not generally flat in the edge coefficient.
SIM_STUDY_INIT <- list(
  mu = 10,
  beta_overall = 1,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  m = 1,
  CS_params = c(0, 0.2, 0, 0)
)
SIM_STUDY_FIXED <- c("CS_params1")

# Helper: simulate one network at a given time window
sim_one <- function(time_window = c(0, 5), seed = 42) {
  set.seed(seed)
  sim_hawkesNet(
    params = SIM_STUDY_PARAMS,
    time_window = time_window,
    PMF_mark = PMF_mark_CS,
    cond_intensity = cond_intensity,
    hashed_edges = TRUE,
    verbose = FALSE,
    mu_multiplier = 3,
    truncation = SIM_STUDY_TRUNCATION,
    formula_RHS = SIM_STUDY_FORMULA,
    mark_decay = "node_entrance",
    growth_only = FALSE
  )
}

# =============================================================================
# TEST 1: Simulation produces a valid network
# =============================================================================
test_that("CS simulation (T=5) produces valid network with edges and nodes", {
  sim <- sim_one(c(0, 5), seed = 42)
  expect_type(sim, "list")
  expect_true(!is.null(sim$net))
  expect_true(network::network.size(sim$net) > 0)
  expect_true(network::network.edgecount(sim$net) > 0)
  expect_true(length(sim$events$t) > 0)
  # All event times within window
  expect_true(all(sim$events$t >= 0 & sim$events$t <= 5))
  # Event times are recoverable via get_times()
  times_obj <- get_times(sim$net)
  expect_true(length(times_obj$times) > 0)
  expect_true(all(is.finite(times_obj$times)))
  expect_equal(times_obj$times, sim$events$t)
  direct_marks <- vapply(sim$events$t, function(t) {
    PMF_mark_CS(t, SIM_STUDY_PARAMS, sim$net, formula_RHS = SIM_STUDY_FORMULA,
                 truncation = SIM_STUDY_TRUNCATION)$log_mark_density
  }, numeric(1))
  expect_equal(sim$events$log_mark_density, direct_marks, tolerance = 1e-10)
})

# =============================================================================
# TEST 2: parscale mapping works with both CS_params and formula name conventions
# =============================================================================
test_that("parscale maps correctly to free flat_par (CS_params1 fixed)", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for fit test")

  # Flatten params_init as fit_hawkesNet does (remove CS_params1 since it's fixed)
  pi <- SIM_STUDY_INIT
  flat_par <- unlist(pi)
  flat_par <- flat_par[names(flat_par) != "CS_params1"]
  # Check that parscale[names(flat_par)] gives no NAs
  mapped <- SIM_STUDY_PSCALE[names(flat_par)]
  expect_true(all(!is.na(mapped)),
    info = paste("parscale NAs for:", paste(names(flat_par)[is.na(mapped)], collapse = ", ")))
  expect_true(all(is.finite(mapped)))
  expect_true(all(mapped > 0))
})

# =============================================================================
# TEST 3: Log-likelihood is finite at the main study init values
# =============================================================================
test_that("whole-update loglik is finite and cached evaluation matches fresh evaluation", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for loglik test")

  ll <- loglik_hawkesNet(
    params = SIM_STUDY_INIT,
    time_window = c(0, 5),
    mark_filtration = sim$net,
    PMF_mark = PMF_mark_CS,
    formula_RHS = SIM_STUDY_FORMULA,
    truncation = SIM_STUDY_TRUNCATION,
    mark_decay = "node_entrance",
    growth_only = FALSE,
    verbose = FALSE
  )
  expect_true(is.finite(ll$loglik),
    info = paste("loglik at init should be finite, got:", ll$loglik))
  changed <- SIM_STUDY_INIT
  changed$m <- 1.4
  changed$node_lambda <- 0.8
  changed$CS_params[2] <- 0.6
  cached <- loglik_hawkesNet(changed, c(0, 5), sim$net, PMF_mark_CS,
    intens_funcs = ll$intens_funcs, formula_RHS = SIM_STUDY_FORMULA,
    truncation = SIM_STUDY_TRUNCATION, verbose = FALSE)
  fresh <- loglik_hawkesNet(changed, c(0, 5), sim$net, PMF_mark_CS,
    formula_RHS = SIM_STUDY_FORMULA, truncation = SIM_STUDY_TRUNCATION,
    verbose = FALSE)
  expect_equal(cached$loglik, fresh$loglik, tolerance = 1e-10)
})

# =============================================================================
# TEST 4: Log-likelihood is finite even with OLD (distant) init values
# =============================================================================
test_that("loglik is finite at old distant init (CS_params = c(-10, 0, 0, 0))", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for loglik test")

  old_init <- list(mu = 10, beta_overall = 1, K = 0.5, beta_edges = 1,
                   node_lambda = 1, m = 1, CS_params = c(-10, 0, 0, 0))
  ll <- loglik_hawkesNet(
    params = old_init,
    time_window = c(0, 5),
    mark_filtration = sim$net,
    PMF_mark = PMF_mark_CS,
    formula_RHS = SIM_STUDY_FORMULA,
    truncation = SIM_STUDY_TRUNCATION,
    mark_decay = "node_entrance",
    growth_only = FALSE,
    verbose = FALSE
  )
  expect_true(is.finite(ll$loglik),
    info = paste("loglik at distant init should be finite, got:", ll$loglik))
})

# =============================================================================
# TEST 5: fit_hawkesNet completes without error (bounded joint study config)
# =============================================================================
test_that("CS fit with main study config completes without error", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for fit test")

  suppressMessages({
    fit <- fit_hawkesNet(
      params_init = SIM_STUDY_INIT,
      time_window = c(0, 5),
      mark_filtration = sim$net,
      PMF_mark = PMF_mark_CS,
      formula_RHS = SIM_STUDY_FORMULA,
      trace = 0,
      maxit = 200,
      truncation = SIM_STUDY_TRUNCATION,
      mark_decay = "node_entrance",
      growth_only = FALSE,
      fixed_params = SIM_STUDY_FIXED,
      method = "Nelder-Mead",
      parscale = SIM_STUDY_PSCALE,
      cores = 1,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      verbose = FALSE
    )
  })

  expect_type(fit, "list")
  expect_true("fit" %in% names(fit))
  expect_true(all(is.finite(fit$fit$par)),
    info = paste("All params should be finite:", paste(round(fit$fit$par, 4), collapse = ", ")))
  # CS_params1 (edges) should NOT be in fitted params (it's fixed)
  expect_false("CS_params1" %in% names(fit$fit$par))
  fitted <- merge_fit_params(fit$fit$par, SIM_STUDY_INIT, fixed_params = SIM_STUDY_FIXED)
  fresh <- loglik_hawkesNet(fitted, c(0, 5), sim$net, PMF_mark_CS,
    formula_RHS = SIM_STUDY_FORMULA, truncation = SIM_STUDY_TRUNCATION,
    verbose = FALSE)
  expect_equal(fit$fit$value, fresh$loglik, tolerance = 1e-8)
})

# =============================================================================
# TEST 6: fit_hawkesNet with formula-name-only parscale (old convention)
#          The guard in fit_hawkesNet should handle the NA parscale issue.
# =============================================================================
test_that("CS fit with formula-name-only parscale does not crash (guard test)", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for fit test")

  # Old parscale that uses ONLY formula names (no CS_params names)
  old_pscale <- c(mu = 1, beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.1,
                  edges = 1, triangles = 0.1, star.2 = 0.1, star.3 = 0.1)

  suppressMessages({
    fit <- fit_hawkesNet(
      params_init = SIM_STUDY_INIT,
      time_window = c(0, 5),
      mark_filtration = sim$net,
      PMF_mark = PMF_mark_CS,
      formula_RHS = SIM_STUDY_FORMULA,
      trace = 0,
      maxit = 50,
      truncation = SIM_STUDY_TRUNCATION,
      mark_decay = "node_entrance",
      growth_only = FALSE,
      fixed_params = SIM_STUDY_FIXED,
      method = "Nelder-Mead",
      parscale = old_pscale,
      cores = 1,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      verbose = FALSE
    )
  })

  expect_type(fit, "list")
  expect_true(all(is.finite(fit$fit$par)),
    info = "Fit should work even with formula-name-only parscale")
})

# =============================================================================
# TEST 7: Consistency study config (near-true init) also works
# =============================================================================
test_that("CS fit with a perturbed start returns finite parameters", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for fit test")

  set.seed(99)
  params_init_near <- list(
    mu = max(0.1, SIM_STUDY_PARAMS$mu * exp(rnorm(1, 0, 0.2))),
    beta_overall = max(0.1, SIM_STUDY_PARAMS$beta_overall * exp(rnorm(1, 0, 0.2))),
    K = min(0.99, max(0.01, SIM_STUDY_PARAMS$K * exp(rnorm(1, 0, 0.2)))),
    beta_edges = max(0.1, SIM_STUDY_PARAMS$beta_edges * exp(rnorm(1, 0, 0.2))),
    node_lambda = max(0.1, SIM_STUDY_PARAMS$node_lambda * exp(rnorm(1, 0, 0.2))),
    m = max(0.1, SIM_STUDY_PARAMS$m * exp(rnorm(1, 0, 0.2))),
    CS_params = c(0, SIM_STUDY_PARAMS$CS_params[-1] + rnorm(3, 0, 0.15))
  )

  suppressMessages({
    fit <- fit_hawkesNet(
      params_init = params_init_near,
      time_window = c(0, 5),
      mark_filtration = sim$net,
      PMF_mark = PMF_mark_CS,
      formula_RHS = SIM_STUDY_FORMULA,
      trace = 0,
      maxit = 500,
      truncation = SIM_STUDY_TRUNCATION,
      mark_decay = "node_entrance",
      growth_only = FALSE,
      fixed_params = SIM_STUDY_FIXED,
      method = "Nelder-Mead",
      parscale = SIM_STUDY_PSCALE,
      cores = 1,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      verbose = FALSE
    )
  })

  expect_type(fit, "list")
  expect_true(all(is.finite(fit$fit$par)))
  # This exercises the optimization path; finite parameters alone do not
  # establish convergence or consistency for a multi-parameter small sample.
})

# =============================================================================
# TEST 8: Temporal Hawkes fit works on CS simulation output (same as sim study)
# =============================================================================
test_that("fit_temporal_hawkes works on CS simulation output", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(length(sim$events$t) < 3, "Need events for temporal Hawkes fit")

  thf <- fit_temporal_hawkes(
    params_init = list(mu = 0.1, beta = 1, K = 0.1),
    realiz = data.frame(t = sim$events$t,
                        n = rep(sim$events$n, length(sim$events$t))),
    windowT = c(0, 5),
    trace = 0,
    maxit = 1000
  )

  expect_type(thf, "list")
  expect_true("par" %in% names(thf))
  expect_true(all(is.finite(thf$par)))
})

# =============================================================================
# TEST 9: Main + consistency inits give same finite results at T=10
# =============================================================================
test_that("Both main and consistency inits give finite fits at T=10", {
  sim <- sim_one(c(0, 10), seed = 123)
  skip_if(network::network.edgecount(sim$net) < 5, "Need edges for T=10 fit test")

  # Main study init
  suppressMessages({
    fit_main <- fit_hawkesNet(
      params_init = SIM_STUDY_INIT,
      time_window = c(0, 10),
      mark_filtration = sim$net,
      PMF_mark = PMF_mark_CS,
      formula_RHS = SIM_STUDY_FORMULA,
      trace = 0,
      maxit = 300,
      truncation = SIM_STUDY_TRUNCATION,
      mark_decay = "node_entrance",
      growth_only = FALSE,
      fixed_params = SIM_STUDY_FIXED,
      method = "Nelder-Mead",
      parscale = SIM_STUDY_PSCALE,
      cores = 1,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      verbose = FALSE
    )
  })

  expect_true(all(is.finite(fit_main$fit$par)),
    info = "Main study init should produce finite params at T=10")

  # Consistency study init (near-true)
  set.seed(7)
  params_init_cons <- list(
    mu = max(0.1, SIM_STUDY_PARAMS$mu * exp(rnorm(1, 0, 0.2))),
    beta_overall = max(0.1, SIM_STUDY_PARAMS$beta_overall * exp(rnorm(1, 0, 0.2))),
    K = min(0.99, max(0.01, SIM_STUDY_PARAMS$K * exp(rnorm(1, 0, 0.2)))),
    beta_edges = max(0.1, SIM_STUDY_PARAMS$beta_edges * exp(rnorm(1, 0, 0.2))),
    node_lambda = max(0.1, SIM_STUDY_PARAMS$node_lambda * exp(rnorm(1, 0, 0.2))),
    m = max(0.1, SIM_STUDY_PARAMS$m * exp(rnorm(1, 0, 0.2))),
    CS_params = c(0, SIM_STUDY_PARAMS$CS_params[-1] + rnorm(3, 0, 0.15))
  )

  suppressMessages({
    fit_cons <- fit_hawkesNet(
      params_init = params_init_cons,
      time_window = c(0, 10),
      mark_filtration = sim$net,
      PMF_mark = PMF_mark_CS,
      formula_RHS = SIM_STUDY_FORMULA,
      trace = 0,
      maxit = 300,
      truncation = SIM_STUDY_TRUNCATION,
      mark_decay = "node_entrance",
      growth_only = FALSE,
      fixed_params = SIM_STUDY_FIXED,
      method = "Nelder-Mead",
      parscale = SIM_STUDY_PSCALE,
      cores = 1,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      verbose = FALSE
    )
  })

  expect_true(all(is.finite(fit_cons$fit$par)),
    info = "Consistency study init should produce finite params at T=10")
})
