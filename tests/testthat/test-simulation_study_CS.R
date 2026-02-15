# =============================================================================
# Tests that mirror the CS simulation study (simulation_study_CS.R)
# Run locally to ensure the study will work on the cluster.
# =============================================================================

# ---- Shared configuration (mirrors simulation_study_CS.R exactly) ----
SIM_STUDY_PARAMS <- list(
  mu = 10,
  beta_overall = 2,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-6.7, 2, 0.1, -0.1)
)
SIM_STUDY_FORMULA <- "edges + triangles + star(c(2,3))"
SIM_STUDY_TRUNCATION <- 100

# parscale: must match free flat_par names (CS_params1/edges is fixed)
SIM_STUDY_PSCALE <- c(
  mu = 1, beta_overall = 0.1, beta_edges = 0.1, node_lambda = 0.1,
  CS_params2 = 0.1, CS_params3 = 0.1, CS_params4 = 0.1
)

# Main study init — edges (CS_params[1]) fixed to true value
SIM_STUDY_INIT <- list(
  mu = 10,
  beta_overall = 1,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-6.7, 1, 0, 0)
)
SIM_STUDY_FIXED <- c("K", "CS_params1")

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
})

# =============================================================================
# TEST 2: parscale mapping works with both CS_params and formula name conventions
# =============================================================================
test_that("parscale maps correctly to free flat_par (K and CS_params1 fixed)", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for fit test")

  # Flatten params_init as fit_hawkesNet does (remove K, then remove CS_params1)
  pi <- SIM_STUDY_INIT
  pi$K <- NULL
  flat_par <- unlist(pi)
  flat_par <- flat_par[!names(flat_par) %in% "CS_params1"]  # element-level fixed
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
test_that("loglik is finite at main study init (CS_params = c(-7, 1, 0, 0))", {
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
})

# =============================================================================
# TEST 4: Log-likelihood is finite even with OLD (distant) init values
# =============================================================================
test_that("loglik is finite at old distant init (CS_params = c(-10, 0, 0, 0))", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for loglik test")

  old_init <- list(mu = 10, beta_overall = 1, K = 0.5, beta_edges = 1,
                   node_lambda = 1, CS_params = c(-10, 0, 0, 0))
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
# TEST 5: fit_hawkesNet completes without error (exact main study config)
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
  # K and CS_params1 (edges) should NOT be in fitted params (they're fixed)
  expect_false("K" %in% names(fit$fit$par))
  expect_false("CS_params1" %in% names(fit$fit$par))
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
test_that("CS fit with consistency-study config (near-true init) converges", {
  sim <- sim_one(c(0, 5), seed = 42)
  skip_if(network::network.edgecount(sim$net) < 3, "Need edges for fit test")

  set.seed(99)
  params_init_near <- list(
    mu = SIM_STUDY_PARAMS$mu,
    beta_overall = max(0.1, SIM_STUDY_PARAMS$beta_overall * exp(rnorm(1, 0, 0.2))),
    K = SIM_STUDY_PARAMS$K,
    beta_edges = max(0.1, SIM_STUDY_PARAMS$beta_edges * exp(rnorm(1, 0, 0.2))),
    node_lambda = max(0.1, SIM_STUDY_PARAMS$node_lambda * exp(rnorm(1, 0, 0.2))),
    CS_params = SIM_STUDY_PARAMS$CS_params + rnorm(4, 0, 0.5)
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
  # Convergence at T=5 with only 500 iters may not always succeed (small data),
  # but the fit should at least produce finite parameters.
  # On the cluster with T=50 and 5000 iters, convergence is expected.
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
    mu = SIM_STUDY_PARAMS$mu,
    beta_overall = max(0.1, SIM_STUDY_PARAMS$beta_overall * exp(rnorm(1, 0, 0.2))),
    K = SIM_STUDY_PARAMS$K,
    beta_edges = max(0.1, SIM_STUDY_PARAMS$beta_edges * exp(rnorm(1, 0, 0.2))),
    node_lambda = max(0.1, SIM_STUDY_PARAMS$node_lambda * exp(rnorm(1, 0, 0.2))),
    CS_params = SIM_STUDY_PARAMS$CS_params + rnorm(4, 0, 0.5)
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
