# Synthetic valued RHEM / timeNet fits must leave the uniform (all-zero) solution.

.rhem_repetition_history <- function(n_events = 24L, n_actors = 4L, seed = 11L) {
    set.seed(seed)
    actors <- seq_len(n_actors)
    # Strong repetition signal: after a burn-in edge, keep replaying the same dyad.
    t <- seq(0.05, 0.95, length.out = n_events)
    i <- integer(n_events)
    j <- integer(n_events)
    i[1] <- 1L
    j[1] <- 2L
    for (k in 2:n_events) {
        if (runif(1) < 0.85) {
            i[k] <- i[k - 1L]
            j[k] <- j[k - 1L]
        } else {
            dyad <- sample(actors, 2L)
            i[k] <- dyad[1L]
            j[k] <- dyad[2L]
        }
    }
    data.frame(t = t, i = i, j = j, weight = rep(1, n_events))
}

test_that("point_process_params_valid allows no-excitation K=0 and beta_edges=0", {
    expect_true(point_process_params_valid(list(
        mu = 1, beta_overall = 1, K = 0, beta_edges = 0
    )))
    expect_false(point_process_params_valid(list(
        mu = 1, beta_overall = 1, K = -0.1, beta_edges = 0
    )))
    expect_false(point_process_params_valid(list(
        mu = 1, beta_overall = 1, K = 0, beta_edges = -0.1
    )))
})

test_that("RHEM PMF exposes log_density_func that tracks RHEM_params", {
    hits <- data.frame(
        t = c(0.1, 0.2, 0.3),
        i = c(1, 1, 2),
        j = c(2, 2, 1),
        weight = c(1, 1, 1)
    )
    params0 <- list(beta_edges = 0,
                    RHEM_params = c(repetition = 0, reciprocity = 0))
    params1 <- list(beta_edges = 0,
                    RHEM_params = c(repetition = 2, reciprocity = 0))

    pmf <- PMF_mark_RHEM(time = 0.2,
                         params = params0,
                         mark_filtration = hits[hits$t <= 0.2, ],
                         actors = 1:3)
    expect_true(is.function(pmf$log_density_func))
    expect_equal(pmf$log_density_func(params0), pmf$log_mark_density, tolerance = 1e-10)

    ll0 <- pmf$log_density_func(params0)
    ll1 <- pmf$log_density_func(params1)
    expect_gt(ll1, ll0)
    expect_gt(ll1, -log(6))
})

test_that("cached intensity closures update when valued RHEM_params change", {
    skip_if_not_installed("timeNet")
    # Deterministic near-pure repetition on one dyad.
    hits <- data.frame(
        t = seq(0.05, 0.95, length.out = 12L),
        i = c(1L, rep(1L, 11L)),
        j = c(2L, rep(2L, 11L)),
        weight = 1
    )
    params0 <- list(mu = 5,
                    beta_overall = 1,
                    K = 0,
                    beta_edges = 0,
                    RHEM_params = c(edgeValue = 0, recipValue = 0))
    params1 <- params0
    params1$RHEM_params <- c(edgeValue = 2.5, recipValue = 0)

    ll <- suppressMessages(loglik_hawkesNet(
        params = params0,
        time_window = c(0, 1),
        mark_filtration = hits,
        PMF_mark = PMF_mark,
        type = "RHEM",
        actors = 1:4,
        formula_RHS = "edgeValue + recipValue",
        verbose = FALSE
    ))
    expect_true(length(ll$intens_funcs) >= 2)
    dens0 <- vapply(ll$intens_funcs, function(f) f(params0), numeric(1))
    dens1 <- vapply(ll$intens_funcs, function(f) f(params1), numeric(1))
    expect_true(all(is.finite(dens0)))
    expect_true(all(is.finite(dens1)))
    # Later events (after history) must improve under a positive edgeValue.
    later <- 4:length(dens0)
    expect_gt(sum(log(dens1[later])), sum(log(dens0[later])))
    # Closures must not be frozen at the uniform init value.
    expect_false(isTRUE(all.equal(dens0, dens1)))
})

test_that("valued RHEM fit moves coefficients away from zero on repetition signal", {
    skip_if_not_installed("timeNet")
    hits <- .rhem_repetition_history(n_events = 36L, n_actors = 4L, seed = 7L)
    n_dyads <- 4L * 3L
    params_init <- list(mu = nrow(hits),
                        beta_overall = 1,
                        K = 0,
                        beta_edges = 0,
                        RHEM_params = c(edgeValue = 0, recipValue = 0))

    fit <- suppressMessages(fit_hawkesNet(
        params_init = params_init,
        time_window = c(0, 1),
        mark_filtration = hits,
        PMF_mark = PMF_mark,
        type = "RHEM",
        actors = 1:4,
        formula_RHS = "edgeValue + recipValue",
        fixed_params = c("mu", "beta_overall", "K", "beta_edges"),
        cache_intensity = TRUE,
        verbose = FALSE,
        maxit = 120,
        get_hessian = FALSE,
        trace = 0
    ))

    est <- fit$fit$par
    expect_true(all(is.finite(est)))
    expect_gt(abs(est[["RHEM_params.edgeValue"]]), 0.15)

    params_hat <- params_init
    params_hat$RHEM_params <- c(edgeValue = unname(est[["RHEM_params.edgeValue"]]),
                                recipValue = unname(est[["RHEM_params.recipValue"]]))

    # Fitted mark model must beat the uniform baseline on later events.
    score_idx <- 8:nrow(hits)
    nll_hat <- mean(vapply(score_idx, function(idx) {
        -PMF_mark_RHEM(time = hits$t[idx],
                       params = params_hat,
                       mark_filtration = hits[seq_len(idx), ],
                       actors = 1:4,
                       formula_RHS = "edgeValue + recipValue")$log_mark_density
    }, numeric(1)))
    expect_lt(nll_hat, log(n_dyads) - 0.05)

    pmf <- PMF_mark_RHEM(time = hits$t[nrow(hits)],
                         params = params_hat,
                         mark_filtration = hits,
                         actors = 1:4,
                         formula_RHS = "edgeValue + recipValue")
    expect_gt(max(pmf$edge_probs), 1 / n_dyads + 0.05)
    expect_false(isTRUE(all.equal(pmf$edge_probs, rep(1 / n_dyads, n_dyads))))
})

test_that("merge_fit_params reconstructs nested RHEM_params from flat optim vector", {
    params_init <- list(mu = 10,
                        beta_overall = 1,
                        K = 0,
                        beta_edges = 0,
                        RHEM_params = c(repetition = 0, reciprocity = 0))
    par_vec <- c(mu = 12,
                 RHEM_params.repetition = 1.5,
                 RHEM_params.reciprocity = 0.4)
    merged <- merge_fit_params(par_vec, params_init)
    expect_equal(merged$mu, 12)
    expect_equal(unname(merged$RHEM_params[["repetition"]]), 1.5)
    expect_equal(unname(merged$RHEM_params[["reciprocity"]]), 0.4)
    expect_equal(merged$K, 0)
})

test_that("free-mu RHEM fit with auto parscale escapes null marks", {
    hits <- .rhem_repetition_history(n_events = 28L, n_actors = 4L, seed = 5L)
    params_init <- list(mu = 2000,  # Ethereum-scale rate vs O(1) marks
                        beta_overall = 1,
                        K = 0,
                        beta_edges = 0,
                        RHEM_params = c(repetition = 0.01,
                                        reciprocity = 0.01,
                                        sender_activity = 0.01,
                                        receiver_activity = 0.01))
    fit <- suppressMessages(fit_hawkesNet(
        params_init = params_init,
        time_window = c(0, 1),
        mark_filtration = hits,
        PMF_mark = PMF_mark,
        type = "RHEM",
        actors = 1:4,
        formula_RHS = NULL,
        fixed_params = c("beta_overall", "K", "beta_edges"),
        cache_intensity = TRUE,
        verbose = FALSE,
        method = "BFGS",
        maxit = 80,
        get_hessian = FALSE,
        trace = 0
    ))
    expect_gt(abs(fit$fit$par[["RHEM_params.repetition"]]), 0.2)
    expect_gt(abs(fit$params$RHEM_params[["repetition"]]), 0.2)
})

test_that("valued timeNet + Hawkes (K free) recovers nonzero mark coeffs", {
    skip_if_not_installed("timeNet")
    hits <- .rhem_repetition_history(n_events = 40L, n_actors = 4L, seed = 9L)
    hits$weight <- runif(nrow(hits), 0.5, 2)
    formula_RHS <- "decayedEdgeValue(beta = 2) + decayedRecipValue(beta = 2)"
    params_init <- list(mu = nrow(hits),
                        beta_overall = 1,
                        K = 0.1,
                        beta_edges = 0,
                        RHEM_params = c(decayedEdgeValue = 0.01,
                                        decayedRecipValue = 0.01))
    fit <- suppressMessages(fit_hawkesNet(
        params_init = params_init,
        time_window = c(0, 1),
        mark_filtration = hits,
        PMF_mark = PMF_mark,
        type = "RHEM",
        actors = 1:4,
        formula_RHS = formula_RHS,
        fixed_params = c("beta_edges"),
        cache_intensity = TRUE,
        verbose = FALSE,
        method = "BFGS",
        maxit = 100,
        get_hessian = FALSE,
        trace = 0
    ))
    expect_true(all(is.finite(fit$fit$par)))
    expect_gt(abs(fit$params$RHEM_params[["decayedEdgeValue"]]), 0.1)
})

test_that("no-excitation HawkesNet RHEM fit with K=0 escapes uniform collapse", {
    hits <- .rhem_repetition_history(n_events = 20L, n_actors = 3L, seed = 3L)
    n_dyads <- 3L * 2L
    params_init <- list(mu = nrow(hits),
                        beta_overall = 1,
                        K = 0,
                        beta_edges = 0,
                        RHEM_params = c(repetition = 0,
                                        reciprocity = 0,
                                        sender_activity = 0,
                                        receiver_activity = 0))

    fit <- suppressMessages(fit_hawkesNet(
        params_init = params_init,
        time_window = c(0, 1),
        mark_filtration = hits,
        PMF_mark = PMF_mark,
        type = "RHEM",
        actors = 1:3,
        formula_RHS = NULL,
        fixed_params = c("mu", "beta_overall", "K", "beta_edges"),
        cache_intensity = TRUE,
        verbose = FALSE,
        maxit = 60,
        get_hessian = FALSE,
        trace = 0
    ))

    est <- fit$fit$par
    expect_true(all(is.finite(est)))
    expect_gt(abs(est[["RHEM_params.repetition"]]), 0.2)

    params_hat <- params_init
    params_hat$RHEM_params <- c(
        repetition = unname(est[["RHEM_params.repetition"]]),
        reciprocity = unname(est[["RHEM_params.reciprocity"]]),
        sender_activity = unname(est[["RHEM_params.sender_activity"]]),
        receiver_activity = unname(est[["RHEM_params.receiver_activity"]])
    )
    pmf <- PMF_mark_RHEM(time = hits$t[nrow(hits)],
                         params = params_hat,
                         mark_filtration = hits,
                         actors = 1:3)
    expect_gt(max(pmf$edge_probs), 1 / n_dyads + 0.05)
})
