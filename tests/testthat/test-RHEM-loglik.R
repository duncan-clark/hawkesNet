# RHEM loglik/score tests

test_that("RHEM repeated-hit loglik and fit run", {
    hits <- data.frame(
        t = c(0.1, 0.2, 0.35, 0.55, 0.8),
        i = c(1, 1, 2, 1, 3),
        j = c(2, 2, 1, 3, 1),
        weight = c(1, 1, 1, 1, 1)
    )
    params <- list(mu = 2,
                   beta_overall = 1,
                   K = 0.2,
                   beta_edges = 0.1,
                   RHEM_params = c(edgeValue = 0.5,
                                   recipValue = 0.25,
                                   senderValueActivity = -0.1,
                                   receiverValueActivity = 0.1))

    loglik <- loglik_hawkesNet(params = params,
                                     time_window = c(0, 1),
                                     mark_filtration = hits,
                                     PMF_mark = PMF_mark,
                                     type = "RHEM",
                                     actors = 1:3,
                                     formula_RHS = "edgeValue + recipValue + senderValueActivity + receiverValueActivity",
                                     verbose = FALSE)
    expect_true(is.finite(loglik$loglik))

    fit <- fit_hawkesNet(params_init = params,
                               time_window = c(0, 1),
                               mark_filtration = hits,
                               PMF_mark = PMF_mark,
                               type = "RHEM",
                               actors = 1:3,
                               formula_RHS = "edgeValue + recipValue + senderValueActivity + receiverValueActivity",
                               fixed_params = c("mu", "beta_overall", "K", "beta_edges"),
                               trace = 0,
                               maxit = 2,
                               get_hessian = FALSE)
    expect_true(all(is.finite(fit$fit$par)))
})

