# RHEM / repeated-hit mark PMF tests (experimental branch)

test_that("RHEM repeated-hit PMF favors repeated edges", {
    hits <- data.frame(
        t = c(0.1, 0.2, 0.3),
        i = c(1, 1, 2),
        j = c(2, 2, 1),
        weight = c(1, 1, 1)
    )
    params <- list(beta_edges = 0.1,
                   RHEM_params = c(intercept = 0,
                                   repetition = 2,
                                   reciprocity = 0))

    pmf <- PMF_mark_RHEM(time = 0.2,
                         params = params,
                         mark_filtration = hits[hits$t <= 0.2, ],
                         actors = 1:3,
                         grad = TRUE)

    expect_equal(length(pmf$edge_probs), 6)
    expect_equal(sum(pmf$edge_probs), 1)
    expect_gt(pmf$mark_density, 1 / 6)
    expect_named(pmf$mark_grad, names(params$RHEM_params))

    pmf_dispatch <- PMF_mark(time = 0.2,
                             params = params,
                             mark_filtration = hits[hits$t <= 0.2, ],
                             type = "RHEM",
                             actors = 1:3)
    expect_equal(pmf_dispatch$mark_density, pmf$mark_density)
})

test_that("RHEM repeated-hit PMF can use ERNM valued change stats", {
    hits <- data.frame(
        t = c(0.1, 0.2, 0.3),
        i = c(1, 1, 2),
        j = c(2, 2, 1),
        weight = c(1, 1, 1)
    )
    params <- list(beta_edges = 0.1,
                   RHEM_params = c(edgeValue = 1,
                                   recipValue = 0.25,
                                   senderValueActivity = -0.1,
                                   receiverValueActivity = 0.1))

    pmf <- PMF_mark_RHEM(time = 0.2,
                         params = params,
                         mark_filtration = hits[hits$t <= 0.2, ],
                         actors = 1:3,
                         formula_RHS = "edgeValue + recipValue + senderValueActivity + receiverValueActivity")

    expect_equal(length(pmf$edge_probs), 6)
    expect_equal(sum(pmf$edge_probs), 1)
    expect_true(is.finite(pmf$log_mark_density))
})

test_that("RHEM repeated-hit PMF can use ERNM valued triadic change stats", {
    hits <- data.frame(
        t = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
        i = c(1, 2, 1, 3, 1, 2, 3),
        j = c(2, 3, 3, 1, 4, 4, 4),
        weight = c(2, 3, 5, 7, 11, 13, 1)
    )
    params <- list(beta_edges = 0,
                   RHEM_params = c(transitiveValue = 0.01,
                                   cycleValue = 0.01,
                                   commonSourceValue = 0.01,
                                   commonTargetValue = 0.01))

    pmf <- PMF_mark_RHEM(time = 0.7,
                         params = params,
                         mark_filtration = hits,
                         actors = 1:4,
                         formula_RHS = paste("transitiveValue",
                                             "cycleValue",
                                             "commonSourceValue",
                                             "commonTargetValue",
                                             sep = " + "),
                         grad = TRUE)

    expect_equal(length(pmf$edge_probs), 12)
    expect_equal(sum(pmf$edge_probs), 1)
    expect_true(is.finite(pmf$log_mark_density))
    expect_named(pmf$mark_grad, names(params$RHEM_params))
})

test_that("RHEM formula change stats can be cached across repeated evaluations", {
    hits <- data.frame(
        t = c(0.1, 0.2, 0.3),
        i = c(1, 1, 2),
        j = c(2, 2, 1),
        weight = c(1, 1, 1)
    )
    params <- list(beta_edges = 0,
                   RHEM_params = c(edgeValue = 1,
                                   recipValue = 0.25))
    cache <- new.env(parent = emptyenv())

    pmf_first <- PMF_mark_RHEM(time = 0.2,
                               params = params,
                               mark_filtration = hits[hits$t <= 0.2, ],
                               actors = 1:3,
                               formula_RHS = "edgeValue + recipValue",
                               rhem_stats_cache = cache)
    cache_size <- length(ls(cache))
    pmf_second <- PMF_mark_RHEM(time = 0.2,
                                params = params,
                                mark_filtration = hits[hits$t <= 0.2, ],
                                actors = 1:3,
                                formula_RHS = "edgeValue + recipValue",
                                rhem_stats_cache = cache)

    expect_equal(cache_size, 1)
    expect_equal(length(ls(cache)), cache_size)
    expect_equal(pmf_second$edge_probs, pmf_first$edge_probs)
})

