test_that("temporal_hawkes", {
    times <- c(0.623, 0.742, 0.788, 1.697, 1.724,
               2.053, 2.154, 2.171, 2.63, 2.731, 2.793,
               3.043, 3.277, 3.418, 3.497, 3.582, 3.74,
               3.797, 3.808, 4.117)
    fit_hgn <- fit_temporal_hawkes(params_init =  c(mu = 10, beta = 0.5, K = 0.1/0.5),
                                   realiz = data.frame(t = times, n = rep(1, length(times))),
                                   windowT = range(times),
                                   maxit = 1000)
    fit_hgn_nd <- fit_temporal_hawkes(params_init =  c(mu = 10, beta = 0.5, K = 0.1/0.5),
                               realiz = data.frame(t = times, n = rep(1, length(times))),
                               windowT = range(times),
                               maxit = 1000, density_approx = FALSE)
    expect_equal(fit_hgn$par[[1]],
                 5.723781,
                 tolerance = 0.1)
     expect_equal(fit_hgn_nd$par[[1]],
                 5.724134,
                 tolerance = 0.1)
}
)
