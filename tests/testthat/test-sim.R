test_that("simulation", {
    time_window <- c(0,10)
    ## BA
    set.seed(1234)
    params_ba <- list(mu = 1,
                      beta_overall = 20,
                      K = 0.1,
                      beta_edges = 0.1)
    sim_ba <- sim_hawkesGrowthNet(params =  params_ba,
                                  time_window = time_window,
                                  PMF_mark = PMF_mark,
                                  cond_intensity = cond_intensity,
                                  hashed_edges = TRUE,
                                  verbose = FALSE,
                                  mu_multiplier = 3,
                                  joint_accept = FALSE,
                                  truncation = 100)
    ## CS
    set.seed(4321)
    params_cs <-list(mu = 1,
                     beta_overall = 0.2,
                     K = 0.5,
                     beta_edges = 0.5,
                     node_lambda = 1,
                     CS_params = c(-6,0.5,0.3,-0.1))
    require(ernm)
    sim_cs <- sim_hawkesGrowthNet(params =  params_cs,
                                  time_window = time_window,
                                  PMF_mark = PMF_mark,
                                  cond_intensity = cond_intensity,
                                  hashed_edges = TRUE,
                                  verbose = FALSE,
                                  mu_multiplier = 1,
                                  joint_accept = FALSE,
                                  truncation = 10,
                                  type = "CS",
                                  formula_RHS = "edges  + triangles() + star(c(2,3))")

    ## expect
    expect_equal(sim_ba$events$t[1],
                 0.09495756,
                 tolerance = 0.01)
    expect_equal(sim_cs$events$t[1],
                 0.4384097,
                 tolerance = 0.01)
})
