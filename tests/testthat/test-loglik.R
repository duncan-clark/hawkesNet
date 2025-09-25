test_that("conditinoal_intensity", {
    data(net, package = "hawkesGrowthNet")
    ## BA
    params_ba <-  list(mu = 10,
                       beta_overall = 0.1,
                       K = 0.1,
                       beta_edges = 0.1,
                       node_lambda = 1)
    time <- get_times(net)$times
    intensity_ba <- cond_intensity(new_net = NULL, t = time[10], 
                                        mark_filtration = net,
                                        PMF_mark = PMF_mark,  
                                        params = params_ba)
    ## CS
    require(ernm)
    params_cs <-  list(mu = 10,
                       beta_overall = 0.1,
                       K = 0.1,
                       beta_edges = 0.1,
                       node_lambda = 1,
                       CS_params =  c(-10,0,0,0))
    intensity_cs <- cond_intensity(new_net = NULL, t = time[10], 
                                mark_filtration = net,
                                PMF_mark = PMF_mark,  
                                params = params_cs,
                                type = "CS",
                                truncation = 1,
                                formula_RHS = "edges + triangles + star(c(2,3))",
                                max_node_time = 1 )
    ## expect
    expect_equal(intensity_ba$result,
                 3.504372,, 
                 tolerance = 0.01)
    expect_equal(intensity_cs$result,
                 4.007338,
                 tolerance = 0.01)
})
test_that("loglik", {
    data(net, package = "hawkesGrowthNet")
    require(network)
    ## BA
    params_ba <-  list(mu = 10,
                       beta_overall = 0.1,
                       K = 0.1,
                       beta_edges = 0.1,
                       node_lambda = 1)
    mark_filtration <-  network::get.inducedSubgraph(net, v = 1:10)
    loglik_ba <- loglik_hawkesGrowthNet(params = params_ba,
                                        time_window = c(0,max(get_times(mark_filtration)$times)),
                                        mark_filtration =  mark_filtration,
                                        PMF_mark = PMF_mark,  
                                        verbose = FALSE)
    ## CS
    require(ernm)
    params_cs <-  list(mu = 10,
                       beta_overall = 0.1,
                       K = 0.1,
                       beta_edges = 0.1,
                       node_lambda = 1,
                       CS_params =  c(-10,0,0,0))
    loglik_cs <- loglik_hawkesGrowthNet(params = params_cs,
                                        time_window = c(0,max(get_times(mark_filtration)$times)),
                                        mark_filtration =  mark_filtration,
                                        PMF_mark = PMF_mark, type = "CS",
                                        truncation = 1,
                                        formula_RHS = "edges + triangles + star(c(2,3))",
                                        max_node_time = 1,
                                        verbose = FALSE)
    ## expect
    expect_equal(loglik_ba$loglik,
                 6.243956, 
                 tolerance = 0.01)
    expect_equal(loglik_cs$loglik,
                 18.67455,
                 tolerance = 0.01)
})
test_that("fit", {
    skip_on_cran() ## takes too long
    data(net, package = "hawkesGrowthNet")
    mark_filtration <-  network::get.inducedSubgraph(net, v = 1:50)
    time <- get_times(mark_filtration)$times
    params_ba <-  list(mu = 10,
                    beta_overall = 0.1,
                    K = 0.1,
                    beta_edges = 0.1)
    fit_ba <- fit_hawkesGrowthNet(
        params_init = params_ba,
        time_window = c(0,max(get_times(mark_filtration)$times)),
        mark_filtration = mark_filtration,
        PMF_mark = PMF_mark,
        grad = FALSE,
        trace = 0,
        maxit = 100,
        truncation = 10,
        get_hessian = FALSE)

    ## CS
    require(ernm)
    params_cs <-  list(mu = 10,
                       beta_overall = 0.1,
                       K = 0.1,
                       beta_edges = 0.1,
                       node_lambda = 1,
                       CS_params =  c(-10,0,0,0))
    fit_cs <- fit_hawkesGrowthNet(params = params_cs,
                                  time_window = c(0,max(get_times(mark_filtration)$times)),
                                  mark_filtration =  mark_filtration,
                                  PMF_mark = PMF_mark, type = "CS",
                                  truncation = 10,
                                  formula_RHS = "edges + triangles + star(c(2,3))",
                                  max_node_time = 10,
                                  grad = FALSE,
                                  trace = 0,
                                  maxit = 100,
                                  get_hessian = FALSE)
    ## expect
    expect_equal(fit_ba$fit$par[[1]],
                 78.48635,
                 tolerance = 0.1)
    expect_equal(fit_cs$fit$par[[1]],
                 9.69961436,
                 tolerance = 0.01)
})
