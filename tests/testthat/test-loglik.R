test_that("loglik", {
    data(net, package = "hawkesGrowthNet")
    time <- get_times(net)$times
    require(network)
    ## BA
    params_ba <-  list(mu = length(time)/max(time),
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
    params_cs <-  list(mu = length(time)/max(time),
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
                 -12725,
                 tolerance = 1)
    expect_equal(loglik_cs$loglik,
                 -12712.57,
                 tolerance = 1)
})
