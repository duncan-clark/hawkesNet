test_that("conditional_intensity", {
    data(net, package = "hawkesGrowthNet")
    # ===========
    ## BA
    # ===========
    params_ba <-  list(mu = 10,
                       beta_overall = 0.1,
                       K = 0.1,
                       beta_edges = 0.1)
    time <- get_times(net)$times
    intensity_ba <- cond_intensity(new_net = NULL, t = time[10], 
                                        mark_filtration = net,
                                        PMF_mark = PMF_mark,  
                                        params = params_ba)
    
    # Construct a BA class based kernel
    PMF_mark_BA_class <- BAKernel$new(params = params_ba)$as_legacy_fun()
    
    intensity_ba_class <- cond_intensity(new_net = NULL, t = time[10], 
                                         mark_filtration = net,
                                         PMF_mark = PMF_mark_BA_class,
                                         params = params_ba)
    
    # ===========
    ## CS
    # ===========
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
    
    PMF_mark_CS_class <- CSKernel$new(params = params_cs,
                                      opts = list(formula_RHS = "edges + triangles + star(c(2,3))",
                                                  mark_decay = 'node_entrance',
                                                  model = NULL,
                                                  max_node_time = 1)
    )$as_legacy_fun()
    
    intensity_cs_class <- cond_intensity(new_net = NULL, t = time[10], 
                                   mark_filtration = net,
                                   PMF_mark = PMF_mark_CS_class,
                                   params = params_cs
                                   )
    
    ## expect
    expect_equal(intensity_ba$result,
                 3.504372, 
                 tolerance = 0.01)
    expect_equal(intensity_cs$result,
                 4.007338,
                 tolerance = 0.01)
    
    # Test new class BA model
    expect_equal(intensity_ba$result,
                 intensity_ba_class$result,
                 tolerance = 0.01)
    # Test new class BA model
    expect_equal(intensity_cs$result,
                 intensity_cs_class$result,
                 tolerance = 0.01)
})
test_that("loglik", {
    data(net, package = "hawkesGrowthNet")
    require(network)
    ## BA
    params_ba <-  list(mu = 10,
                       beta_overall = 0.1,
                       K = 0.1,
                       beta_edges = 0.1)
    mark_filtration <-  network::get.inducedSubgraph(net, v = 1:10)
    loglik_ba <- loglik_hawkesGrowthNet(params = params_ba,
                                        time_window = c(0,max(get_times(mark_filtration)$times)),
                                        mark_filtration =  mark_filtration,
                                        PMF_mark = PMF_mark,  
                                        verbose = FALSE)
    # Construct a BA class based kernel
    PMF_mark_BA_class <- BAKernel$new(params = params_ba)$as_legacy_fun()
    loglik_ba_class <- loglik_hawkesGrowthNet(params = params_ba,
                                              time_window = c(0,max(get_times(mark_filtration)$times)),
                                              mark_filtration =  mark_filtration,
                                              PMF_mark = PMF_mark_BA_class,  
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
    
    PMF_mark_CS_class <- CSKernel$new(params = params_cs,
                                      opts = list(formula_RHS = "edges + triangles + star(c(2,3))",
                                                  mark_decay = 'node_entrance',
                                                  model = NULL,
                                                  max_node_time = 1)
    )$as_legacy_fun()
    loglik_cs_class <- loglik_hawkesGrowthNet(params = params_cs,
                                              time_window = c(0,max(get_times(mark_filtration)$times)),
                                              mark_filtration =  mark_filtration,
                                              PMF_mark = PMF_mark_CS_class,
                                              verbose = FALSE)
    
    ## ==========================================================
    ## Comment on breaking changes
    ## - Filtration -> network was not deleting edges properly
    ## - This means the conditionaly intensities for early in the processs events was wrong
    ## - Didn't seem to effect late process points
    ## - Much more confident this is correcet
    ## - Also tightend the equal or not equal version so equal == FALSE gives you the strictly less than version
    ## =========================================================
    
    ## expect
    expect_equal(loglik_ba$loglik,
                 7.680491, 
                 tolerance = 0.01)
    expect_equal(loglik_cs$loglik,
                 8.289209,
                 tolerance = 0.01)
    
    # Test new class BA model
    expect_equal(-21.22668,
                 loglik_ba_class$loglik,
                 tolerance = 0.01)
    
    expect_equal(-53.41084,
                 loglik_cs_class$loglik,
                 tolerance = 0.01)
})

