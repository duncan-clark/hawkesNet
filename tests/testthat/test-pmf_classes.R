test_that("BA PMF", {
    data(net, package = "hawkesGrowthNet")
    time <- get_times(net)$times
    
    # Construct a BA kernel
    ker_ba <- BAKernel$new(
      params     = list(beta_edges = 0.1),
      filtration = net,
      opts       = list(max_node_time = 10)
    )
    
    # Get a drop-in legacy function
    PMF_mark_BA <- ker_ba$as_legacy_fun()
    
    pmf <- PMF_mark_BA(time[10],
                       params,
                       mark_filtration = net,
                       mark = NULL,
                       generate_mark = FALSE,
                       generate_density = TRUE,
                       grad = FALSE,
                       new_edge_hash = NULL,
                       truncation = NULL)
    expect_equal(pmf$mark_density,
                 0.3215021, 
                 tolerance = 0.01)

})