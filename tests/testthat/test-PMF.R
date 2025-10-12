test_that("BA PMF", {
    data(net, package = "hawkesGrowthNet")
    time <- get_times(net)$times
    params <-  list(beta_edges = 0.1)
    pmf <- PMF_mark_BA(time[10],  params, mark_filtration = net, mark = NULL,
                       generate_mark = FALSE,  generate_density = TRUE,
                       grad = FALSE, new_edge_hash = NULL, truncation = NULL)
    PMF_mark_BA_class <- BAKernel$new(params = params)$as_legacy_fun()
    pmf_class <- PMF_mark_BA_class(time[10],  params, mark_filtration = net, mark = NULL,
                       generate_mark = FALSE,  generate_density = TRUE,
                       grad = FALSE, new_edge_hash = NULL, truncation = NULL)
    
    expect_equal(pmf$mark_density,
                 0.3215021,
                 tolerance = 0.01)
    
    expect_equal(pmf$mark_density,
                 pmf_class$mark_density,
                 tolerance = 0.01)
})
test_that("BA-bip PMF", {
    events_list <- list(
        i = c("A", "B", "C", "D","D"),
        j = c(2, 3, 1, 2, 1),
        t = c(5, 10, 15, 20, 21)
    )
    net <- events_to_bipartite_net(events_list)
    time <- get_times(net)$times
    params <-  list(beta_edges = 0.1, lambda_new = 1.5)
    pmf <- PMF_mark_BA_bipartite(time[3],
                                 params, mark_filtration = net, mark = NULL,
                       generate_mark = FALSE,  generate_density = TRUE,
                       grad = FALSE, new_edge_hash = NULL, truncation = NULL)
    expect_equal(pmf$mark_density,
                 1,
                 tolerance = 0.01)
})
test_that("CS PMF", {
    require(ernm)
    data(net, package = "hawkesGrowthNet")
    time <- get_times(net)$times
    params <-  list(beta_edges = 0.1,
                    node_lambda = 1,
                    CS_params =  c(-10,0,0,0))
    mark_filtration <-  filtration_to_net(net,10)
    pmf <- PMF_mark_CS(time = time[10],  params = params,
                       mark_filtration = mark_filtration,
                       mark = NULL,
                       generate_mark = FALSE,  generate_density = TRUE,
                       grad = FALSE, new_edge_hash = NULL, truncation = 1,
                       formula_RHS = "edges + triangles + star(c(2,3))",
                       mark_decay = 'node_entrance',
                       model = NULL, max_node_time = 1)
    PMF_mark_CS_class <- CSKernel$new(params = params,
                                      opts = list(formula_RHS = "edges + triangles + star(c(2,3))",
                                                  mark_decay = 'node_entrance',
                                                  model = NULL,
                                                  max_node_time = 1)
                                      )$as_legacy_fun()
    pmf_class <- PMF_mark_CS_class(time[10],
                                   params,
                                   mark_filtration = net,
                                   mark = NULL,
                                   generate_mark = FALSE,
                                   generate_density = TRUE,
                                   grad = FALSE,
                                   new_edge_hash = NULL,
                                   truncation = NULL)
    
    expect_equal(pmf$mark_density,
                 0.3678293,
                 tolerance = 0.01)
    expect_equal(pmf$mark_density,
                 pmf_class$mark_density,
                 tolerance = 0.01)
    
})
test_that("PMF", {
    data(net, package = "hawkesGrowthNet")
    time <- get_times(net)$times
    ## BA
    params_ba <-  list(beta_edges = 0.1)
    mark_filtration <-  filtration_to_net(net,10)
    pmf_ba <- PMF_mark(time[10],  params_ba, mark_filtration)
    ## CS
    require(ernm)
    params_cs <-  list(beta_edges = 0.1,
                       node_lambda = 1,
                       CS_params =  c(-10,0,0,0))
    pmf_cs <- PMF_mark(time = time[10],  params = params_cs,
                       mark_filtration = mark_filtration,
                       type = "CS",  truncation = 1,
                       formula_RHS = "edges + triangles + star(c(2,3))",
                       max_node_time = 1)
    ## expect
    expect_equal(pmf_ba$mark_density,
                 0.3215021,
                 tolerance = 0.01)
    expect_equal(pmf_cs$mark_density,
                 0.3678293,
                 tolerance = 0.01)
    expect_error(PMF_mark(time[10],  params_ba, mark_filtration, type = "CA"))
})
