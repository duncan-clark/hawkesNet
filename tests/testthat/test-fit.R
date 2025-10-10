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
  
  # Construct a BA class based kernel
  PMF_mark_BA_class <- BAKernel$new(params = params_ba)$as_legacy_fun()
  fit_ba_class <- fit_hawkesGrowthNet(
    params_init = params_ba,
    time_window = c(0,max(get_times(mark_filtration)$times)),
    mark_filtration = mark_filtration,
    PMF_mark = PMF_mark_BA_class,
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
                     CS_params =  c(-10,0,0,0)
                     )
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
  PMF_mark_CS_class <- CSKernel$new(params = params_cs,
                                    opts = list(formula_RHS = "edges + triangles + star(c(2,3))",
                                                mark_decay = 'node_entrance',
                                                model = NULL,
                                                max_node_time = 1)
                                    )
  fit_cs_class <- fit_hawkesGrowthNet(params = params_cs,
                                      time_window = c(0,max(get_times(mark_filtration)$times)),
                                      mark_filtration =  mark_filtration,
                                      PMF_mark = PMF_mark_CS_class,
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
  
  # Test new class BA model
  expect_equal(fit_ba$fit$par,
               fit_ba_class$fit$par,
               tolerance = 0.01)
  # Test new class CS model
  expect_equal(fit_cs$fit$par,
               fit_cs_class$fit$par,
               tolerance = 0.01)
})