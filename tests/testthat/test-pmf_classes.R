test_that("BA PMF", {
    data(net, package = "hawkesGrowthNet")
    time <- get_times(net)$times
    params <-  list(beta_edges = 0.1)
    PMF_mark_BA <- function(time,
                            params,
                            mark_filtration,
                            mark = NULL,
                            generate_mark = FALSE,
                            generate_density = TRUE,
                            grad = FALSE,
                            new_edge_hash = NULL,
                            truncation = NULL, ...) {
      
      ker <- BAKernel$new(
        params     = params,
        filtration = mark_filtration,
        opts       = list(generate_mark = generate_mark,
                          generate_density = generate_density,
                          truncation = truncation,
                          max_node_time = 10,
                          new_edge_hash = new_edge_hash)
      )
      out <- ker$compute(time = time, mark = mark)
      
      # keep legacy field names
      list(
        mark_density            = out$mark_density,
        log_mark_density        = out$log_mark_density,
        edge_probs              = out$edge_probs,
        mark_sample             = out$mark_sample,
        mark_sample_density     = out$mark_sample_density,
        log_mark_sample_density = out$log_mark_sample_density
      )
    }
    
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