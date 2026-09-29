test_that("CS marks cannot change attributes or birth times of existing vertices", {
  old <- network::network.initialize(2, directed = FALSE)
  network::set.vertex.attribute(old, "time", c(0, .1))
  network::set.vertex.attribute(old, "group", c("a", "b"))
  proposed <- network::network.copy(old)
  network::add.vertices(proposed, 1)
  network::set.vertex.attribute(proposed, "time", c(0, .1, 1))
  network::set.vertex.attribute(proposed, "group", c("a", "b", "a"))
  params <- list(m = 0, node_lambda = 1, beta_edges = 0, CS_params = 0,
                 vertex_categorical = list(group = c(a = .4)),
                 vertex_categorical_levels = list(group = c("a", "b")))
  for (mode in c("independent", "size_conditional")) {
    evaluate <- function(mark) PMF_mark_CS(1, params, old, mark = mark,
      formula_RHS = "edges", truncation = 3, cs_mode = mode)
    expect_gt(evaluate(proposed)$mark_density, 0)
    changed <- network::network.copy(proposed)
    network::set.vertex.attribute(changed, "group", c("b", "b", "a"))
    q <- evaluate(changed)
    expect_identical(q$log_mark_density, -Inf)
    expect_identical(q$log_density_func(params), -Inf)
    changed <- network::network.copy(proposed)
    network::set.vertex.attribute(changed, "time", c(.05, .1, 1))
    expect_identical(evaluate(changed)$log_mark_density, -Inf)
    expect_identical(network::get.vertex.attribute(old, "group"), c("a", "b"))
  }
})
