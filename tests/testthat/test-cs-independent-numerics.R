ci_numeric_graph <- function(n = 3L) {
  g <- network::network.initialize(n, directed = FALSE)
  network::set.vertex.attribute(g, "time", rep(0, n))
  g
}

ci_numeric_pmf <- function(g, params, mark = NULL, generate = FALSE,
                           nonempty = TRUE, window = 4L, rhs = "edges") {
  .pmf_cs_independent(1, params, g, mark, generate, rhs, window,
                     "node_entrance", FALSE, NULL, nonempty)
}

test_that("CS-1 tiny attempt means retain nonempty support without rejection", {
  g <- ci_numeric_graph()
  params <- list(m = 1e-300, beta_edges = .4, node_lambda = 0, CS_params = -800)
  set.seed(182)
  draw <- ci_numeric_pmf(g, params, generate = TRUE)
  expect_equal(network::network.edgecount(draw$mark_sample), 1)
  expect_equal(draw$mark_density, 1/3, tolerance = 1e-12)
  expect_equal(draw$edge_probs, rep(1/3, 3), tolerance = 1e-12)
  expect_equal(draw$density_func(params), draw$mark_density)
  expect_identical(parent.env(environment(draw$log_density_func)), baseenv())
  expect_false(any(c("old", "mark", "model", "zero_local") %in%
                     ls(environment(draw$log_density_func))))
  expect_false("states" %in% names(environment(draw$log_density_func)$table_local))
})

test_that("CS-1 normalized log weights survive extreme common scores", {
  params <- list(m = 1, beta_edges = 0, node_lambda = 0, CS_params = -1e308)
  tab <- list(tails = 1:3, changes = matrix(1, 3, 1), ages = rep(1, 3))
  dist <- .cs_independent_distribution(tab, params)
  expect_equal(dist$selection_probs, rep(1/3, 3), tolerance = 1e-14)
  expect_equal(dist$probs, rep(-expm1(-1/3), 3), tolerance = 1e-14)
  # Keep the log mass of an inclusion whose ordinary probability underflows.
  tab$changes <- matrix(c(0, 1), 2, 1)
  tab$tails <- 1:2; tab$ages <- c(0, 0); params$CS_params <- -1000
  dist <- .cs_independent_distribution(tab, params)
  expect_true(is.finite(.cs_independent_log_prob(dist, 2L)))
  expect_lt(.cs_independent_log_prob(dist, 2L), -900)
})

test_that("CS-1 supports more than twelve candidates without enumeration", {
  g <- ci_numeric_graph(7)
  params <- list(m = 1.4, beta_edges = .4, node_lambda = 0, CS_params = 0)
  set.seed(183)
  draw <- ci_numeric_pmf(g, params, generate = TRUE, window = 7L)
  expect_equal(draw$n_candidate_edges, 21L)
  expect_equal(draw$n_mark_states, 2^21)
  expect_equal(draw$reference_edge_probs, rep(-expm1(-1.4/21), 21))
  expect_equal(nrow(draw$edge_change_stats), 21L)
  expect_true(is.finite(draw$log_mark_density))
})

test_that("CS-1 handles zero attempts, birth-only updates and empty supports", {
  g <- ci_numeric_graph()
  params <- list(m = 0, beta_edges = .4, node_lambda = 0, CS_params = 0)
  expect_error(ci_numeric_pmf(g, params, generate = TRUE), "No nonempty")
  expect_equal(ci_numeric_pmf(g, params, mark = g, nonempty = FALSE)$mark_density, 1)
  params$node_lambda <- 1e-300
  draw <- ci_numeric_pmf(g, params, generate = TRUE)
  expect_equal(network::network.size(draw$mark_sample), 4)
  expect_equal(network::network.edgecount(draw$mark_sample), 0)
  expect_equal(draw$mark_density, 1, tolerance = 1e-12)
  params$node_lambda <- 5e-324
  draw <- ci_numeric_pmf(g, params, generate = TRUE)
  expect_equal(network::network.size(draw$mark_sample), 4)
  expect_equal(draw$mark_density, 1, tolerance = 1e-12)
  full <- network::network(matrix(1, 3, 3) - diag(3), directed = FALSE)
  network::set.vertex.attribute(full, "time", rep(0, 3))
  network::set.edge.attribute(full, "time", rep(0, 3))
  params$m <- 1; params$node_lambda <- 0
  expect_error(ci_numeric_pmf(full, params, generate = TRUE), "No nonempty")
  expect_equal(ci_numeric_pmf(full, params, mark = full, nonempty = FALSE)$mark_density, 1)
  params$node_lambda <- 1
  seed <- ci_numeric_graph(0)
  expect_true(is.finite(ci_numeric_pmf(seed, params, generate = TRUE)$log_mark_density))
})

test_that("CS-1 freezes single-edge statistics and reports the whole change", {
  g <- ci_numeric_graph()
  network::add.edge(g, 1, 2)
  network::set.edge.attribute(g, "time", 0)
  network::set.vertex.attribute(g, "unused_metadata", letters[1:3])
  original <- serialize(g, NULL)
  mark <- network::network.copy(g)
  network::add.edges(mark, c(1, 2), c(3, 3))
  network::set.edge.attribute(mark, "time", c(0, 1, 1))
  params <- list(m = 1, beta_edges = 0, node_lambda = 0, CS_params = c(0, 5))
  out <- ci_numeric_pmf(g, params, mark = mark, rhs = "edges + triangles")
  expect_equal(unname(out$edge_change_stats[, "triangles"]), c(0, 0))
  expect_equal(unname(out$mark_change_stats["triangles"]), 1)
  expect_equal(out$edge_selection_probs, c(.5, .5))
  expect_identical(serialize(g, NULL), original)
  expect_equal(network::get.vertex.attribute(mark, "unused_metadata"), letters[1:3])
})

test_that("CS-1 rejects observed multiedges instead of silently collapsing them", {
  g <- ci_numeric_graph()
  multiedges <- network::network.initialize(3, directed = FALSE, multiple = TRUE)
  network::set.vertex.attribute(multiedges, "time", rep(0, 3))
  network::add.edges(multiedges, c(1, 1), c(2, 2))
  network::set.edge.attribute(multiedges, "time", c(1, 1))
  params <- list(m = 1, beta_edges = 0, node_lambda = 0, CS_params = 0)
  expect_equal(ci_numeric_pmf(g, params, mark = multiedges)$log_mark_density, -Inf)
})
