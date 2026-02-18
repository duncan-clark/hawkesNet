test_that("expand_vertex_categorical_probs works", {
  p <- expand_vertex_categorical_probs(c(male = 0.4, female = 0.4),
                                        c("male", "female", "unknown"))
  expect_length(p, 3)
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_equal(names(p), c("male", "female", "unknown"))
  expect_equal(unname(p[3]), 0.2)
})

test_that("expand_vertex_categorical_probs returns NULL on bad input", {
  expect_null(expand_vertex_categorical_probs(NULL, c("a", "b")))
  expect_null(expand_vertex_categorical_probs(c(0.5), NULL))
  expect_null(expand_vertex_categorical_probs(c(0.5, 0.3), c("a", "b")))
})

test_that("get_truncated_candidates returns correct structure", {
  el <- matrix(c(1, 2, 2, 3, 3, 4, 4, 5), ncol = 2, byrow = TRUE)
  net <- network::network(el, directed = FALSE)

  cands <- get_truncated_candidates(net, 6, 5, 3, "node_entrance")
  expect_type(cands, "list")
  expect_true("tails" %in% names(cands))
  expect_true("heads" %in% names(cands))
  expect_true(all(cands$tails > 0))
  expect_true(all(cands$heads > 0))
})

test_that("get_truncated_candidates with growth_only restricts edges", {
  el <- matrix(c(1, 2, 2, 3, 3, 4, 4, 5), ncol = 2, byrow = TRUE)
  net <- network::network(el, directed = FALSE)

  cands_go <- get_truncated_candidates(net, 6, 5, 3, "node_entrance", growth_only = TRUE)
  cands_all <- get_truncated_candidates(net, 6, 5, 3, "node_entrance", growth_only = FALSE)

  # growth_only should have a subset of candidates
  expect_lte(length(cands_go$tails), length(cands_all$tails))
})

test_that("PMF_mark_BA runs and returns log_mark_density", {
  params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
  el <- matrix(c(1, 2, 2, 3, 3, 4), ncol = 2, byrow = TRUE)
  net <- network::network(el, directed = FALSE)
  network::set.vertex.attribute(net, "time", c(0.1, 0.2, 0.3, 0.4))

  pmf <- PMF_mark_BA(0.5, params, net)
  expect_type(pmf, "list")
  expect_true("log_mark_density" %in% names(pmf))
  expect_true(is.finite(pmf$log_mark_density))
})
