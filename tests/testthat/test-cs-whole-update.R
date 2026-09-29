cs_test_graph <- function() {
  g <- network::network.initialize(4, directed = FALSE)
  network::add.edge(g, 1, 2)
  network::set.vertex.attribute(g, "time", rep(0, 4))
  network::set.edge.attribute(g, "time", 0)
  g
}

cs_test_params <- function() list(mu = 2, K = .2, beta_overall = 1,
  m = 1.3, beta_edges = .4, node_lambda = .7, CS_params = c(0, .8, .05, -.02))

test_that("sequential calculations recover the unordered whole-update change", {
  g <- cs_test_graph()
  edges <- rbind(c(1, 3), c(2, 3))
  rhs <- "edges + triangles + star(c(2,3))"
  a <- multi_edge_change_stats(g, edges, rhs)
  b <- multi_edge_change_stats(g, edges[2:1, ], rhs)
  expect_equal(unname(a$increments[, "triangles"]), c(0, 1))
  expect_equal(a$total, b$total)
  expect_equal(unname(a$total), c(2, 1, 3, 0))
  expect_equal(colSums(a$increments), a$total)
  expect_equal(network::network.edgecount(g), 1)
  expect_error(multi_edge_change_stats(g, rbind(c(1, 3), c(3, 1)), rhs), "duplicate")
})

test_that("joint state enumeration equals whole-graph recomputation", {
  g <- cs_test_graph()
  tails <- c(1L, 2L, 3L); heads <- c(3L, 3L, 4L)
  rhs <- "edges + triangles + star(c(2,3))"
  tab <- .cs_state_table(g, tails, heads, rhs, 12)
  before <- .cs_change_model(g, rhs)$statistics()
  for (i in seq_len(nrow(tab$states))) {
    h <- network::network.copy(g)
    on <- which(tab$states[i, ])
    if (length(on)) network::add.edges(h, tails[on], heads[on])
    after <- .cs_change_model(h, rhs)$statistics()
    expect_equal(unname(tab$changes[i, ]), unname(after - before))
  }
  expect_error(.cs_state_table(g, tails, heads, rhs, 2), "exceeds max_candidates")
})

test_that("whole-update CS normalizes and reduces to collapsed Poisson at theta zero", {
  g <- cs_test_graph()
  params <- cs_test_params()
  tab <- .cs_state_table(g, c(1L, 2L, 3L), c(3L, 3L, 4L),
                         "edges + triangles + star(c(2,3))", 12)
  tab$ages <- c(.2, .7, 1.1)
  dist <- .cs_table_distribution(tab, params)
  expect_equal(sum(dist$probs), 1, tolerance = 1e-12)
  params$CS_params[] <- 0
  dist <- .cs_table_distribution(tab, params)
  r <- dist$reference_probs
  expected <- apply(tab$states, 1, function(z) prod(ifelse(z, r, 1-r)))
  expect_equal(dist$probs, expected, tolerance = 1e-12)
  params$CS_params[2] <- 1.4
  tilted <- .cs_table_distribution(tab, params)
  # A joint triangle feature changes odds beyond the independent-edge product.
  both <- 4L; only_a <- 2L; only_b <- 3L; neither <- 1L
  odds <- tilted$probs[both] * tilted$probs[neither] /
    (tilted$probs[only_a] * tilted$probs[only_b])
  expect_equal(odds, exp(1.4), tolerance = 1e-10)
})

test_that("CS generated, direct and cached probabilities agree including new params", {
  old <- cs_test_graph()
  params <- cs_test_params()
  rhs <- "edges + triangles + star(c(2,3))"
  set.seed(508)
  for (i in seq_len(12)) {
    drawn <- PMF_mark_CS(1, params, old, generate_mark = TRUE,
                         formula_RHS = rhs, truncation = 4)
    direct <- PMF_mark_CS(1, params, old, mark = drawn$mark_sample,
                          formula_RHS = rhs, truncation = 4)
    expect_equal(drawn$log_mark_sample_density, direct$log_mark_density, tolerance = 1e-11)
    expect_equal(direct$log_density_func(params), direct$log_mark_density)
    expect_true(network::network.size(drawn$mark_sample) > 4 ||
                  network::network.edgecount(drawn$mark_sample) > 1)
    changed <- params; changed$m <- .8; changed$node_lambda <- 1.1
    changed$CS_params[2] <- .2; changed$beta_edges <- .9
    fresh <- PMF_mark_CS(1, changed, old, mark = drawn$mark_sample,
                         formula_RHS = rhs, truncation = 4)
    expect_equal(direct$log_density_func(changed), fresh$log_mark_density, tolerance = 1e-11)
  }
  expect_equal(network::network.edgecount(old), 1)
})

test_that("nonempty conditioning normalizes a fixed-node CS support", {
  old <- cs_test_graph(); params <- cs_test_params(); params$node_lambda <- 0
  rhs <- "edges + triangles + star(c(2,3))"
  pairs <- t(utils::combn(1:4, 2))
  pairs <- pairs[!(pairs[, 1] == 1 & pairs[, 2] == 2), , drop = FALSE]
  d <- nrow(pairs); probs <- numeric(2^d)
  for (mask in 0:(2^d-1)) {
    on <- which(as.logical(intToBits(mask)[seq_len(d)]))
    h <- network::network.copy(old)
    if (length(on)) {
      network::add.edges(h, pairs[on,1], pairs[on,2])
      network::set.edge.attribute(h, "time", c(0, rep(1, length(on))))
    }
    q <- PMF_mark_CS(1, params, old, mark=h, formula_RHS=rhs, truncation=4)
    probs[mask+1] <- q$mark_density
  }
  expect_equal(probs[1], 0)
  expect_equal(sum(probs), 1, tolerance=1e-10)
  full <- network::network(matrix(1,4,4)-diag(4), directed=FALSE)
  network::set.vertex.attribute(full,"time",rep(0,4))
  network::set.edge.attribute(full,"time",rep(0,6))
  expect_error(PMF_mark_CS(1, params, full, generate_mark=TRUE,
                           formula_RHS=rhs, truncation=4), "No nonempty")
})

test_that("rare nonempty updates keep their support and cached environments stay small", {
  old <- network::network.initialize(3, directed = FALSE)
  network::set.vertex.attribute(old, "time", rep(0, 3))
  params <- cs_test_params(); params$node_lambda <- 0; params$CS_params <- -800
  set.seed(281)
  # Omitting truncation selects the joint-model default node window of four.
  q <- PMF_mark_CS(1, params, old, generate_mark = TRUE, formula_RHS = "edges")
  expect_equal(network::network.edgecount(q$mark_sample), 1)
  expect_equal(q$mark_density, 1/3, tolerance = 1e-12)
  expect_equal(q$edge_probs, rep(1/3, 3), tolerance = 1e-12)
  expect_equal(q$density_func(params), q$mark_density)
  expect_identical(parent.env(environment(q$log_density_func)), baseenv())
  expect_false(any(c("old", "mark", "model") %in% ls(environment(q$log_density_func))))
})

test_that("whole marks include birth statistics and require declared birth attributes", {
  old <- cs_test_graph(); params <- cs_test_params()
  params$vertex_categorical <- list(group = c(a = .4, b = .6))
  network::set.vertex.attribute(old, "group", rep("a", 4))
  born <- network::network.copy(old)
  network::add.vertices(born, 1)
  network::set.vertex.attribute(born, "time", c(rep(0,4), 1))
  network::set.vertex.attribute(born, "group", c(rep("a",4), "b"))
  rhs <- "edges + triangles + star(c(2,3))"
  q <- PMF_mark_CS(1, params, old, mark = born, formula_RHS = rhs, truncation = 4)
  before <- .cs_change_model(strip_vertex_attrs_for_ernm(old, rhs, params), rhs)$statistics()
  after <- .cs_change_model(strip_vertex_attrs_for_ernm(born, rhs, params), rhs)$statistics()
  expect_equal(q$mark_change_stats, after - before)
  network::delete.vertex.attribute(born, "group")
  expect_error(PMF_mark_CS(1, params, old, mark = born, formula_RHS = rhs), "categorical attribute")
})
