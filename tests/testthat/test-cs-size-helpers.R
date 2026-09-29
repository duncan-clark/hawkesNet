test_that("conditioned Poisson counts are exact and stable at boundaries", {
  for (m in c(0, 1e-300, .7, 30, 1e300)) {
    lp <- .cs_size_log_count(m, 5L)
    expect_equal(sum(exp(lp)), 1, tolerance = 1e-12)
    expect_false(anyNA(lp))
  }
  expect_equal(.cs_size_log_count(0, 3), c(0, -Inf, -Inf, -Inf))
  expect_equal(.cs_size_log_count(100, 0), 0)
  expect_equal(exp(.cs_size_log_count(2.3, 5)),
               stats::dpois(0:5, 2.3) / stats::ppois(5, 2.3))
  expect_error(.cs_size_log_count(-1, 3), "nonnegative")
  expect_error(.cs_size_log_count(1, 1.5), "integer")
})

test_that("fixed-size tables contain only the requested unordered sets", {
  g <- network::network.initialize(4, directed = FALSE)
  network::add.edge(g, 1, 2)
  network::set.vertex.attribute(g, "time", c(0, .1, .2, .3))
  network::set.edge.attribute(g, "time", 0)
  context <- .cs_size_context(g, g, 1, list(),
    "edges + triangles + star(c(2,3))", 4, "node_entrance", FALSE, 12)
  rhs <- "edges + triangles + star(c(2,3))"
  before <- .cs_change_model(g, rhs)$statistics()
  for (k in c(0, 2, 5)) {
    tab <- .cs_size_state_table(context, k, rhs)
    expect_equal(nrow(tab$states), choose(5, k))
    expect_equal(rowSums(tab$states), rep(k, choose(5, k)))
    for (i in seq_len(nrow(tab$states))) {
      copy <- network::network.copy(g)
      selected <- which(tab$states[i, ])
      if (length(selected)) {
        network::add.edges(copy, tab$tails[selected], tab$heads[selected])
      }
      expect_equal(unname(tab$changes[i, ]),
                   unname(.cs_change_model(copy, rhs)$statistics() - before))
    }
  }
  expect_equal(network::network.edgecount(g), 1L)
  expect_equal(network::get.vertex.attribute(g, "time"), c(0, .1, .2, .3))
  expect_equal(network::get.edge.attribute(g, "time"), 0)
})

test_that("CS-2 normalizes tied enormous finite whole-update scores", {
  tab <- list(changes = matrix(c(0, 1, 1), ncol = 1),
              summed_ages = c(0, 0, 0))
  dist <- .cs_size_distribution(tab, list(CS_params = 1e308, beta_edges = 0))
  expect_equal(dist$probs, c(0, .5, .5))
  expect_equal(sum(dist$probs), 1)
})

test_that("CS-2 handles empty candidates and impossible nonempty states", {
  old <- network::network.initialize(1, directed = FALSE)
  network::set.vertex.attribute(old, "time", 0)
  params <- list(m = 1, node_lambda = 0, beta_edges = .5, CS_params = 0)
  q <- PMF_mark_CS(1, params, old, mark = old, formula_RHS = "edges",
                   cs_mode = "size_conditional", condition_nonempty = FALSE)
  expect_equal(q$mark_density, 1)
  expect_equal(q$edge_probs, numeric(0))
  expect_equal(q$edge_count_probs, 1)
  expect_equal(q$n_mark_states, 1)
  expect_error(PMF_mark_CS(1, params, old, generate_mark = TRUE,
    formula_RHS = "edges", cs_mode = "size_conditional"), "No nonempty")
  impossible <- PMF_mark_CS(1, params, old, mark = old,
    formula_RHS = "edges", cs_mode = "size_conditional")
  expect_equal(impossible$log_mark_density, -Inf)
  params$node_lambda <- .2
  params$m <- 0
  set.seed(1443)
  born <- PMF_mark_CS(1, params, old, generate_mark = TRUE,
    formula_RHS = "edges", cs_mode = "size_conditional")
  expect_gt(network::network.size(born$mark_sample), 1)
  expect_equal(born$edge_count, 0)
  expect_equal(born$mark_density, born$density_func(params))
})

test_that("CS-2 caches support changed parameters and rare nonempty updates", {
  old <- network::network.initialize(3, directed = FALSE)
  network::set.vertex.attribute(old, "time", rep(0, 3))
  params <- list(m = 1e-300, node_lambda = 0, beta_edges = .5, CS_params = 0)
  set.seed(17)
  q <- PMF_mark_CS(1, params, old, generate_mark = TRUE,
    formula_RHS = "edges", cs_mode = "size_conditional")
  expect_equal(q$edge_count, 1)
  expect_equal(q$mark_density, 1/3, tolerance = 1e-12)
  changed <- params
  changed$m <- 2
  changed$node_lambda <- .3
  changed$CS_params <- 1e300
  fresh <- PMF_mark_CS(1, changed, old, mark = q$mark_sample,
    formula_RHS = "edges", cs_mode = "size_conditional")
  expect_equal(q$log_density_func(changed), fresh$log_mark_density, tolerance = 1e-12)
  expect_identical(parent.env(environment(q$log_density_func)), baseenv())
  expect_false(any(vapply(as.list(environment(q$log_density_func)),
                         network::is.network, logical(1))))
  expect_false(any(vapply(as.list(environment(q$log_density_func))$table_local,
                         network::is.network, logical(1))))
  params$m <- 0
  params$node_lambda <- 1e-300
  born <- PMF_mark_CS(1, params, old, generate_mark = TRUE,
    formula_RHS = "edges", cs_mode = "size_conditional")
  expect_equal(network::network.size(born$mark_sample), 4)
  expect_equal(born$edge_count, 0)
  expect_equal(born$mark_density, 1, tolerance = 1e-12)
})

test_that("CS-2 rejects deletions and marks outside the candidate support", {
  old <- network::network.initialize(4, directed = FALSE)
  network::set.vertex.attribute(old, "time", rep(0, 4))
  network::add.edge(old, 1, 2)
  network::set.edge.attribute(old, "time", 0)
  params <- list(m = 1, node_lambda = 0, beta_edges = .5, CS_params = 0)
  deleted <- network::network.copy(old)
  network::delete.edges(deleted, network::get.edgeIDs(deleted, 1, 2))
  q <- PMF_mark_CS(1, params, old, mark = deleted, formula_RHS = "edges",
                   cs_mode = "size_conditional", condition_nonempty = FALSE)
  expect_equal(q$log_mark_density, -Inf)
  outside <- network::network.copy(old)
  network::add.edge(outside, 1, 3)
  q <- PMF_mark_CS(1, params, old, mark = outside, formula_RHS = "edges",
                   truncation = 2, cs_mode = "size_conditional")
  expect_equal(q$log_mark_density, -Inf)
  multiple <- network::network.initialize(4, directed = FALSE, multiple = TRUE)
  network::set.vertex.attribute(multiple, "time", rep(0, 4))
  network::add.edges(multiple, c(1, 1, 1), c(2, 3, 3))
  q <- PMF_mark_CS(1, params, old, mark = multiple, formula_RHS = "edges",
                   cs_mode = "size_conditional")
  expect_equal(q$log_mark_density, -Inf)
})
