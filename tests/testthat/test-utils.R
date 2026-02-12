# Unit tests for utility functions: get_times, has_edge, events_to_net, filtration_to_net

test_that("get_times returns list with node_times, edge_times, times", {
  net <- network(matrix(c(1, 2), nrow = 1), directed = FALSE)
  set.vertex.attribute(net, "time", c(0.1, 0.2))
  set.edge.attribute(net, "time", 0.15)
  out <- get_times(net)
  expect_type(out, "list")
  expect_true("node_times" %in% names(out))
  expect_true("edge_times" %in% names(out))
  expect_true("times" %in% names(out))
  expect_equal(sort(unique(c(out$node_times, out$edge_times))), out$times)
})

test_that("get_times with single vertex has times from vertex only", {
  net <- network.initialize(1, directed = FALSE)
  set.vertex.attribute(net, "time", 0.5)
  out <- get_times(net)
  expect_equal(out$times, 0.5)
})

test_that("has_edge returns TRUE for existing edge in hash", {
  edge_hash <- hash::hash(keys = c("1-2", "2-3"), values = c(TRUE, TRUE))
  expect_true(has_edge(1L, 2L, edge_hash))
  expect_true(has_edge(2L, 3L, edge_hash))
  expect_false(has_edge(1L, 3L, edge_hash))
})

test_that("events_to_net builds network from event list", {
  events_list <- list(
    i = c(1L, 1L, 2L),
    j = c(2L, 3L, 3L),
    t = c(0.1, 0.2, 0.3)
  )
  net <- events_to_net(events_list)
  expect_s3_class(net, "network")
  expect_equal(network.size(net), 3)
  expect_equal(network.edgecount(net), 3)
  expect_equal(net %n% "n", 3)
})

test_that("events_to_net with directed = TRUE creates directed network", {
  events_list <- list(i = 1L, j = 2L, t = 0.1)
  net <- events_to_net(events_list, directed = TRUE)
  expect_true(is.directed(net))
})

test_that("filtration_to_net subsets by time", {
  events_list <- list(
    i = c(1L, 1L, 2L),
    j = c(2L, 3L, 3L),
    t = c(0.1, 0.2, 0.5)
  )
  net <- events_to_net(events_list)
  net_t <- filtration_to_net(net, 0.25, equals = FALSE)
  expect_true(network.size(net_t) <= 3)
  expect_true(network.edgecount(net_t) <= 3)
})

test_that("normalize_times_01 scales network times to [0, 1]", {
  events_list <- list(
    i = c(1L, 1L, 2L),
    j = c(2L, 3L, 3L),
    t = c(10, 20, 30)
  )
  net <- events_to_net(events_list)
  out <- normalize_times_01(net)
  expect_s3_class(out, "network")
  times_obj <- get_times(out)
  expect_true(all(times_obj$times >= 0 & times_obj$times <= 1))
  expect_equal(range(times_obj$times), c(0, 1))
})

test_that("normalize_times_01 handles single time value", {
  net <- network.initialize(1, directed = FALSE)
  set.vertex.attribute(net, "time", 5)
  out <- normalize_times_01(net, constant_value = 0)
  expect_equal(get.vertex.attribute(out, "time"), 0)
})
