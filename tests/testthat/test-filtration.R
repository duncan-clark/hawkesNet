# tests/testthat/test-filtration_to_net.R

test_that("edge_ids_where maps to active EIDs after deletions", {
  suppressMessages({
    library(network)
  })
  
  # Build a simple directed network: 1->2 (t=1), 2->3 (t=2), 3->4 (t=3)
  net <- network.initialize(4, directed = TRUE)
  add.edges(net, tail = c(1, 2, 3), head = c(2, 3, 4))
  set.edge.attribute(net, "time", c(1, 2, 3))
  set.vertex.attribute(net, "time", c(0, 1, 2, 3))
  
  # Create a hole in net$mel: delete the middle edge (EID 2)
  delete.edges(net, 2L)
  
  # Now active EIDs should be {1, 3}. edge_ids_where must return TRUE EIDs.
  e_ge3 <- edge_ids_where(net, "time", function(x) x >= 3)
  e_le1 <- edge_ids_where(net, "time", function(x) x <= 1)
  
  expect_type(e_ge3, "integer")
  expect_type(e_le1, "integer")
  
  expect_equal(e_ge3, 3L)   # only the (original) 3rd edge remains with time 3
  expect_equal(e_le1, 1L)   # only the (original) 1st edge remains with time 1
})

test_that("filtration_to_net(equals = TRUE) keeps <= t and drops later times", {
  suppressMessages({
    library(network)
  })
  
  # Build fresh network: edges times 1,2,3; vertex times 0,1,2,3
  net <- network.initialize(4, directed = TRUE)
  add.edges(net, tail = c(1, 2, 3), head = c(2, 3, 4))
  set.edge.attribute(net, "time", c(1, 2, 3))
  set.vertex.attribute(net, "time", c(0, 1, 2, 3))
  set.vertex.attribute(net, "vertex.names", letters[1:4])
  
  t <- 1
  out <- filtration_to_net(net, t = t, equals = TRUE)
  
  # Should keep vertices with time <= 1: {1,2}; remove others
  expect_equal(network.size(out), 2L)
  expect_true(all(get.vertex.attribute(out, "time") <= t))
  
  # Should keep edges with time <= 1: only 1->2 (time = 1)
  expect_equal(network.edgecount(out), 1L)
  expect_equal(get.edge.attribute(out, "time"), 1)
  
  # vertex.names should be removed
  expect_false("vertex.names" %in% list.vertex.attributes(out))
})

test_that("filtration_to_net(equals = FALSE) also removes the latest remaining time slice", {
  skip_on_cran()
  suppressMessages({
    library(network)
  })
  
  # Build fresh network again
  net <- network.initialize(4, directed = TRUE)
  add.edges(net, tail = c(1, 2, 3), head = c(2, 3, 4))
  set.edge.attribute(net, "time", c(1, 2, 3))
  set.vertex.attribute(net, "time", c(0, 1, 2, 3))
  
  t <- 1
  out <- filtration_to_net(net, t = t, equals = FALSE)
  
  # First pass removes > t (edges at 2,3 and vertices at 2,3),
  # then removes the max remaining time (= 1) among edges/vertices.
  # So we expect only the vertex with time 0 to remain, and no edges.
  expect_equal(network.size(out), 1L)
  expect_equal(network.edgecount(out), 0L)
  expect_true(all(get.vertex.attribute(out, "time") == 0))
  
  # vertex.names should be removed if present
  expect_false("vertex.names" %in% list.vertex.attributes(out))
})

test_that("edge_ids_where returns integer(0) when attribute missing or empty", {
  skip_on_cran()
  suppressMessages({
    library(network)
  })
  
  net <- network.initialize(3, directed = TRUE)
  add.edges(net, tail = 1, head = 2)
  
  # No 'time' attribute set yet
  expect_equal(edge_ids_where(net, "time", function(x) x > 0), integer(0))
  
  # Set time but filter matches nothing
  set.edge.attribute(net, "time", 0)
  expect_equal(edge_ids_where(net, "time", function(x) x > 0), integer(0))
})