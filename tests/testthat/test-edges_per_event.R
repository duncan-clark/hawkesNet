test_that("get_edges_per_event returns correct counts", {
  el <- matrix(c(1, 2, 1, 3, 2, 3), ncol = 2, byrow = TRUE)
  net <- network::network(el, directed = FALSE)
  network::set.edge.attribute(net, "time", c(0.1, 0.1, 0.5))

  df <- get_edges_per_event(net)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("time", "n_edges"))
  expect_equal(nrow(df), 2)
  expect_equal(df$n_edges[df$time == 0.1], 2L)
  expect_equal(df$n_edges[df$time == 0.5], 1L)
})

test_that("get_edges_per_event handles empty network", {
  net <- network::network.initialize(3, directed = FALSE)
  df <- get_edges_per_event(net)
  expect_equal(nrow(df), 0)
})

test_that("get_edges_per_event rejects non-network", {
  expect_error(get_edges_per_event(data.frame()), "network object")
})
