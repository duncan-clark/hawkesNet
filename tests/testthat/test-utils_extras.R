test_that("get_latest_times returns per-node latest activity", {
  el <- matrix(c(1, 2, 2, 3, 1, 3), ncol = 2, byrow = TRUE)
  net <- network::network(el, directed = FALSE)
  network::set.vertex.attribute(net, "time", c(0.1, 0.2, 0.3))
  network::set.edge.attribute(net, "time", c(0.1, 0.2, 0.5))

  lt <- get_latest_times(net)
  expect_length(lt, 3)
  # Node 1 has edges at 0.1 and 0.5; latest = 0.5

  expect_equal(lt[1], 0.5)
  # Node 3 has edges at 0.2 and 0.5; latest = 0.5
  expect_equal(lt[3], 0.5)
})

test_that("safe_parallel_lapply with 1 core falls back to lapply", {
  res <- hawkesNet:::safe_parallel_lapply(1:5, function(x) x^2, mc.cores = 1)
  expect_equal(res, as.list((1:5)^2))
})
