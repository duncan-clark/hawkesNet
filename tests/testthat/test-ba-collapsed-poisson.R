ba_test_history <- function(n = 3L, edges = matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE)) {
  net <- network::network.initialize(n, directed = FALSE)
  network::set.vertex.attribute(net, "time", seq_len(n) / 10)
  if (nrow(edges)) {
    network::add.edges(net, edges[, 1L], edges[, 2L], names.eval = "time",
                       vals.eval = rep(0.5, nrow(edges)))
  }
  net
}

ba_test_update <- function(net, targets = integer(0), time = 1) {
  out <- network::network.copy(net)
  n <- network::network.size(net)
  network::add.vertices(out, 1L)
  network::set.vertex.attribute(out, "time", c(network::get.vertex.attribute(net, "time"), time))
  if (length(targets)) {
    network::add.edges(out, rep(n + 1L, length(targets)), targets,
                       names.eval = "time", vals.eval = rep(time, length(targets)))
  }
  out
}

test_that("BA complete-update probabilities normalize and match Poisson splitting", {
  history <- ba_test_history()
  params <- list(m = 2, beta_edges = 0.4)
  probabilities <- vapply(0:7, function(mask) {
    targets <- which(as.logical(intToBits(mask)[1:3]))
    ans <- PMF_mark_BA(1, params, history, mark = ba_test_update(history, targets))
    expected <- prod(ifelse(seq_len(3) %in% targets, ans$edge_probs, 1 - ans$edge_probs))
    expect_equal(ans$mark_density, expected, tolerance = 1e-12)
    expect_equal(ans$log_density_func(params), ans$log_mark_density, tolerance = 1e-12)
    expect_equal(ans$edge_probs, -expm1(-params$m * ans$attachment_probs))
    ans$mark_density
  }, numeric(1))
  expect_equal(sum(probabilities), 1, tolerance = 1e-12)
  expect_equal(probabilities[1], exp(-params$m), tolerance = 1e-12)
})

test_that("BA generated, direct, and cached densities agree including attributes", {
  history <- ba_test_history()
  params <- list(m = 1.7, beta_edges = 0.3,
                 vertex_categorical = list(group = c(0.25, 0.75)),
                 vertex_categorical_levels = list(group = c("A", "B")))
  set.seed(42)
  for (i in seq_len(20)) {
    generated <- PMF_mark_BA(1, params, history, generate_mark = TRUE)
    direct <- PMF_mark_BA(1, params, history, mark = generated$mark_sample)
    expect_equal(generated$log_mark_sample_density, direct$log_mark_density)
    expect_equal(generated$log_mark_density, direct$log_mark_density)
    expect_equal(generated$log_density_func(params), direct$log_mark_density)
    expect_equal(network::network.size(generated$mark_sample), 4)
    expect_equal(network::network.size(history), 3)
    alt <- params
    alt$m <- 0.8
    alt$beta_edges <- 1.2
    alt$vertex_categorical$group <- c(0.6, 0.4)
    alt_direct <- PMF_mark_BA(1, alt, history, mark = generated$mark_sample)
    expect_equal(direct$log_density_func(alt), alt_direct$log_mark_density)
  }
})

test_that("BA sampler frequencies follow the collapsed law", {
  history <- ba_test_history()
  params <- list(m = 2, beta_edges = 0.4)
  reference <- PMF_mark_BA(1, params, history, mark = ba_test_update(history))
  set.seed(321)
  samples <- replicate(600, {
    out <- PMF_mark_BA(1, params, history, generate_mark = TRUE)$mark_sample
    vapply(seq_len(3), function(target) length(network::get.edgeIDs(out, 4, target)) > 0, logical(1))
  })
  expect_lt(max(abs(rowMeans(samples) - reference$edge_probs)), 0.07)
  expect_lt(abs(mean(colSums(samples) == 0) - exp(-params$m)), 0.05)
})

test_that("BA zero, one, and zero-weight candidate cases use the same law", {
  empty_edges <- matrix(integer(0), ncol = 2)
  params <- list(m = 2, beta_edges = 0.4)
  empty <- ba_test_history(0L, empty_edges)
  seed <- PMF_mark_BA(1, params, empty, generate_mark = TRUE)
  expect_equal(network::network.size(seed$mark_sample), 1)
  expect_equal(seed$log_mark_sample_density, 0)
  expect_length(seed$edge_probs, 0)
  expect_equal(PMF_mark_BA(1, params, empty, mark = seed$mark_sample)$log_mark_density, 0)

  single <- ba_test_history(1L, empty_edges)
  absent <- PMF_mark_BA(1, params, single, mark = ba_test_update(single))
  present <- PMF_mark_BA(1, params, single, mark = ba_test_update(single, 1L))
  expect_equal(absent$mark_density, exp(-2))
  expect_equal(present$mark_density, -expm1(-2))
  expect_equal(absent$mark_density + present$mark_density, 1)

  isolated <- ba_test_history(3L, empty_edges)
  uniform <- PMF_mark_BA(1, params, isolated, mark = ba_test_update(isolated))
  expect_equal(uniform$attachment_probs, rep(1 / 3, 3))
  expect_equal(uniform$log_density_func(params), -2)
  no_targets <- PMF_mark_BA(1, params, isolated, generate_mark = TRUE, truncation = 0)
  expect_equal(no_targets$log_mark_sample_density, 0)
  expect_equal(network::network.edgecount(no_targets$mark_sample), 0)
  zero_m <- PMF_mark_BA(1, list(m = 0, beta_edges = 0.4), isolated, generate_mark = TRUE)
  expect_equal(zero_m$log_mark_sample_density, 0)
  expect_equal(network::network.edgecount(zero_m$mark_sample), 0)
})

test_that("BA truncation, activity decay, and impossible marks have consistent support", {
  history <- ba_test_history()
  params <- list(m = 2, beta_edges = 1)
  allowed <- ba_test_update(history, 3L)
  out <- PMF_mark_BA(1, params, history, mark = allowed, truncation = 1)
  expect_equal(out$candidate_heads, 3L)
  expect_equal(out$attachment_probs, 1)
  expect_equal(out$mark_density, -expm1(-2))
  forbidden <- PMF_mark_BA(1, params, history, mark = ba_test_update(history, 1L), truncation = 1)
  expect_equal(forbidden$log_mark_density, -Inf)
  expect_equal(forbidden$log_density_func(params), -Inf)
  missing_node <- PMF_mark_BA(1, params, history, mark = history)
  expect_equal(missing_node$log_mark_density, -Inf)

  activity <- PMF_mark_BA(1, params, history, mark = allowed, mark_decay = "activity")
  expect_equal(activity$attachment_probs, c(0.25, 0.5, 0.25))
  expect_equal(activity$log_density_func(params), activity$log_mark_density)
  far_future <- PMF_mark_BA(1e7, params, history, generate_mark = TRUE)
  expect_equal(sum(far_future$attachment_probs), 1)
  expect_true(is.finite(far_future$log_mark_sample_density))
})
