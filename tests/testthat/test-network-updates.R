test_that("turnover updates revise past amplitudes before the next event", {
  graph <- matrix(0, 3, 3); graph[1, 2] <- graph[2, 1] <- 1
  path <- network_update_data(c(1, 2), list(
    list(births = 0, anchor = 1, add = 3, delete = integer()),
    list(births = 0, anchor = 2, add = integer(), delete = 1)), graph, 3)
  for (timing in c("M0", "M1", "M2")) {
    model <- network_update_model(timing = timing, threshold = 2)
    result <- loglik_network_updates(path, model)
    gb <- if (timing == "M0") 0 else .7
    gd <- if (timing == "M2") 1.2 else 0
    expected_rates <- c(2, 2 + .25 * exp(gb - 1))
    expected_integral <- 6 + .25 * (exp(gb) + exp(gb - gd - 1) + 1) * (1 - exp(-1))
    expect_equal(unname(result$intensities), expected_rates)
    expect_equal(unname(result$compensator), expected_integral)
    expect_equal(unname(result$timing), sum(log(expected_rates)) - expected_integral)
  }
})

test_that("CS candidate changes match independent graph statistic differences", {
  stats <- function(g) {
    d <- rowSums(g)
    c(sum(diag(g %*% g %*% g)) / 6, sum(choose(d, 2)) / 6, sum(choose(d, 3)) / 18)
  }
  for (mask in 0:63) {
    graph <- matrix(FALSE, 4, 4)
    graph[upper.tri(graph)] <- as.logical(intToBits(mask)[1:6]); graph <- graph | t(graph)
    candidates <- hawkesNet:::.update_candidates(graph, 1, "CS-1")
    for (j in seq_along(candidates$targets)) {
      after <- graph; v <- candidates$targets[j]
      after[1, v] <- after[v, 1] <- !graph[1, v]
      expect_equal(as.numeric(candidates$x[j, ]), stats(after) - stats(graph))
    }
  }
})

test_that("reference simulation replays and timing nests exactly", {
  for (mark in c("PA", "CS-1")) {
    model <- network_update_model(mark, "M2")
    set.seed(812); path <- simulate_network_updates(model, 8)
    result <- loglik_network_updates(path, model)
    expect_true(is.finite(result$loglik))
    expect_equal(result$loglik, result$timing + result$marks)
    expect_true(all(vapply(path$updates, function(u) length(c(u$add, u$delete)) > 0, logical(1))))
    zero <- model; zero$parameters[c("gamma_B", "gamma_Delta")] <- 0
    m0 <- network_update_model(mark, "M0")
    expect_equal(loglik_network_updates(path, zero), loglik_network_updates(path, m0))
    m1 <- model; m1$timing <- "M1"; m1$parameters <- m1$parameters[names(m1$parameters) != "gamma_Delta"]
    m2zero <- model; m2zero$parameters["gamma_Delta"] <- 0
    expect_equal(loglik_network_updates(path, m1), loglik_network_updates(path, m2zero))
    fit <- fit_network_updates(path, m0, fixed = setdiff(names(m0$parameters), "mu"))
    expect_equal(fit$convergence, 0L)
    expect_gte(fit$loglik + 1e-8, loglik_network_updates(path, m0)$loglik)
    expect_true(is.finite(fit$standard_errors["mu"]))
  }
})

test_that("path validation and guards fail visibly", {
  m <- network_update_model()
  expect_error(network_update_data(c(1, 1), list(list(), list()), diag(2), 2), "initial")
  expect_error(network_update_data(c(1, 1), list(list(), list()), matrix(c(0,1,1,0),2), 2), "increasing")
  p <- network_update_data(1, list(list(births=0,anchor=1,add=2,delete=integer())), matrix(c(0,1,1,0),2), 2)
  expect_error(loglik_network_updates(p, m), "inadmissible")
  set.seed(1)
  expect_error(simulate_network_updates(m, 100, max_events=1), "max_events")
  expect_error(simulate_network_updates(network_update_model(timing="M2",summary="log_degree"), 1), "non-explosion")
})

test_that("unchanging anchor summaries flag unidentified later feedback", {
  m <- network_update_model(timing="M1", threshold=10000)
  set.seed(82); path <- simulate_network_updates(m, 8)
  result <- feedback_score(path, m)
  expect_equal(unname(result$score), 0)
  expect_match(result$reason, "rank-deficient")
  empty <- network_update_data(numeric(), list(), matrix(c(0,1,1,0),2), 3)
  expect_equal(unname(loglik_network_updates(empty, m)$loglik), -6)
})
