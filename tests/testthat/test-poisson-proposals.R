test_that("overlapping proposals have a normalized unordered union law", {
  groups <- list("a", "b", "c", c("a", "b"), c("b", "c"))
  rates <- c(.2, .3, .4, .5, .6)
  for (nonempty in c(FALSE, TRUE)) {
    law <- poisson_proposals(groups, rates, nonempty)
    sets <- lapply(0:7, function(mask) letters[1:3][as.logical(intToBits(mask)[1:3])])
    actual <- vapply(sets, function(s) exp(logprob_update(law, s)), numeric(1))
    # Independently enumerate all 32 candidate-group selections.
    expected <- numeric(8)
    for (mask in 0:31) {
      selected <- as.logical(intToBits(mask)[1:5])
      edits <- unique(unlist(groups[selected], use.names = FALSE))
      state <- 1L + sum(2L^(match(edits, letters[1:3]) - 1L))
      probability <- prod(ifelse(selected, -expm1(-rates), exp(-rates)))
      expected[state] <- expected[state] + probability
    }
    if (nonempty) { expected[1] <- 0; expected <- expected / sum(expected) }
    expect_equal(actual, expected, tolerance = 1e-12)
    set.seed(812)
    sampled <- replicate(10000, paste(sample_update(law), collapse = ""))
    observed <- vapply(sets, function(s) mean(sampled == paste(s, collapse = "")), numeric(1))
    expect_lt(max(abs(observed - expected)), .018)
  }
})

test_that("rare unions stay finite and computation guards do not alter support", {
  law <- poisson_proposals(list(c("a", "b"), c("b", "c")), c(1e-200, 2e-200))
  expect_equal(logprob_update(law, c("a", "b", "c")), log(2) + 2 * log(1e-200))
  expect_error(logprob_update(law, c("a", "b", "c"), max_component = 2), "refused")
  independent <- poisson_proposals(as.list(as.character(1:40)), rep(.1, 40))
  expect_true(is.finite(logprob_update(independent, as.character(1:40))))
  tiny <- poisson_proposals(list("a", "b"), c(1e-200, 1e-200), TRUE)
  expect_length(sample_update(tiny), 1)
  expect_error(poisson_proposals(list("a"), 0, TRUE), "positive")
  expect_equal(logprob_update(independent, "absent"), -Inf)
})
