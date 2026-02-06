# Unit tests for normalize_vertex_categorical_probs

test_that("normalize_vertex_categorical_probs returns vector that sums to 1", {
  p <- normalize_vertex_categorical_probs(c(a = 0.25, b = 0.75))
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0))
  expect_named(p, c("a", "b"))
})

test_that("normalize_vertex_categorical_probs normalizes unnormalized positive vector", {
  p <- normalize_vertex_categorical_probs(c(1, 2, 3))
  expect_equal(sum(p), 1)
  expect_equal(p[1], 1/6)
  expect_equal(p[2], 2/6)
  expect_equal(p[3], 3/6)
})

test_that("normalize_vertex_categorical_probs handles zeros and small values with eps", {
  p <- normalize_vertex_categorical_probs(c(a = 0, b = 1, c = 0), eps = 1e-10)
  expect_equal(sum(p), 1)
  expect_true(all(p > 0))
  expect_true(p["b"] > p["a"])
  expect_true(p["b"] > p["c"])
})

test_that("normalize_vertex_categorical_probs returns NULL for NULL or empty input", {
  expect_null(normalize_vertex_categorical_probs(NULL))
  expect_null(normalize_vertex_categorical_probs(numeric(0)))
})

test_that("normalize_vertex_categorical_probs clamps NA/zero/negative and renormalizes", {
  # NA and 0 are replaced by eps then renormalized; negative clamped to eps
  p_na <- normalize_vertex_categorical_probs(c(NA, NA))
  expect_equal(sum(p_na), 1)
  expect_true(all(p_na > 0))
  p_zero <- normalize_vertex_categorical_probs(c(0, 0))
  expect_equal(sum(p_zero), 1)
  expect_true(all(p_zero > 0))
  p_neg <- normalize_vertex_categorical_probs(c(-1, 2))
  expect_equal(sum(p_neg), 1)
  expect_true(all(p_neg > 0))
})

test_that("normalize_vertex_categorical_probs preserves names", {
  p <- c(x = 1, y = 1, z = 2)
  out <- normalize_vertex_categorical_probs(p)
  expect_equal(names(out), c("x", "y", "z"))
})
