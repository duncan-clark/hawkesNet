test_that("compensator_temporal_hawkes returns increasing values", {
  set.seed(1)
  realiz <- data.frame(t = sort(runif(50, 0, 10)))
  comp <- compensator_temporal_hawkes(c(1, 0.5, 0.2), realiz, c(0, 10), kernel = "exp")
  expect_type(comp, "double")
  expect_length(comp, 50)
  expect_true(all(diff(comp) >= 0))
})

test_that("ks_test_pval_temporal returns valid p-value", {
  set.seed(1)
  realiz <- data.frame(t = sort(runif(50, 0, 10)))
  pval <- ks_test_pval_temporal(realiz, c(0, 10), c(1, 0.5, 0.2), kernel = "exp")
  expect_type(pval, "double")
  expect_gte(pval, 0)
  expect_lte(pval, 1)
})

test_that("compensator_temporal_hawkes works with powerlaw kernel", {
  set.seed(1)
  realiz <- data.frame(t = sort(runif(50, 0, 10)))
  comp <- compensator_temporal_hawkes(c(1, 0.5, 1.5), realiz, c(0, 10), kernel = "powerlaw")
  expect_type(comp, "double")
  expect_length(comp, 50)
})
