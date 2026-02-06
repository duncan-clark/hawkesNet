# Unit tests for validate_point_process_params and related behavior

test_that("validate_point_process_params accepts valid minimal params", {
  expect_invisible(validate_point_process_params(list(mu = 0.1, beta_overall = 1, K = 0.5)))
  expect_invisible(validate_point_process_params(list(
    mu = 1, beta_overall = 2, K = 0.3, beta_edges = 0.5
  )))
})

test_that("validate_point_process_params accepts NULL or empty params", {
  expect_invisible(validate_point_process_params(NULL))
  expect_invisible(validate_point_process_params(list()))
})

test_that("validate_point_process_params errors when mu <= 0", {
  expect_error(
    validate_point_process_params(list(mu = 0, beta_overall = 1, K = 0.5)),
    "mu must be"
  )
  expect_error(
    validate_point_process_params(list(mu = -0.1, beta_overall = 1, K = 0.5)),
    "mu must be"
  )
})

test_that("validate_point_process_params errors when beta_overall <= 0", {
  expect_error(
    validate_point_process_params(list(mu = 0.1, beta_overall = 0, K = 0.5)),
    "beta_overall must be"
  )
})

test_that("validate_point_process_params errors when K <= 0", {
  expect_error(
    validate_point_process_params(list(mu = 0.1, beta_overall = 1, K = 0)),
    "K must be"
  )
})

test_that("validate_point_process_params errors when beta_edges <= 0 when present", {
  expect_error(
    validate_point_process_params(list(mu = 0.1, beta_overall = 1, K = 0.5, beta_edges = 0)),
    "beta_edges must be"
  )
})

test_that("validate_point_process_params errors when node_lambda <= 0 when present", {
  expect_error(
    validate_point_process_params(list(mu = 0.1, beta_overall = 1, K = 0.5, node_lambda = 0)),
    "node_lambda must be"
  )
})

test_that("validate_point_process_params errors when mu is NA or non-finite", {
  expect_error(
    validate_point_process_params(list(mu = NA_real_, beta_overall = 1, K = 0.5)),
    "finite"
  )
  expect_error(
    validate_point_process_params(list(mu = Inf, beta_overall = 1, K = 0.5)),
    "finite"
  )
})

test_that("validate_point_process_params errors when mu is not a scalar", {
  expect_error(
    validate_point_process_params(list(mu = c(0.1, 0.2), beta_overall = 1, K = 0.5)),
    "numeric scalar"
  )
})

test_that("validate_point_process_params accepts valid vertex_categorical", {
  expect_invisible(validate_point_process_params(list(
    mu = 0.1, beta_overall = 1, K = 0.5,
    vertex_categorical = list(attr1 = c(a = 0.5, b = 0.5))
  )))
  expect_invisible(validate_point_process_params(list(
    mu = 0.1, beta_overall = 1, K = 0.5,
    vertex_categorical = list(x = c(1, 0, 0))
  )))
})

test_that("validate_point_process_params errors when vertex_categorical has negative probs", {
  expect_error(
    validate_point_process_params(list(
      mu = 0.1, beta_overall = 1, K = 0.5,
      vertex_categorical = list(attr1 = c(a = 0.7, b = -0.1, c = 0.4))
    )),
    "non-negative"
  )
})

test_that("validate_point_process_params errors when vertex_categorical does not sum to 1", {
  expect_invisible(validate_point_process_params(list(
    mu = 0.1, beta_overall = 1, K = 0.5,
    vertex_categorical = list(attr1 = c(a = 0.5, b = 0.5))
  )))
  expect_error(
    validate_point_process_params(list(
      mu = 0.1, beta_overall = 1, K = 0.5,
      vertex_categorical = list(attr1 = c(a = 0.5, b = 0.6))
    )),
    "sum to 1"
  )
})

test_that("validate_point_process_params errors when vertex_categorical is empty or non-numeric", {
  expect_error(
    validate_point_process_params(list(
      mu = 0.1, beta_overall = 1, K = 0.5,
      vertex_categorical = list(attr1 = numeric(0))
    )),
    "non-empty numeric"
  )
})
