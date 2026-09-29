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

test_that("validate_point_process_params allows K = 0 but rejects negative excitation", {
  expect_invisible(validate_point_process_params(list(mu = 0.1, beta_overall = 1, K = 0)))
  expect_error(
    validate_point_process_params(list(mu = 0.1, beta_overall = 1, K = -0.1)),
    "K must be"
  )
})

test_that("validate_point_process_params allows zero edge decay but rejects negative decay", {
  expect_invisible(validate_point_process_params(list(mu = 0.1, beta_overall = 1, K = 0.5, beta_edges = 0)))
  expect_error(
    validate_point_process_params(list(mu = 0.1, beta_overall = 1, K = 0.5, beta_edges = -0.1)),
    "beta_edges must be"
  )
})

test_that("both validators accept zero attempts and reject negative attempts", {
  params <- list(mu = 1, beta_overall = 1, K = 0, beta_edges = 0, m = 0)
  expect_invisible(validate_point_process_params(params))
  expect_true(point_process_params_valid(params))
  params$m <- -0.1
  expect_error(validate_point_process_params(params), "m must be")
  expect_false(point_process_params_valid(params))
})

test_that("validators allow no node births but reject negative node rates", {
  params <- list(mu = 0.1, beta_overall = 1, K = 0.5, node_lambda = 0)
  expect_invisible(validate_point_process_params(params))
  expect_true(point_process_params_valid(params))
  params$node_lambda <- -0.1
  expect_error(
    validate_point_process_params(params),
    "node_lambda must be"
  )
  expect_false(point_process_params_valid(params))
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

test_that("validate_point_process_params accepts valid vertex_categorical (n-1 parametrization)", {
  # n-1 params: 2 levels, 1 param, sum < 1
  expect_invisible(validate_point_process_params(list(
    mu = 0.1, beta_overall = 1, K = 0.5,
    vertex_categorical = list(attr1 = c(a = 0.4))
  )))
  # 3 levels, 2 params, sum < 1
  expect_invisible(validate_point_process_params(list(
    mu = 0.1, beta_overall = 1, K = 0.5,
    vertex_categorical = list(x = c(a = 0.3, b = 0.3))
  )))
})

test_that("validate_point_process_params errors when vertex_categorical has negative probs", {
  expect_error(
    validate_point_process_params(list(
      mu = 0.1, beta_overall = 1, K = 0.5,
      vertex_categorical = list(attr1 = c(a = 0.7, b = -0.1))
    )),
    "non-negative"
  )
})

test_that("validate_point_process_params errors when vertex_categorical sum >= 1 (n-1)", {
  # sum == 1 is invalid (reference level would get 0)
  expect_error(
    validate_point_process_params(list(
      mu = 0.1, beta_overall = 1, K = 0.5,
      vertex_categorical = list(attr1 = c(a = 0.5, b = 0.5))
    )),
    "sum in \\(0, 1\\)"
  )
  # sum > 1 is invalid
  expect_error(
    validate_point_process_params(list(
      mu = 0.1, beta_overall = 1, K = 0.5,
      vertex_categorical = list(attr1 = c(a = 0.5, b = 0.6))
    )),
    "sum in \\(0, 1\\)"
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
