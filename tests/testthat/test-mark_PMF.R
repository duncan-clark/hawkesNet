# Unit tests for mark PMF helpers and expected params

test_that("expected_params_PMF_mark_BA returns required beta_edges and m", {
  out <- expected_params_PMF_mark_BA()
  expect_type(out, "list")
  expect_true("required" %in% names(out))
  expect_true("beta_edges" %in% out$required)
  expect_true("m" %in% out$required)
})

test_that("expected_params_PMF_mark_CS returns required and CS_params_length", {
  out <- expected_params_PMF_mark_CS(NULL, NULL)
  expect_true("required" %in% names(out))
  expect_true("CS_params_length" %in% names(out))
  expect_true("node_lambda" %in% out$required)
  expect_true("CS_params" %in% out$required)
  expect_true("beta_edges" %in% out$required)
})

test_that("expected_params_PMF_mark_CS with formula and net returns CS_params_length when possible", {
  net <- network::network(matrix(c(1, 2), nrow = 1), directed = FALSE)
  network::set.vertex.attribute(net, "time", c(0, 0))
  network::set.edge.attribute(net, "time", 0)
  out <- expected_params_PMF_mark_CS(net, "edges")
  expect_true("CS_params_length" %in% names(out))
})

test_that("validate_params_for_PMF errors when required param missing for BA", {
  expect_error(
    validate_params_for_PMF(
      list(mu = 1, beta_overall = 1, K = 0.5),
      PMF_mark_BA,
      mark_filtration = NULL
    ),
    "beta_edges|m"
  )
})

test_that("validate_params_for_PMF passes when required params present for BA", {
  expect_invisible(validate_params_for_PMF(
    list(mu = 1, beta_overall = 1, K = 0.5, beta_edges = 0.5, m = 1),
    PMF_mark_BA,
    mark_filtration = NULL
  ))
})
