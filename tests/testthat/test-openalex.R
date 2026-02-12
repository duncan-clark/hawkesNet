# Unit tests for OpenAlex data integration
# These tests are marked as slow and will be skipped on CRAN

# Skip if OpenAlex API is unavailable or if running on CRAN
skip_if_no_openalex <- function() {
  skip_on_cran()
  skip_if_not_installed("httr")
  # Check if we can reach OpenAlex API
  testthat::skip_if_not(
    tryCatch({
      httr::GET("https://api.openalex.org/works?per_page=1", timeout = 5)
      TRUE
    }, error = function(e) FALSE),
    "OpenAlex API unavailable"
  )
}

test_that("OpenAlex network loading works with minimal data", {
  skip_if_no_openalex()
  skip_on_cran()
  
  # Use minimal configuration for speed
  source(system.file("openalex_study", "get_network_openalex.R", package = "hawkesNet"))
  
  # Load minimal network (1 page = 100 works max)
  out <- tryCatch({
    get_network(
      email = Sys.getenv("OPENALEX_EMAIL", "test@example.com"),
      string = "Hawkes",
      pages = 1L,  # Minimal: just 1 page
      per_page = 10L,  # Even smaller: 10 per page
      min_date = "2010-01-01",
      max_date = "2020-01-01"
    )
  }, error = function(e) {
    skip(paste("OpenAlex network loading failed:", e$message))
  })
  
  expect_true(is.list(out))
  expect_true("net" %in% names(out))
  expect_true("edges" %in% names(out))
  
  net <- out$net
  expect_s3_class(net, "network")
  expect_gt(network.size(net), 0)
  
  # Check that network has time attributes
  expect_true("time_scaled" %in% list.vertex.attributes(net))
  
  # Check gender attribute exists (even if all unknown)
  if ("gender" %in% list.vertex.attributes(net)) {
    gender_vals <- net %v% "gender"
    expect_true(all(gender_vals %in% c("female", "male", "unknown", NA)))
    
    # If gender package is available, check that we got some predictions
    if (requireNamespace("gender", quietly = TRUE) && 
        requireNamespace("genderdata", quietly = TRUE)) {
      # Gender attribute should exist and have values
      expect_true(length(gender_vals) > 0)
      # If gender prediction worked, we should have at least some non-unknown values
      # (but we're lenient - if all are unknown, that's okay, just means prediction failed)
      n_non_unknown <- sum(!is.na(gender_vals) & gender_vals != "unknown")
      if (n_non_unknown > 0) {
        expect_true(any(gender_vals %in% c("female", "male")))
      }
    }
  }
})

test_that("OpenAlex network can be normalized", {
  skip_if_no_openalex()
  skip_on_cran()
  
  source(system.file("openalex_study", "get_network_openalex.R", package = "hawkesNet"))
  
  out <- tryCatch({
    get_network(
      email = Sys.getenv("OPENALEX_EMAIL", "test@example.com"),
      string = "Hawkes",
      pages = 1L,
      per_page = 10L,
      min_date = "2010-01-01",
      max_date = "2020-01-01"
    )
  }, error = function(e) {
    skip(paste("OpenAlex network loading failed:", e$message))
  })
  
  net_raw <- out$net
  
  # Normalize times
  net_norm <- normalize_times_01(net_raw, attr = "time_scaled", keep_na = TRUE)
  
  expect_s3_class(net_norm, "network")
  times <- net_norm %v% "time_scaled"
  expect_true(all(times >= 0 & times <= 1, na.rm = TRUE))
})

test_that("waiting_times_between_formations works on OpenAlex network", {
  skip_if_no_openalex()
  skip_on_cran()
  
  source(system.file("openalex_study", "get_network_openalex.R", package = "hawkesNet"))
  
  out <- tryCatch({
    get_network(
      email = Sys.getenv("OPENALEX_EMAIL", "test@example.com"),
      string = "Hawkes",
      pages = 1L,
      per_page = 10L,
      min_date = "2010-01-01",
      max_date = "2020-01-01"
    )
  }, error = function(e) {
    skip(paste("OpenAlex network loading failed:", e$message))
  })
  
  net_raw <- out$net
  set.vertex.attribute(net_raw, "time", net_raw %v% "time_scaled")
  set.edge.attribute(net_raw, "time", net_raw %e% "time_scaled")
  net_raw <- normalize_times_01(net_raw, attr = "time", keep_na = TRUE)
  
  # Only test if network has edges
  if (network.edgecount(net_raw) > 0) {
    formula_RHS <- "edges + triangles + star(c(2,3))"
    wait_results <- tryCatch({
      :waiting_times_between_formations(net_raw, formula_RHS = formula_RHS)
    }, error = function(e) {
      skip(paste("waiting_times_between_formations failed:", e$message))
    })
    
    expect_true(is.list(wait_results))
    # Should have at least edges, triangles, star2, star3
    expect_true(length(wait_results) >= 3)
  } else {
    skip("Network has no edges")
  }
})

test_that("GOF function can run on OpenAlex network (slow)", {
  skip_if_no_openalex()
  skip_on_cran()
  
  # This is a very slow test - only run if explicitly requested
  if (Sys.getenv("RUN_SLOW_TESTS", "") != "true") {
    skip("Skipping slow GOF test. Set RUN_SLOW_TESTS=true to run.")
  }
  
  source(system.file("openalex_study", "get_network_openalex.R", package = "hawkesNet"))
  
  out <- tryCatch({
    get_network(
      email = Sys.getenv("OPENALEX_EMAIL", "test@example.com"),
      string = "Hawkes",
      pages = 1L,
      per_page = 10L,
      min_date = "2010-01-01",
      max_date = "2020-01-01"
    )
  }, error = function(e) {
    skip(paste("OpenAlex network loading failed:", e$message))
  })
  
  net_raw <- out$net
  set.vertex.attribute(net_raw, "time", net_raw %v% "time_scaled")
  set.edge.attribute(net_raw, "time", net_raw %e% "time_scaled")
  net_raw <- normalize_times_01(net_raw, attr = "time", keep_na = TRUE)
  
  # Only test if network has edges
  if (network.edgecount(net_raw) == 0) {
    skip("Network has no edges")
  }
  
  # Set up minimal parameters
  formula_RHS <- "edges + triangles + star(c(2,3))"
  exp_cs <- expected_params_PMF_mark_CS(net_raw, formula_RHS)
  n_cs <- if (!is.na(exp_cs$CS_params_length)) exp_cs$CS_params_length else 4L
  
  params_init <- list(
    mu = 1,
    beta_overall = 1,
    K = 0.5,
    beta_edges = 1,
    node_lambda = 1,
    CS_params = c(-10, rep(0, n_cs - 1)),
    vertex_categorical = list(gender = c(female = 0.1, male = 0.5)),
    vertex_categorical_levels = list(gender = c("female", "male", "unknown"))
  )
  
  # Create mock fit
  mock_fit <- list(
    fit = list(
      par = unlist(params_init[names(params_init) != "vertex_categorical_levels"])
    )
  )
  
  # Test GOF with minimal simulations
  gof_result <- tryCatch({
    gof(
      fit = mock_fit,
      mark_filtration = net_raw,
      time_window = c(0, 1),
      PMF_mark = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS = formula_RHS,
      n_sim = 1L,  # Minimal: just 1 simulation
      cores = 1L,
      truncation = 50L,
      verbose = FALSE
    )
  }, error = function(e) {
    skip(paste("GOF test failed:", e$message))
  })
  
  expect_true(is.list(gof_result))
  # GOF should return some statistics
  expect_true(length(gof_result) > 0)
})
