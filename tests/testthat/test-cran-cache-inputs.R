write_test_npy <- function(dtype, shape, bytes = raw(), version = 1L,
                           fortran = FALSE, header_length = NULL) {
  file <- tempfile(fileext = ".npy")
  header <- sprintf("{'descr': '%s', 'fortran_order': %s, 'shape': (%s), }\n",
                    dtype, if (fortran) "True" else "False", shape)
  con <- file(file, "wb")
  on.exit(close(con))
  writeBin(c(as.raw(0x93), charToRaw("NUMPY"), as.raw(c(version, 0))), con)
  writeBin(as.integer(if (is.null(header_length)) nchar(header, type = "bytes") else header_length),
           con, size = if (version == 1L) 2L else 4L, endian = "little")
  writeBin(charToRaw(header), con)
  writeBin(bytes, con)
  file
}

test_that("NPY import preserves numeric row order and packed uint64 bytes", {
  for (version in 1:3) {
    path <- write_test_npy("<f8", "2, 3", writeBin(as.double(1:6), raw(), size = 8, endian = "little"), version)
    expect_equal(cran_read_npy(path), matrix(as.double(1:6), 2, 3, byrow = TRUE))
    unlink(path)
  }
  path <- write_test_npy("<i4", "3,", writeBin(c(-2L, 0L, 3L), raw(), size = 4, endian = "little"))
  expect_equal(cran_read_npy(path), c(-2, 0, 3))
  unlink(path)
  bytes <- as.raw(c(255, 0, 17, 64, 128, 1, 200, 255))
  path <- write_test_npy("<u8", "1,", bytes)
  expect_identical(cran_read_npy(path), structure(bytes, npy_shape = 1L))
  unlink(path)
  path <- write_test_npy("<f8", "0, 3")
  expect_equal(dim(cran_read_npy(path)), c(0L, 3L))
  unlink(path)
})

test_that("NPY import rejects malformed versions, types, shapes and truncated data", {
  cases <- list(
    list(dtype = "<f1", shape = "1,", bytes = as.raw(0)),
    list(dtype = ">f8", shape = "1,", bytes = raw(8)),
    list(dtype = "<f8", shape = "1,", version = 4L),
    list(dtype = "<f8", shape = "1,", header_length = 9999L),
    list(dtype = "<f8", shape = "2,", bytes = raw(8)),
    list(dtype = "<f8", shape = "2,,3"),
    list(dtype = "<f8", shape = "2 3"),
    list(dtype = "<f8", shape = "-1,"),
    list(dtype = "<f8", shape = "2147483648,"),
    list(dtype = "<f8", shape = "1,", fortran = TRUE))
  for (args in cases) {
    path <- do.call(write_test_npy, args)
    expect_error(cran_read_npy(path), "NPY|Fortran")
    unlink(path)
  }
  path <- tempfile(); writeBin(as.raw(c(0x93, 0x4e)), path)
  expect_error(cran_read_npy(path), "Invalid NPY")
  unlink(path)
})

test_candidate_input <- function() {
  path <- tempfile(); writeBin(raw(), path)
  list(a = list(initial_degree = c(1, 1), initial_birth = c(0, 0),
    group_batch = 0, group_outdegree = 1,
    update_batch = numeric(), update_kind = numeric(), update_old_degree = numeric(),
    update_new_degree = numeric(), update_birth = numeric(),
    existing_degree = numeric(), existing_triangle = numeric(), existing_outdegree = numeric(),
    existing_birth = numeric(), existing_group = numeric(),
    group_source_degree = 1, group_source_birth = 0),
    index = list(offsets = c(0, 0), birth = numeric()), codes = raw(), rawpath = path,
    history = list(group = numeric(), features = matrix(numeric(), 0, 5), flags = matrix(numeric(), 0, 2)),
    time = 1, undirected = FALSE)
}

test_that("candidate cache validates integer arrays before native casts and indexing", {
  input <- test_candidate_input(); on.exit(unlink(input$rawpath))
  expect_type(do.call(cran_native_cache, input), "externalptr")
  historical <- input
  historical$undirected <- TRUE
  historical$history$features <- matrix(numeric(), 0, 7)
  expect_type(do.call(cran_native_cache, historical), "externalptr")
  historical$undirected <- FALSE
  expect_error(do.call(cran_native_cache, historical), "history feature dimensions")
  for (invalid in c(NA_real_, NaN, Inf, .5, 1e30)) {
    bad <- input; bad$index$offsets[2] <- invalid
    expect_error(do.call(cran_native_cache, bad), "finite|fractional|range")
  }
  bad <- input; bad$a$group_batch <- 1
  expect_error(do.call(cran_native_cache, bad), "batch index")
  terminal <- input
  terminal$a$update_batch <- 1
  terminal$a$update_kind <- 0
  terminal$a$update_old_degree <- 1
  terminal$a$update_new_degree <- 2
  terminal$a$update_birth <- 0
  expect_type(do.call(cran_native_cache, terminal), "externalptr")
  terminal$a$update_batch <- 2
  expect_error(do.call(cran_native_cache, terminal), "degree-update batch index")
  bad <- input; bad$time <- Inf
  expect_error(do.call(cran_native_cache, bad), "times must be finite")
  bad <- input
  for (name in c("existing_degree", "existing_triangle", "existing_outdegree", "existing_birth")) bad$a[[name]] <- 0
  bad$a$existing_group <- 1
  expect_error(do.call(cran_native_cache, bad), "existing-candidate group")
  expect_error(cran_group_sum(c(0, .5), matrix(1, 2, 1), 2L), "group-sum index")
  expect_error(cran_group_sum(NA_real_, matrix(1, 1, 1), 2L), "group-sum index")
})

test_release_input <- function() list(time = 1, new_count = 0,
  initial = matrix(c(0, 0, 0, 0, 0, 0), 1, 6),
  changes = matrix(numeric(), 0, 6), change_ptr = c(0, 0),
  selected = matrix(c(0, 1, 1), 1, 3), selected_ptr = c(0, 1))

test_that("direct and cached release evaluation reject malformed native inputs", {
  data <- test_release_input(); p <- c(rep(log(.1), 3), rep(0, 12))
  expect_true(is.finite(cran_release_evaluate(data, p, 1L)$loglik))
  mutations <- list(
    function(x) {x$initial <- matrix(0, 1, 5); x},
    function(x) {x$changes <- matrix(0, 1, 5); x},
    function(x) {x$change_ptr <- c(0, 1); x},
    function(x) {x$selected_ptr <- c(0, 2); x},
    function(x) {x$selected_ptr <- c(0, NA_real_); x},
    function(x) {x$selected_ptr <- c(0, .5); x},
    function(x) {x$selected[1, 1] <- 1; x},
    function(x) {x$selected[1, 2] <- .5; x},
    function(x) {x$selected[1, 3] <- 0; x},
    function(x) {x$new_count <- Inf; x},
    function(x) {x$new_count <- .5; x},
    function(x) {x$time <- NA_real_; x},
    function(x) {x$initial[1, 1] <- 1e30; x},
    function(x) {x$initial[1, 2] <- NaN; x})
  for (f in mutations) {
    bad <- f(data)
    expect_error(cran_release_evaluate(bad, p, 1L))
    expect_error(cran_release_mixture_cache(bad, p, 1L))
  }
  for (bad in list(p[-1], replace(p, 1, Inf), replace(p, 2, NA_real_))) {
    expect_error(cran_release_evaluate(data, bad, 1L), "parameters")
    expect_error(cran_release_mixture_cache(data, bad, 1L), "parameters")
  }
  expect_error(cran_release_evaluate(data, p, 0L), "family count")
  inactive <- data; inactive$selected[1, 1] <- 1
  expect_error(cran_release_evaluate(inactive, p, 2L), "not yet active")
})

test_that("cached release evaluation checks all dimensions and index arrays", {
  p <- c(rep(log(.1), 3), rep(0, 12))
  cache <- cran_release_mixture_cache(test_release_input(), p, 1L)
  expect_equal(cran_release_mixture_cached(cache, p, p)[[1]],
               cran_release_evaluate(test_release_input(), p, 1L))
  mutations <- list(
    function(x) {x$selection <- numeric(); x},
    function(x) {x$centered <- matrix(0, 1, 5); x},
    function(x) {x$fulln <- Inf; x},
    function(x) {x$upto <- 2; x},
    function(x) {x$batch <- .5; x},
    function(x) {x$batch <- 1L; x},
    function(x) {x$channel <- NA_integer_; x},
    function(x) {x$channel <- 2L; x},
    function(x) {x$logprob <- Inf; x},
    function(x) {x$centered[1, 1] <- NA_real_; x})
  for (f in mutations) expect_error(cran_release_mixture_cached(f(cache), p, p))
})
