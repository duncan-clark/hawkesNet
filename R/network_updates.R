#' Birth and turnover models with network-dependent Hawkes timing
#'
#' A compact reference implementation of the undirected anchor model in the
#' paper's recovery study. At each event Poisson births are added, an anchor is
#' selected, and its incident edges are added/deleted by independent Poisson
#' proposals conditional on at least one edit. PA uses log one plus target
#' degree; CS-1 uses signed changes in triangles, two-stars/6 and three-stars/18,
#' all evaluated before the complete update. There are no vertex deletions.
#'
#' Timing is `mu + sum(K * exp(-beta * lag) * amplitude)`, with amplitude
#' `exp(gamma_B * initial_summary + gamma_Delta * clipped_summary_change)`.
#' M0 fixes both coefficients to zero; M1 fixes `gamma_Delta` to zero. The
#' anchor threshold and degree exposure summaries are bounded. Log degree is
#' unbounded: simulation requires `gamma_B <= 0` and finite clipping to certify
#' non-explosion here. This conservative restriction does not limit evaluation
#' of a finite observed path's likelihood.
#'
#' Labels are sequential positive integers and persist throughout a path.
#' Matrices and candidate sets are explicit: this reference implementation is
#' intended for small networks, with quadratic timing storage for fitting.
#' The specialised `cran_*` interfaces provide the scalable CRAN implementation.
#' @param mark Update law, `"PA"` or `"CS-1"`.
#' @param timing Timing model, `"M0"`, `"M1"` or `"M2"`.
#' @param summary Anchor degree summary: threshold indicator, degree divided
#'   by `scale + degree`, or log one plus degree.
#' @param parameters Named numeric vector of model parameters, on their natural
#'   scales. Missing parameters use the defaults printed in `model$parameters`.
#'   Strictly positive parameters are `mu`, `K`, `beta`, `omega_B`, `omega_A`,
#'   `omega_D`; boundary submodels with zero rates are not fitted by this API.
#'   PA has `theta_A`, `theta_D`; CS-1 has `theta_A_T`, `theta_A_2`, `theta_A_3`
#'   and the analogous deletion coefficients. Only active timing coefficients
#'   appear in the vector.
#' @param threshold Positive integer degree threshold.
#' @param scale Positive degree-exposure scale.
#' @param clip Positive bound on later summary changes; `Inf` disables clipping.
#' @param repeat_probability Probability of reusing the previous anchor;
#'   otherwise an anchor is uniform over all vertices after births. Known, not fitted.
#' @return A model specification containing parameters and fixed design choices.
#' @examples
#' model <- network_update_model(timing = "M2")
#' set.seed(12)
#' path <- simulate_network_updates(model, horizon = 2)
#' loglik_network_updates(path, model)$loglik
#' @export
network_update_model <- function(mark = c("PA", "CS-1"),
                                 timing = c("M0", "M1", "M2"),
                                 summary = c("anchor_threshold", "degree_exposure", "log_degree"),
                                 parameters = NULL, threshold = 3L, scale = 1,
                                 clip = log(4), repeat_probability = 0.9) {
  mark <- match.arg(mark); timing <- match.arg(timing); summary <- match.arg(summary)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
  if (!scalar(threshold) || !is.finite(threshold) || threshold < 1 || threshold != floor(threshold))
    stop("threshold must be a positive integer")
  if (!scalar(scale) || !is.finite(scale) || scale <= 0 ||
      !scalar(clip) || clip <= 0) stop("scale must be finite positive; clip must be positive")
  if (!scalar(repeat_probability) || !is.finite(repeat_probability) ||
      repeat_probability < 0 || repeat_probability > 1) stop("repeat_probability must be in [0, 1]")
  p <- c(mu = 2, K = .25, beta = 1,
         if (timing != "M0") c(gamma_B = .7),
         if (timing == "M2") c(gamma_Delta = 1.2),
         omega_B = 1, omega_A = 2, omega_D = 2,
         if (mark == "PA") c(theta_A = .8, theta_D = .4) else
           c(theta_A_T = .6, theta_A_2 = .3, theta_A_3 = -.15,
             theta_D_T = .4, theta_D_2 = .2, theta_D_3 = -.1))
  if (!is.null(parameters)) {
    if (!is.numeric(parameters) || is.null(names(parameters)) || anyDuplicated(names(parameters)) ||
        any(!names(parameters) %in% names(p))) stop("parameters must have unique active parameter names")
    p[names(parameters)] <- parameters
  }
  positive <- p[c("mu", "K", "beta", "omega_B", "omega_A", "omega_D")]
  if (any(!is.finite(p)) || any(positive <= 0) || !is.finite(sum(positive)))
    stop("parameters must be finite, with positive baseline, kernel and proposal parameters")
  structure(list(mark = mark, timing = timing, summary = summary, parameters = p,
                 threshold = threshold, scale = scale, clip = clip,
                 repeat_probability = repeat_probability), class = "network_update_model")
}

.validate_update_model <- function(model) {
  if (!inherits(model, "network_update_model")) stop("model must come from network_update_model()")
  fields <- c("mark", "timing", "summary", "parameters", "threshold", "scale", "clip", "repeat_probability")
  if (!all(fields %in% names(model))) stop("incomplete network_update_model")
  checked <- do.call(network_update_model, model[fields])
  if (!setequal(names(checked$parameters), names(model$parameters)))
    stop("model parameters do not match its active specification")
  invisible(model)
}

.update_summary <- function(degree, model) switch(model$summary,
  anchor_threshold = as.numeric(degree >= model$threshold),
  degree_exposure = degree / (model$scale + degree), log_degree = log1p(degree))

.update_amplitudes <- function(initial, current, model) {
  p <- model$parameters
  if (model$timing == "M0") return(rep(1, length(initial)))
  predictor <- p["gamma_B"] * initial
  if (model$timing == "M2") predictor <- predictor + p["gamma_Delta"] *
    pmax(-model$clip, pmin(model$clip, current - initial))
  result <- exp(predictor)
  if (any(!is.finite(result))) stop("non-finite excitation amplitude")
  result
}

.update_initial <- function(initial) {
  initial <- as.matrix(initial)
  if (length(dim(initial)) != 2L || nrow(initial) != ncol(initial) || nrow(initial) < 2L ||
      anyNA(initial) || !all(initial %in% c(0, 1)) ||
      !isTRUE(all.equal(initial, t(initial), check.attributes = FALSE)) || any(diag(initial) != 0))
    stop("initial must be a symmetric 0/1 adjacency matrix with at least two vertices and no loops")
  storage.mode(initial) <- "logical"; initial
}

.update_births <- function(graph, births) {
  if (!births) return(graph)
  n <- nrow(graph); out <- matrix(FALSE, n + births, n + births)
  out[seq_len(n), seq_len(n)] <- graph; out
}

.update_candidates <- function(graph, anchor, mark) {
  targets <- setdiff(seq_len(nrow(graph)), anchor)
  degree <- rowSums(graph); present <- graph[anchor, targets]
  if (mark == "PA") x <- matrix(log1p(degree[targets]), ncol = 1L) else {
    du <- degree[anchor]; dv <- degree[targets]
    common <- as.numeric(graph[targets, , drop = FALSE] %*% graph[anchor, ])
    choose2 <- function(x) x * (x - 1) / 2
    x <- cbind(ifelse(present, -common, common),
               ifelse(present, -(du + dv - 2), du + dv) / 6,
               ifelse(present, -(choose2(du - 1) + choose2(dv - 1)),
                      choose2(du) + choose2(dv)) / 18)
  }
  list(targets = targets, present = present, x = x)
}

.update_rates <- function(candidates, model) {
  p <- model$parameters; log_rates <- numeric(length(candidates$targets))
  for (operation in c("A", "D")) {
    index <- which(candidates$present == (operation == "D"))
    if (!length(index)) next
    names <- if (model$mark == "PA") paste0("theta_", operation) else
      paste0("theta_", operation, c("_T", "_2", "_3"))
    weights <- as.numeric(candidates$x[index, , drop = FALSE] %*% p[names])
    if (any(!is.finite(weights))) stop("non-finite candidate score; reduce coefficient magnitudes")
    maximum <- max(weights)
    log_rates[index] <- log(p[paste0("omega_", operation)]) +
      weights - maximum - log(sum(exp(weights - maximum)))
  }
  log_rates
}

.update_mark_logprob <- function(candidates, model, births, anchor_logprob) {
  lr <- .update_rates(candidates, model); rates <- exp(lr)
  chosen <- candidates$selected
  # log(1-exp(-exp(lr))) stays finite when exp(lr) underflows.
  included <- lr[chosen]
  ordinary <- included > -30
  included[ordinary] <- log(-expm1(-exp(included[ordinary])))
  sum(included) - sum(rates[!chosen]) - log(-expm1(-sum(rates))) +
    stats::dpois(births, model$parameters["omega_B"], log = TRUE) + anchor_logprob
}

#' Record a path of complete birth-and-turnover updates
#' @param times Strictly increasing positive event times before `horizon`.
#' @param updates List of updates. Each is a list with `births` (nonnegative
#'   count), `anchor` (vertex label), and integer vectors `add` and `delete`
#'   of targets of edges incident to the anchor. Birth labels are the next
#'   consecutive integers. Every update must edit at least one edge.
#' @param initial Symmetric 0/1 adjacency matrix before the first event.
#' @param horizon Positive observation endpoint, including terminal survival.
#' @return A `network_update_data` object. Structural validity is checked on replay
#'   by the likelihood; the simulated object also has this format.
#' @export
network_update_data <- function(times, updates, initial, horizon) {
  initial <- .update_initial(initial)
  if (!is.numeric(horizon) || length(horizon) != 1L || !is.finite(horizon) || horizon <= 0 ||
      !is.numeric(times) || any(!is.finite(times)) || any(times <= 0 | times >= horizon) ||
      any(diff(times) <= 0) || !is.list(updates) || length(updates) != length(times))
    stop("provide one update per strictly increasing time in (0, horizon)")
  structure(list(times = times, updates = updates, initial = initial,
                 horizon = horizon), class = "network_update_data")
}

#' Simulate a birth-and-turnover network by exact thinning
#' @param model Model returned by `network_update_model()`.
#' @param horizon Positive simulation endpoint.
#' @param initial Symmetric initial adjacency matrix; defaults to one edge.
#' @param max_events,max_vertices Positive integer computation guards. Reaching
#'   a guard raises an error rather than returning a truncated path.
#' @return A path for `loglik_network_updates()` or `fit_network_updates()`.
#' @export
simulate_network_updates <- function(model, horizon,
                                     initial = matrix(c(0, 1, 1, 0), 2L),
                                     max_events = 100000L, max_vertices = 10000L) {
  .validate_update_model(model)
  network_update_data(numeric(), list(), initial, horizon)
  for (guard in list(max_events, max_vertices))
    if (length(guard) != 1L || !is.numeric(guard) || !is.finite(guard) ||
        guard < 1 || guard != floor(guard)) stop("computation guards must be positive integers")
  if (model$summary == "log_degree" && model$timing != "M0" &&
      (model$parameters["gamma_B"] > 0 || (model$timing == "M2" && !is.finite(model$clip))))
    stop("log-degree simulation requires gamma_B <= 0 and finite M2 clipping for non-explosion")
  graph <- .update_initial(initial)
  if (nrow(graph) > max_vertices) stop("initial graph exceeds max_vertices")
  p <- model$parameters; times <- born_summary <- numeric(); anchors <- integer(); updates <- list(); t <- 0
  while (t < horizon) {
    current <- .update_summary(rowSums(graph)[anchors], model)
    amplitudes <- .update_amplitudes(born_summary, current, model)
    excitation <- p["K"] * sum(amplitudes * exp(-p["beta"] * (t - times)))
    bound <- p["mu"] + excitation
    if (!is.finite(bound)) stop("non-finite thinning bound")
    wait <- stats::rexp(1, bound); candidate <- t + wait
    if (candidate >= horizon) break
    if (candidate <= t) stop("event time cannot advance at machine precision")
    t <- candidate
    if (stats::runif(1) > (p["mu"] + excitation * exp(-p["beta"] * wait)) / bound) next
    if (length(times) >= max_events) stop("max_events exceeded; no partial path returned")
    births <- stats::rpois(1, p["omega_B"])
    if (nrow(graph) + births > max_vertices) stop("max_vertices exceeded; no partial path returned")
    graph <- .update_births(graph, births)
    anchor <- if (length(anchors) && stats::runif(1) < model$repeat_probability)
      tail(anchors, 1) else sample.int(nrow(graph), 1L)
    candidates <- .update_candidates(graph, anchor, model$mark)
    chosen <- .proposal_indicators(exp(.update_rates(candidates, model)), TRUE)
    added <- candidates$targets[chosen & !candidates$present]
    deleted <- candidates$targets[chosen & candidates$present]
    graph[anchor, added] <- graph[added, anchor] <- TRUE
    graph[anchor, deleted] <- graph[deleted, anchor] <- FALSE
    times <- c(times, t); anchors <- c(anchors, anchor)
    born_summary <- c(born_summary, .update_summary(sum(graph[anchor, ]), model))
    updates[[length(times)]] <- list(births = births, anchor = anchor, add = added, delete = deleted)
  }
  network_update_data(times, updates, initial, horizon)
}

.prepare_updates <- function(data, model) {
  if (!inherits(data, "network_update_data")) stop("data must come from network_update_data()")
  network_update_data(data$times, data$updates, data$initial, data$horizon)
  graph <- data$initial; n <- length(data$times)
  anchors <- integer(n); births <- integer(n); initial <- anchor_lp <- numeric(n)
  features <- candidates <- vector("list", n)
  integer_vector <- function(x) is.numeric(x) && !anyNA(x) && all(is.finite(x)) && all(x == floor(x))
  for (i in seq_len(n)) {
    u <- data$updates[[i]]
    if (!is.list(u) || !all(c("births", "anchor", "add", "delete") %in% names(u)) ||
        !integer_vector(u$births) || length(u$births) != 1L || u$births < 0 ||
        !integer_vector(u$anchor) || length(u$anchor) != 1L ||
        !integer_vector(u$add) || !integer_vector(u$delete)) stop("malformed update at index ", i)
    graph <- .update_births(graph, u$births)
    if (u$anchor < 1 || u$anchor > nrow(graph)) stop("invalid anchor at update ", i)
    edits <- c(u$add, u$delete)
    if (!length(edits) || anyDuplicated(edits) || any(edits < 1 | edits > nrow(graph) | edits == u$anchor) ||
        any(graph[u$anchor, u$add]) || any(!graph[u$anchor, u$delete])) stop("inadmissible edits at update ", i)
    anchors[i] <- u$anchor; births[i] <- u$births
    probability <- if (i == 1L) 1 / nrow(graph) else
      (1 - model$repeat_probability) / nrow(graph) +
        model$repeat_probability * (u$anchor == anchors[i - 1L])
    anchor_lp[i] <- log(probability)
    candidates[[i]] <- .update_candidates(graph, u$anchor, model$mark)
    candidates[[i]]$selected <- candidates[[i]]$targets %in% edits
    graph[u$anchor, u$add] <- graph[u$add, u$anchor] <- TRUE
    graph[u$anchor, u$delete] <- graph[u$delete, u$anchor] <- FALSE
    features[[i]] <- .update_summary(rowSums(graph)[anchors[seq_len(i)]], model)
    initial[i] <- tail(features[[i]], 1L)
  }
  list(times = data$times, horizon = data$horizon, candidates = candidates, births = births,
       initial = initial, features = features, anchor_logprob = anchor_lp)
}

.evaluate_updates <- function(cache, model) {
  times <- cache$times; n <- length(times); p <- model$parameters
  rates <- rep(p["mu"], n); integrated <- p["mu"] * cache$horizon
  marks <- numeric(n)
  for (i in seq_len(n)) {
    history <- seq_len(i)
    amplitudes <- .update_amplitudes(cache$initial[history], cache$features[[i]], model)
    excitation <- p["K"] * sum(amplitudes * exp(-p["beta"] * (times[i] - times[history])))
    endpoint <- if (i == n) cache$horizon else times[i + 1L]
    integrated <- integrated + excitation * (-expm1(-p["beta"] * (endpoint - times[i]))) / p["beta"]
    if (i < n) rates[i + 1L] <- p["mu"] + excitation * exp(-p["beta"] * (endpoint - times[i]))
    marks[i] <- .update_mark_logprob(cache$candidates[[i]], model, cache$births[i], cache$anchor_logprob[i])
  }
  timing <- sum(log(rates)) - integrated
  list(loglik = as.numeric(timing + sum(marks)), timing = as.numeric(timing), marks = sum(marks),
       intensities = unname(rates), compensator = as.numeric(integrated), mark_contributions = marks)
}

#' Evaluate the exact likelihood of a birth-and-turnover path
#' @param data Path from `network_update_data()` or `simulate_network_updates()`.
#' @param model Model from `network_update_model()`.
#' @return List with total, timing and mark log likelihoods, pre-event
#'   intensities, compensator and per-update mark contributions.
#' @export
loglik_network_updates <- function(data, model) {
  .validate_update_model(model)
  .evaluate_updates(.prepare_updates(data, model), model)
}

#' Fit a birth-and-turnover network model
#'
#' Maximises the unpenalised complete likelihood. Positive parameters are
#' optimised on the log scale. Numerical bounds are returned; boundary fits
#' have no reported Wald standard errors. Multiple starting values can be
#' supplied as a list of named vectors. Convergence is local, not a guarantee
#' of a global maximum or identifiability. The model's design choices are fixed.
#' @param data Path from `network_update_data()`.
#' @param model Model whose parameters supply the default initial values.
#' @param starts Optional list of complete named parameter vectors on natural scales.
#' @param fixed Character vector of parameter names held at model values.
#' @param control Controls passed to `stats::optim()`.
#' @param hessian Compute local observed-information standard errors and
#'   approximate 95 percent Wald intervals, only for an interior optimum with
#'   positive definite, well-conditioned information. These are descriptive
#'   approximations, without a general asymptotic guarantee for growing networks.
#' @return List containing the fitted model, parameters, log likelihood,
#'   convergence, numerical bounds and all optimiser results.
#' @export
fit_network_updates <- function(data, model, starts = NULL, fixed = character(),
                                control = list(), hessian = TRUE) {
  .validate_update_model(model)
  cache <- .prepare_updates(data, model)
  if (!length(cache$times)) stop("an empty path cannot identify the model parameters")
  p <- model$parameters
  if (!is.character(fixed) || any(!fixed %in% names(p))) stop("fixed must name model parameters")
  positive <- names(p) %in% c("mu", "K", "beta", "omega_B", "omega_A", "omega_D")
  lower <- upper <- setNames(rep(3, length(p)), names(p)); lower[] <- -3
  lower[positive] <- log(c(mu = 1e-4, K = 1e-5, beta = .001, omega_B = .01, omega_A = .01, omega_D = .01)[names(p)[positive]])
  upper[positive] <- log(c(mu = 30, K = 15, beta = 15, omega_B = 10, omega_A = 20, omega_D = 20)[names(p)[positive]])
  if (model$mark == "PA") lower[c("theta_A", "theta_D")] <- -2
  encode <- function(p) { p[positive] <- log(p[positive]); p }
  decode <- function(x) { x[positive] <- exp(x[positive]); x }
  active <- which(!names(p) %in% fixed); base <- encode(p)
  if (is.null(starts)) starts <- list(p)
  if (!is.list(starts) || !length(starts)) stop("starts must be a nonempty list of parameter vectors")
  objective <- function(x) {
    full <- base; full[active] <- x; candidate <- model; candidate$parameters <- decode(full)
    value <- -.evaluate_updates(cache, candidate)$loglik
    if (!is.finite(value)) stop("non-finite likelihood; check data, model and starting values")
    value
  }
  fits <- lapply(starts, function(start) {
    if (!is.numeric(start) || anyDuplicated(names(start)) || !setequal(names(start), names(p)) || any(!is.finite(start)))
      stop("each start must name every active model parameter exactly once")
    start <- start[names(p)]
    if (any(start[positive] <= 0)) stop("positive parameters must have positive starts")
    x <- encode(start); x[names(p) %in% fixed] <- base[names(p) %in% fixed]
    if (any(x[active] < lower[active] | x[active] > upper[active])) stop("starting values fall outside numerical bounds")
    if (!length(active)) return(list(par = numeric(), value = objective(numeric()), convergence = 0L))
    stats::optim(x[active], objective, method = "L-BFGS-B", lower = lower[active], upper = upper[active],
                 control = utils::modifyList(list(maxit = 500), control))
  })
  best <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]
  base[active] <- best$par; model$parameters <- decode(base)
  boundary <- names(p)[active][abs(base[active] - lower[active]) < 1e-4 | abs(base[active] - upper[active]) < 1e-4]
  standard_errors <- setNames(rep(NA_real_, length(p)), names(p))
  intervals <- matrix(NA_real_, length(p), 2L, dimnames = list(names(p), c("lower", "upper")))
  if (isTRUE(hessian) && length(active) && !length(boundary) && best$convergence == 0L) {
    information <- stats::optimHess(base[active], objective)
    eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
    if (all(is.finite(eigenvalues)) && min(eigenvalues) > 1e-7 && max(eigenvalues) / min(eigenvalues) < 1e10) {
      transformed_se <- sqrt(diag(solve(information)))
      multiplier <- ifelse(positive[active], model$parameters[active], 1)
      standard_errors[active] <- transformed_se * multiplier
      intervals[active, ] <- cbind(model$parameters[active] - 1.96 * standard_errors[active],
                                    model$parameters[active] + 1.96 * standard_errors[active])
    }
  }
  list(model = model, parameters = model$parameters, loglik = -best$value,
       convergence = best$convergence, boundary = boundary,
       standard_errors = standard_errors, intervals = intervals,
       bounds = list(lower = decode(lower), upper = decode(upper)), fits = fits)
}

#' Check information for later-topology feedback at a fitted M1 model
#' @param data Observed path.
#' @param model Fitted M1 model from `fit_network_updates()$model`.
#' @param step Positive finite-difference step on natural parameter scales.
#' @return Score, efficient outer-product information, predicted standard error
#'   and score statistic for gamma_Delta=0. A rank-deficient nuisance information
#'   matrix returns NA for the last three, with an explanation. The score test
#'   assumes an interior M1 likelihood maximum and suitable asymptotic regularity.
#' @export
feedback_score <- function(data, model, step = 1e-5) {
  .validate_update_model(model)
  if (!inherits(model, "network_update_model") || model$timing != "M1") stop("supply a fitted M1 model")
  if (length(step) != 1L || !is.finite(step) || step <= 0) stop("step must be positive and finite")
  model$timing <- "M2"; model$parameters <- c(model$parameters, gamma_Delta = 0)
  cache <- .prepare_updates(data, model)
  names <- c("mu", "K", "beta", "gamma_B", "gamma_Delta")
  gradients <- matrix(0, length(cache$times), length(names)); score <- numeric(length(names))
  for (j in seq_along(names)) {
    name <- names[j]; h <- step * max(1, abs(model$parameters[name]))
    if (name %in% c("mu", "K", "beta")) h <- min(h, model$parameters[name] / 2)
    plus <- minus <- model; plus$parameters[name] <- plus$parameters[name] + h
    minus$parameters[name] <- minus$parameters[name] - h
    a <- .evaluate_updates(cache, plus); b <- .evaluate_updates(cache, minus)
    gradients[, j] <- (log(a$intensities) - log(b$intensities)) / (2 * h)
    score[j] <- (a$timing - b$timing) / (2 * h)
  }
  information <- crossprod(gradients); nuisance <- information[1:4, 1:4, drop = FALSE]
  if (qr(nuisance)$rank < 4L) return(list(score = score[5], information = NA_real_,
    predicted_se = NA_real_, statistic = NA_real_, reason = "rank-deficient M1 nuisance information"))
  efficient <- as.numeric(information[5, 5] - information[5, 1:4] %*% solve(nuisance, information[1:4, 5]))
  efficient <- max(0, efficient)
  list(score = score[5], information = efficient,
       predicted_se = if (efficient > 0) 1 / sqrt(efficient) else Inf,
       statistic = if (efficient > 0) score[5] / sqrt(efficient) else NA_real_)
}
