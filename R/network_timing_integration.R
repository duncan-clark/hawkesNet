# Integration of network-dependent timing with the existing whole-update PMFs.

.feedback_cond_intensity <- function(tmp, t, params, mark_filtration, timing,
                                     mu_at_t = NULL) {
  cache <- .prepare_network_timing(mark_filtration, timing)
  ldf <- tmp$log_density_func
  if (!is.function(ldf)) stop("PMF_mark must return a callable log_density_func")
  force(t); force(mu_at_t)
  ground <- function(p) {
    if (!is.null(mu_at_t)) p$mu <- 0
    value <- hawkes_ground_intensity(t, p, NULL, timing, cache)
    if (!is.null(mu_at_t)) value <- value + mu_at_t
    value
  }
  log_func <- function(p) ldf(p) + log(ground(p))
  func <- function(p) exp(log_func(p))
  baseline <- if (is.null(mu_at_t)) params$mu else mu_at_t
  list(result = exp(tmp$log_mark_density + log(ground(params))), func = func,
       log_func = log_func, lambda = baseline,
       kernel_sum = ground(params) - baseline,
       diffs = t - cache$times[cache$times < t])
}

# Cache mark log densities separately from the parameter-independent graph path.
# In particular, changing gamma or beta during fitting reweights every interval.
.loglik_network_timing <- function(params, time_window, mark_filtration, PMF_mark,
                                   timing, mu_vec = NULL, integral_bg = NULL,
                                   intens_funcs = NULL, edge_hash_list = NULL, ...) {
  .network_timing_validate(timing)
  use_inhom <- !is.null(mu_vec) || !is.null(integral_bg)
  if (use_inhom && (is.null(mu_vec) || is.null(integral_bg)))
    stop("Supply both mu_vec and integral_bg for an inhomogeneous likelihood")
  if (length(time_window) != 2L || any(!is.finite(time_window)) ||
      time_window[2L] < time_window[1L]) stop("time_window must contain finite ordered endpoints")
  if (is.null(intens_funcs)) {
    cache <- .prepare_network_timing(mark_filtration, timing)
    times <- cache$times
    if (any(times < time_window[1L] | times > time_window[2L]))
      stop("realization has points outside time window")
    if (use_inhom && (length(mu_vec) != length(times) || any(!is.finite(mu_vec)) ||
                     any(mu_vec < 0) || length(integral_bg) != 1L ||
                     !is.finite(integral_bg) || integral_bg < 0))
      stop("mu_vec and integral_bg must be finite nonnegative background values matching event times")
    extra <- list(...)
    extra[c("cores", "combine_intensity", "parallel_type", "cache_intensity",
            "return_combined_inputs")] <- NULL
    mark_funcs <- lapply(seq_along(times), function(i) {
      current <- filtration_to_net(mark_filtration, times[i], equals = TRUE)
      args <- c(list(time = times[i], params = params, mark_filtration = current,
                     mark = current, generate_mark = FALSE,
                     new_edge_hash = if (is.null(edge_hash_list)) NULL else edge_hash_list[[i]]), extra)
      mark <- do.call(PMF_mark, args)
      if (!is.function(mark$log_density_func))
        stop("PMF_mark must return a callable log_density_func for feedback fitting")
      mark$log_density_func
    })
    # Force values retained by closures instead of retaining unevaluated calls.
    force(time_window); force(mu_vec); force(integral_bg); force(timing)
    log_values <- function(p) {
      ground_params <- p
      if (use_inhom) ground_params$mu <- 0
      lambda <- hawkes_ground_intensity(times, ground_params, NULL, timing, cache)
      if (use_inhom) lambda <- lambda + mu_vec
      marks <- vapply(mark_funcs, function(f) f(p), numeric(1))
      if (any(!is.finite(lambda)) || any(lambda <= 0)) return(rep(-Inf, length(times)))
      marks + log(lambda)
    }
    integral <- function(p) {
      if (use_inhom) p$mu <- 0
      hawkes_ground_compensator(p, time_window, NULL, timing, cache, integral_bg = integral_bg)
    }
    evaluate <- function(p) {
      values <- log_values(p)
      if (any(!is.finite(values))) return(-Inf)
      sum(values) - integral(p)
    }
    intens_funcs <- list(function(p) exp(log_values(p)))
    attr(intens_funcs, "timing") <- timing
    attr(intens_funcs, "feedback_loglik") <- evaluate
    attr(intens_funcs, "log_intensities") <- log_values
    attr(intens_funcs, "ground_compensator") <- integral
    attr(intens_funcs, "time_window") <- time_window
    attr(intens_funcs, "mu_vec") <- mu_vec
    attr(intens_funcs, "integral_bg") <- integral_bg
  } else if (!identical(attr(intens_funcs, "timing"), timing) ||
             !is.function(attr(intens_funcs, "feedback_loglik")) ||
             !identical(attr(intens_funcs, "time_window"), time_window) ||
             !identical(attr(intens_funcs, "mu_vec"), mu_vec) ||
             !identical(attr(intens_funcs, "integral_bg"), integral_bg)) {
    stop("intensity cache does not match the timing model, window, or background")
  }
  value <- attr(intens_funcs, "feedback_loglik")(params)
  list(loglik = value, intens_funcs = intens_funcs, timing = timing)
}

# The two built-in structural features remain constant between network events.
# Thus the current post-update rate bounds the excitation until the next event.
.simulate_network_timing <- function(params, time_window, PMF_mark, timing,
                                     seed_net = NULL, seed_times = NULL,
                                     inhom_bg = NULL, stop_on_full_network = TRUE,
                                     verbose = FALSE, ...) {
  started <- proc.time()[3L]
  net <- if (is.null(seed_net)) network::network.initialize(0, directed = FALSE) else
    network::network.copy(seed_net)
  cache <- .prepare_network_timing(net, timing)
  .network_timing_bound(params, timing)
  if (!is.null(seed_times) && (!is.numeric(seed_times) || any(!is.finite(seed_times)) ||
      !isTRUE(all.equal(sort(unique(seed_times)), cache$times, check.attributes = FALSE))))
    stop("M1/M2 seed_times must match the complete updates recorded in seed_net")
  t <- max(c(time_window[1L], cache$times))
  if (t > time_window[2L]) stop("seed history extends beyond the simulation window")
  has_background <- !is.null(inhom_bg)
  if (has_background) {
    mu_fun <- inhom_bg$mu_fit$mu_fun
    mu_bound <- inhom_bg$mu_fit$mu_bound
    if (is.null(mu_bound)) mu_bound <- inhom_bg$mu_bound
    if (!is.function(mu_fun) || length(mu_bound) != 1L || !is.finite(mu_bound) || mu_bound < 0)
      stop("M1/M2 inhomogeneous simulation requires mu_fit$mu_fun and a finite nonnegative mu_bound")
  } else {
    mu_fun <- function(t) params$mu
    mu_bound <- params$mu
  }
  events <- densities <- log_densities <- acceptance <- numeric(0)
  while (t < time_window[2L]) {
    n <- length(cache$times)
    features <- if (timing$model == "M1") cache$birth_features else
      if (n) cache$features_after[[n]] else numeric(0)
    weights <- .network_timing_alpha(features, params, timing)
    excitation <- params$K * sum(weights * exp(-params$beta_overall * (t - cache$times)))
    bound <- mu_bound + excitation
    if (!is.finite(bound)) stop("non-finite network timing thinning bound")
    if (bound == 0) break
    candidate <- t + stats::rexp(1L, bound)
    if (candidate > time_window[2L]) break
    if (candidate <= t) stop("Hawkes event spacing is below numerical precision")
    mu <- mu_fun(candidate)
    if (length(mu) != 1L || !is.finite(mu) || mu < 0 || mu > mu_bound)
      stop("background rate violates the declared mu_bound")
    rate <- mu + excitation * exp(-params$beta_overall * (candidate - t))
    probability <- rate / bound
    if (!is.finite(probability) || probability < 0 || probability > 1 + 1e-12)
      stop("network timing intensity exceeded its thinning bound")
    probability <- min(1, probability)
    acceptance <- c(acceptance, probability)
    t <- candidate
    if (stats::runif(1L) > probability) next
    sampled <- PMF_mark(time = t, params = params, mark_filtration = net, mark = NULL,
                        generate_mark = TRUE, generate_density = FALSE,
                        new_edge_hash = NULL, stop_on_full_network = stop_on_full_network, ...)
    if (is.null(sampled$mark_sample)) stop("mark sampler did not return mark_sample")
    log_q <- sampled$log_mark_sample_density
    if (is.null(log_q) && !is.null(sampled$mark_sample_density)) log_q <- log(sampled$mark_sample_density)
    if (length(log_q) != 1L || !is.finite(log_q) || log_q > 1e-8)
      stop("mark sampler returned an invalid generated log probability")
    net <- sampled$mark_sample
    next_cache <- .prepare_network_timing(net, timing)
    if (!isTRUE(all.equal(next_cache$times, c(cache$times, t), check.attributes = FALSE)))
      stop("mark sampler must record exactly one nonempty complete update at the accepted time")
    cache <- next_cache
    events <- c(events, t)
    log_densities <- c(log_densities, log_q)
    densities <- c(densities, exp(log_q))
  }
  if (verbose) message("Simulated ", length(events), " ", timing$model,
                       " events in ", round(proc.time()[3L] - started, 2), " seconds")
  list(events = list(n = length(events), t = events, mark_density = densities,
                     log_mark_density = log_densities), net = net,
       accept_probs = acceptance, timing = timing)
}
