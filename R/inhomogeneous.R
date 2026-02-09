# =============================================================================
# HawkesGrowthNet with inhomogeneous background rate (KDE)
# =============================================================================
# Replaces constant mu with time-varying mu(t) from KDE. Intensity:
# lambda(t) = mu(t) + K * sum(decays). Integral: int_0^T mu(s) ds + (1/beta)*K*sum(pieces).
# =============================================================================

#' Conditional intensity at one event time with inhomogeneous background
#'
#' Same as cond_intensity but uses mu_at_t (scalar) instead of params$mu.
#' @param new_net Network state at current time
#' @param t Current event time
#' @param mark_filtration Full mark filtration
#' @param PMF_mark Mark PMF function (e.g. PMF_mark_CS)
#' @param params Parameter list (beta_overall, K, CS_params, etc.; mu not used)
#' @param mu_at_t Scalar background rate at time t (e.g. from KDE)
#' @param new_edge_hash Optional hash of existing edges
#' @param times Precomputed times from get_times; if NULL, taken from mark_filtration
#' @param ... Passed to PMF_mark (e.g. formula_RHS, truncation)
#' @return List with result (intensity), func (function to evaluate intensity at new params), lambda, kernel_sum, decays, diffs
#' @export
cond_intensity_inhom <- function(new_net,
                                 t,
                                 mark_filtration = NULL,
                                 PMF_mark,
                                 params,
                                 mu_at_t,
                                 new_edge_hash = NULL,
                                 times = NULL,
                                 ...) {
  tmp <- PMF_mark(time = t,
                  params = params,
                  mark_filtration = mark_filtration,
                  mark = new_net,
                  generate_mark = FALSE,
                  new_edge_hash = new_edge_hash,
                  ...)
  if (is.null(times)) times <- get_times(mark_filtration)
  tt <- times$times
  tt <- tt[tt < t]
  diffs <- t - tt

  log_mark_density0 <- tmp$log_mark_density
  log_density_func  <- tmp$log_density_func

  decays0 <- exp(-params$beta_overall * diffs)
  log_result0 <- log_mark_density0 + log(mu_at_t + params$K * sum(decays0))
  result0 <- exp(log_result0)

  e_tiny <- new.env(parent = baseenv())
  e_tiny$diffs_local <- diffs
  e_tiny$ldf_local   <- log_density_func
  e_tiny$mu_at_t     <- mu_at_t

  func_template <- function(params) {
    decays <- exp(-params$beta_overall * diffs_local)
    exp(ldf_local(params) + log(mu_at_t + params$K * sum(decays)))
  }
  func <- func_template
  environment(func) <- e_tiny

  list(
    result = result0,
    func   = func,
    lambda = mu_at_t,
    kernel_sum = params$K * sum(decays0),
    decays = decays0,
    diffs  = diffs
  )
}


#' Log-likelihood for HawkesGrowthNet with inhomogeneous background
#'
#' mu_vec: length(times) vector of background rate at each event time.
#' integral_bg: scalar, integral of mu(s) from 0 to time_window[2].
#' @param params Parameter list (relist format)
#' @param time_window c(t0, t1)
#' @param mark_filtration Network with vertex/edge times
#' @param PMF_mark Mark PMF function
#' @param mu_vec Background rate at each event time (same length as times)
#' @param integral_bg Integral of background over time window
#' @param edge_hash_list Optional
#' @param verbose Print progress
#' @param intens_funcs Pre-computed intensity closures (for caching)
#' @param ... Passed to PMF_mark (formula_RHS, truncation, etc.)
#' @return List with loglik, intens_funcs
#' @noRd
loglik_hawkesNet_inhom <- function(params,
                                         time_window,
                                         mark_filtration,
                                         PMF_mark,
                                         mu_vec,
                                         integral_bg,
                                         edge_hash_list = NULL,
                                         verbose = FALSE,
                                         intens_funcs = NULL,
                                         ...) {
  times <- get_times(mark_filtration)
  times <- times$times

  tval <- time_window[2] - time_window[1]
  if (tval < max(times) - min(times)) {
    stop("realization has points outside time window")
  }

  if (length(mu_vec) != length(times)) {
    stop("mu_vec must have same length as event times")
  }

  if (is.null(intens_funcs)) {
    times_precalc <- get_times(mark_filtration)
    intens_func <- function(i) {
      current_net <- filtration_to_net(mark_filtration, times[i], equals = TRUE)
      if ("formula_RHS" %in% names(list(...))) {
        model <- createCppModel(as.formula(paste("current_net ~ ", list(...)$formula_RHS)))
      } else {
        model <- NULL
      }
      intensity <- cond_intensity_inhom(
        new_net = current_net,
        t = times[i],
        mark_filtration = current_net,
        PMF_mark = PMF_mark,
        params = params,
        mu_at_t = mu_vec[i],
        model = model,
        times = times_precalc,
        ...
      )
      list(result = intensity$result, func = intensity$func)
    }

    cores <- list(...)$cores
    if (!is.null(cores) && is.numeric(cores) && cores > 1 && requireNamespace("pbmcapply", quietly = TRUE)) {
      if (verbose) message("Using ", cores, " cores for intensity list")
      intens_list <- pbmcapply::pbmclapply(
        seq_along(times), intens_func,
        mc.cores = cores, mc.preschedule = FALSE
      )
    } else {
      intens_list <- lapply(seq_along(times), intens_func)
    }
    intens_vec <- sapply(intens_list, function(x) x$result)
    intens_funcs <- lapply(intens_list, function(x) x$func)
  } else {
    intens_vec <- numeric(length(intens_funcs))
    for (i in seq_along(intens_funcs)) intens_vec[i] <- intens_funcs[[i]](params)
  }

  tmp <- intens_vec
  if (any(is.na(tmp))) tmp[is.na(tmp)] <- min(tmp[!is.na(tmp)], na.rm = TRUE) / 2
  if (any(tmp <= 0)) tmp[tmp <= 0] <- min(tmp[tmp > 0]) / 2
  intens_sum <- sum(log(tmp))

  # Guard against division by zero or very small beta_overall
  if (!is.finite(params$beta_overall) || params$beta_overall <= 0 || params$beta_overall < 1e-10) {
    return(list(loglik = -1e10, intens_funcs = intens_funcs))
  }
  
  pieces <- 1 - exp(-params$beta_overall * (tval - times))
  integral <- integral_bg + (1 / params$beta_overall) * params$K * sum(pieces)
  loglik <- intens_sum - integral

  # Guard: ensure loglik is finite for L-BFGS-B
  if (!is.finite(loglik) || !is.finite(intens_sum) || !is.finite(integral)) {
    if (verbose) {
      warning("loglik_hawkesNet_inhom: non-finite loglik detected. intens_sum=", intens_sum,
              ", integral=", integral, ", returning -1e10")
    }
    return(list(loglik = -1e10, intens_funcs = intens_funcs))
  }

  list(
    loglik = loglik,
    intens_funcs = intens_funcs
  )
}


#' Fit HawkesGrowthNet with inhomogeneous (KDE) background
#'
#' mu_vec and integral_bg come from prepare_inhomogeneous_background (which uses estimate_mu_kde).
#' params_init can include mu but it is not used in the model.
#'
#' @param params_init List of initial parameters (same shape as for fit_hawkesNet)
#' @param time_window c(t0, t1)
#' @param mark_filtration Network with vertex/edge times
#' @param PMF_mark Mark PMF (e.g. PMF_mark_CS)
#' @param mu_vec Background rate at each event time (from prepare_inhomogeneous_background)
#' @param integral_bg Integral of background over time window
#' @param maxit Maximum iterations for optimizer
#' @param trace Trace level
#' @param reltol Relative tolerance
#' @param parscale Parameter scaling vector
#' @param fixed_params Names of parameters to fix
#' @param cache_intensity Pre-compute intensity closures for speed
#' @param method Optimization method: \code{"Nelder-Mead"} (default) or \code{"L-BFGS-B"}
#'   (gradient-based with box constraints; constrains mu, beta_overall, beta_edges,
#'   K, node_lambda > 0 automatically).
#' @param ... Passed to PMF_mark (formula_RHS, truncation, cores, etc.)
#' @return List with fit (optim result), intens_funcs, params_init_old, fit_table (parameter estimates and standard errors), and hessian (numerical Hessian of negative log-likelihood at MLE, if numDeriv available).
#' @export
fit_hawkesNet_inhom <- function(params_init,
                                      time_window,
                                      mark_filtration,
                                      PMF_mark,
                                      mu_vec,
                                      integral_bg,
                                      maxit,
                                      trace = 0,
                                      reltol = 1e-8,
                                      parscale = NULL,
                                      fixed_params = NULL,
                                      cache_intensity = TRUE,
                                      method = "Nelder-Mead",
                                      ...) {
  params_init_old <- params_init
  if (!is.null(fixed_params)) {
    for (k in fixed_params) params_init[[k]] <- NULL
  }
  # Strip vertex_categorical_levels (character metadata, not numeric parameters)
  # so unlist() yields a purely numeric vector for optim.
  params_init$vertex_categorical_levels <- NULL

  if (is.null(parscale)) {
    flat <- unlist(params_init)
    parscale <- rep(1, length(flat))
  }

  if (exists("validate_params_for_PMF", mode = "function")) {
    validate_params_for_PMF(params_init_old, PMF_mark, mark_filtration, ...)
  }

  cached_funcs <- NULL
  if (cache_intensity) {
    message("Pre-calculating intensity closures (inhomogeneous background)...")
    init_lik <- loglik_hawkesNet_inhom(
      params = params_init_old,
      time_window = time_window,
      mark_filtration = mark_filtration,
      PMF_mark = PMF_mark,
      mu_vec = mu_vec,
      integral_bg = integral_bg,
      ...
    )
    cached_funcs <- init_lik$intens_funcs
  }

  dot_args <- list(...)
  optim_func <- function(params) {
    params_curr <- relist(params, skeleton = params_init)
    # Restore vertex_categorical_levels from the original params (not optimized)
    params_curr$vertex_categorical_levels <- params_init_old$vertex_categorical_levels
    if (!is.null(fixed_params)) {
      for (k in fixed_params) params_curr[[k]] <- params_init_old[[k]]
    }
    if (!point_process_params_valid(params_curr)) return(-1e10)
    result <- tryCatch({
      do.call(loglik_hawkesNet_inhom, c(
        list(params = params_curr,
             time_window = time_window,
             mark_filtration = mark_filtration,
             PMF_mark = PMF_mark,
             mu_vec = mu_vec,
             integral_bg = integral_bg,
             intens_funcs = cached_funcs),
        dot_args
      ))
    }, error = function(e) {
      return(list(loglik = -1e10, intens_funcs = cached_funcs))
    })
    ll <- result$loglik
    # Final guard: ensure return value is finite for L-BFGS-B
    if (!is.finite(ll)) return(-1e10)
    ll
  }

  flat_par <- unlist(params_init)
  optim_args <- list(
    par = flat_par,
    fn = optim_func,
    method = method,
    control = list(fnscale = -1, trace = trace, maxit = maxit, reltol = reltol, parscale = parscale),
    hessian = TRUE
  )
  if (method == "L-BFGS-B") {
    bounds <- build_optim_bounds(names(flat_par))
    optim_args$lower <- bounds$lower
    optim_args$upper <- bounds$upper
  }
  fit <- do.call(optim, optim_args)

  message("Fitting (inhomogeneous) took ", round(proc.time()[3], 2), " seconds")

  # Results table: estimate and standard error (from numerical Hessian)
  par_names <- names(fit$par)
  # Debug: check parameter count
  if (length(par_names) != length(fit$par)) {
    warning("fit$par has ", length(fit$par), " elements but ", length(par_names), " names")
  }
  fit_table <- data.frame(
    parameter = par_names,
    estimate  = fit$par,
    std.error = NA_real_,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  hessian <- fit$hessian
  if (!is.null(hessian)) {
    # optim returns hessian of fn (loglik), which is negative definite at a maximum.
    # The variance-covariance matrix is the inverse of the *negative* hessian (observed information).
    vcov <- tryCatch(solve(-hessian), error = function(e) NULL)
    if (!is.null(vcov)) {
      se <- sqrt(pmax(diag(vcov), 0))
      fit_table$std.error <- se
    }
  }
  # Replace CS_params1, CS_params2, ... with actual ERNM statistic names
  fit_table <- rename_CS_params_in_table(fit_table, mark_filtration, list(...))
  cat("Inhomogeneous fit results:\n")
  print(fit_table, max = NULL, row.names = TRUE)
  if (all(is.na(fit_table$std.error))) {
    message("(Standard errors not available; install numDeriv for SEs.)")
  }

  list(
    fit = fit,
    intens_funcs = cached_funcs,
    params_init_old = params_init_old,
    fit_table = fit_table,
    hessian = hessian
  )
}


#' Build mu_vec and integral_bg from mark_filtration (network) using KDE
#'
#' Returns mu_vec aligned with get_times(mark_filtration)$times for use in
#' loglik_hawkesNet_inhom and fit_hawkesNet_inhom.
#'
#' @param mark_filtration Network with vertex and edge time attributes
#' @param time_attr Name of the time attribute (default "time")
#' @param bw Bandwidth for KDE; NULL = default
#' @param grid_n Number of grid points for KDE
#' @return List with mu_vec, integral_bg, times, mu_fit, Lambda_fun
#' @export
prepare_inhomogeneous_background <- function(mark_filtration, time_attr = "time", bw = NULL, grid_n = 2048) {
  times_obj <- get_times(mark_filtration, time_name = time_attr)
  times <- times_obj$times
  times <- times[!is.na(times)]
  if (length(times) < 5) stop("Need more event times for KDE")

  windowT <- c(min(times), max(times))
  mu_fit <- estimate_mu_kde(times, windowT = windowT, bw = bw, grid_n = grid_n)
  Lambda_obj <- make_cumhaz_fun(mu_fit)
  mu_vec_at_events <- mu_fit$mu_fun(times)
  mu_vec_at_events <- pmax(mu_vec_at_events, 1e-12)
  integral_bg <- Lambda_obj$Lambda_fun(windowT[2])

  list(
    mu_vec = mu_vec_at_events,
    integral_bg = integral_bg,
    times = times,
    mu_fit = mu_fit,
    Lambda_fun = Lambda_obj$Lambda_fun
  )
}
