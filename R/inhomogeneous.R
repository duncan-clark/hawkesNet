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

  out <- list(
    result = result0,
    func   = func,
    lambda = mu_at_t,
    kernel_sum = params$K * sum(decays0),
    decays = decays0,
    diffs  = diffs
  )
  if (!is.null(tmp$combined_inputs)) out$combined_inputs <- tmp$combined_inputs
  out
}


#' Build one combined intensity closure from per-event CS ingredients
#'
#' Used when \code{combine_intensity = TRUE}: stacks change_stats into one matrix,
#' one matrix multiply per evaluation (vectorized); per-event steps (plogis, decay,
#' in_mark, node_dens, kernel) stay in a loop. Faster than N separate multiplies.
#' @param combined_inputs_list List of length N of combined_inputs from PMF_mark_CS (each has change_stats, in_mark, diffs, ...).
#' @param diffs_kernel_list List of length N of kernel diffs (t_i - times[times < t_i]).
#' @param mu_vec Length-N background rate at each event time.
#' @param times Length-N event times (used only for length check).
#' @return List of one function \code{f(params)} returning numeric vector of length N, or NULL on error.
#' @noRd
build_combined_intensity_funcs <- function(combined_inputs_list, diffs_kernel_list, mu_vec, times) {
  N <- length(times)
  if (N == 0L || length(combined_inputs_list) != N || length(diffs_kernel_list) != N || length(mu_vec) != N) {
    return(NULL)
  }
  n_cs <- NA_integer_
  for (i in seq_len(N)) {
    inp <- combined_inputs_list[[i]]
    if (!is.null(inp) && !isTRUE(inp$degenerate_edges) && nrow(inp$change_stats) > 0L) {
      n_cs <- ncol(inp$change_stats)
      break
    }
  }
  if (!is.finite(n_cs) || n_cs < 1L) return(NULL)
  nrows <- vapply(combined_inputs_list, function(x) {
    if (is.null(x) || isTRUE(x$degenerate_edges)) 0L else nrow(x$change_stats)
  }, 0L)
  row_start <- c(0L, cumsum(nrows))
  change_stats_stacked <- do.call(rbind, lapply(seq_len(N), function(i) {
    inp <- combined_inputs_list[[i]]
    if (nrows[i] == 0L) matrix(0, 0, n_cs) else inp$change_stats
  }))
  total_rows <- row_start[N + 1L]
  if (total_rows == 0L) return(NULL)
  # Free duplicate change_stats from per-event lists (stacked copy is used instead)
  for (i in seq_len(N)) combined_inputs_list[[i]]$change_stats <- NULL
  eps <- 1e-10
  expand_probs <- tryCatch(
    getFromNamespace("expand_vertex_categorical_probs", "hawkesNet"),
    error = function(e) NULL
  )
  if (is.null(expand_probs)) return(NULL)
  combined_closure <- function(params) {
    eta_all <- as.vector(change_stats_stacked %*% params$CS_params)
    out <- numeric(N)
    for (i in seq_len(N)) {
      inp <- combined_inputs_list[[i]]
      if (is.null(inp) || isTRUE(inp$degenerate_edges) || nrows[i] == 0L) {
        # Degenerate: mark density = 1, so intensity = mu + K * sum(decays)
        diffs_k <- diffs_kernel_list[[i]]
        decays <- if (length(diffs_k) > 0L) exp(-params$beta_overall * diffs_k) else numeric(0)
        out[i] <- max(mu_vec[i] + params$K * sum(decays), 1e-10)
        next
      }
      idx_i <- (row_start[i] + 1L):row_start[i + 1L]
      eta_i <- eta_all[idx_i]
      p_base <- stats::plogis(eta_i)
      p <- p_base * exp(-params$beta_edges * inp$diffs)
      p[p > (1 - eps)] <- 1 - eps
      p[p < eps] <- eps
      if (anyNA(p)) { out[i] <- 1e-10; next }
      log_edge_part <- sum(log(p[inp$in_mark])) + sum(log1p(-p[!inp$in_mark]))
      node_dens <- if (!is.null(inp$max_node_time) && inp$time > inp$max_node_time) {
        0
      } else {
        dval <- stats::dpois(inp$new_nodes - inp$old_nodes, params$node_lambda)
        if (!is.finite(dval) || dval <= 0) -1e10 else log(dval)
      }
      vcat <- params$vertex_categorical
      if (!is.null(vcat) && is.list(vcat) && length(inp$observed_categorical) > 0) {
        for (attr_name in names(inp$observed_categorical)) {
          if (attr_name %in% names(vcat)) {
            levs_attr <- inp$level_names_by_attr[[attr_name]]
            if (is.null(levs_attr) && !is.null(params$vertex_categorical_levels) && attr_name %in% names(params$vertex_categorical_levels))
              levs_attr <- params$vertex_categorical_levels[[attr_name]]
            if (!is.null(levs_attr)) {
              p_attr <- expand_probs(vcat[[attr_name]], levs_attr, eps = eps)
              if (!is.null(p_attr)) {
                obs <- inp$observed_categorical[[attr_name]]
                idx <- match(obs, names(p_attr))
                idx[is.na(idx)] <- match("unknown", names(p_attr))
                idx[is.na(idx)] <- 1L
                node_dens <- node_dens + sum(log(pmax(p_attr[idx], eps)))
              }
            }
          }
        }
      }
      log_dens <- log_edge_part + node_dens
      diffs_k <- diffs_kernel_list[[i]]
      decays <- if (length(diffs_k) > 0L) exp(-params$beta_overall * diffs_k) else numeric(0)
      intensity <- exp(log_dens) * (mu_vec[i] + params$K * sum(decays))
      out[i] <- max(intensity, 1e-10)
    }
    out
  }
  list(combined_closure)
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
    times_precalc <- list(times = times)
    dot_args <- list(...)
    formula_rhs <- dot_args$formula_RHS
    cores <- dot_args$cores
    use_parallel <- !is.null(cores) && is.numeric(cores) && cores > 1 && requireNamespace("pbmcapply", quietly = TRUE)
    shared_model <- NULL
    if (!use_parallel && !is.null(formula_rhs)) {
      g0 <- network::network.initialize(0L, directed = FALSE)
      if ("na" %in% network::list.vertex.attributes(g0)) network::delete.vertex.attribute(g0, "na")
      shared_model <- createCppModel(as.formula(paste("g0 ~ ", formula_rhs)))
      shared_model$setNetwork(ernm::as.BinaryNet(g0))
    }
    combine_intensity <- isTRUE(dot_args$combine_intensity)
    if (use_parallel) message("Intensity cache: using ", cores, " cores (pbmclapply)")
    else message("Intensity cache: using 1 core")
    if (combine_intensity) message("Intensity cache: will combine CS closures into one (saves closure envs)")
    intens_func <- function(i) {
      current_net <- filtration_to_net(mark_filtration, times[i], equals = TRUE)
      model <- if (!is.null(shared_model)) shared_model else if (!is.null(formula_rhs)) {
        createCppModel(as.formula(paste("current_net ~ ", formula_rhs)))
      } else NULL
      intensity <- cond_intensity_inhom(
        new_net = current_net,
        t = times[i],
        mark_filtration = current_net,
        PMF_mark = PMF_mark,
        params = params,
        mu_at_t = mu_vec[i],
        model = model,
        times = times_precalc,
        ...,
        return_combined_inputs = combine_intensity
      )
      out <- list(result = intensity$result, func = intensity$func)
      if (combine_intensity) {
        out$combined_inputs <- intensity$combined_inputs
        out$diffs_kernel <- intensity$diffs
      }
      out
    }
    if (use_parallel) {
      intens_list <- pbmcapply::pbmclapply(
        rev(seq_along(times)), intens_func,
        mc.cores = cores, mc.preschedule = FALSE
      )
      intens_list <- rev(intens_list)
    } else {
      intens_list <- lapply(seq_along(times), intens_func)
    }
    intens_vec <- sapply(intens_list, function(x) x$result)
    intens_funcs <- lapply(intens_list, function(x) x$func)
    # If combine_intensity and all events returned combined_inputs, build one closure (same data, one env).
    if (combine_intensity) {
      ci <- lapply(intens_list, function(x) x$combined_inputs)
      dk <- lapply(intens_list, function(x) x$diffs_kernel)
      if (all(vapply(ci, function(x) !is.null(x), NA))) {
        combined <- build_combined_intensity_funcs(
          combined_inputs_list = ci,
          diffs_kernel_list = dk,
          mu_vec = mu_vec,
          times = times
        )
        if (!is.null(combined)) {
          intens_funcs <- combined
          message("Intensity cache: combined ", length(times), " closures into 1 (memory: one list of change_stats, no per-event envs)")
        }
      }
    }
  } else {
    # Use cached closures. Do NOT parallelize here: each optim iteration would fork (e.g. 100)
    # processes; with many closures and large envs, fork overhead dominates and optimization
    # can be 10–50x slower than sequential (especially for nodeMatch with many events).
    # Parallelism is used only for the one-off intensity cache build above.
    if (length(intens_funcs) == 1L) {
      v <- intens_funcs[[1]](params)
      if (length(v) > 1L) intens_vec <- v else { intens_vec <- numeric(1); intens_vec[1] <- v }
    } else {
      intens_vec <- numeric(length(intens_funcs))
      for (i in seq_along(intens_funcs)) intens_vec[i] <- intens_funcs[[i]](params)
    }
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
#' @param ... Passed to PMF_mark (formula_RHS, truncation, cores, etc.) and to the loglik.
#'   \code{combine_intensity}: if \code{TRUE}, combine per-event CS intensity closures into one
#'   (saves closure envs; same total change_stats data in one list, so memory does not blow up).
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
  # Shallow copy so stripping levels does not modify params_init_old (needed for loglik and relist restore)
  params_init <- as.list(params_init_old)
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
  t_fit_start <- proc.time()[3]
  if (cache_intensity) {
    message("Pre-calculating intensity closures (inhomogeneous background)...")
    t_cache_start <- proc.time()[3]
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
    message("Intensity cache: ", round(proc.time()[3] - t_cache_start, 1), " s")
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
  n_par <- length(flat_par)
  times_for_diag <- get_times(mark_filtration)$times
  n_events_actual <- length(times_for_diag)
  if (length(cached_funcs) == 1L) {
    message("Optimizing ", n_par, " parameters, 1 combined closure (", n_events_actual, " events, vectorized eval)")
  } else {
    message("Optimizing ", n_par, " parameters, ", length(cached_funcs), " cached intensity closures (", n_events_actual, " events, sequential eval)")
  }
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
  t_optim_start <- proc.time()[3]
  fit <- do.call(optim, optim_args)
  t_optim_elapsed <- round(proc.time()[3] - t_optim_start, 1)
  n_iter <- if (!is.null(fit$counts)) fit$counts[1L] else NA_integer_
  n_fneval <- if (!is.null(fit$counts) && length(fit$counts) >= 2L) fit$counts[2L] else NA_integer_
  message("Optimization: ", t_optim_elapsed, " s | iterations: ", n_iter, if (is.finite(n_fneval)) paste0(" | fn evals: ", n_fneval) else "", " | s/iter: ", if (is.finite(n_iter) && n_iter > 0) round(t_optim_elapsed / n_iter, 2) else "n/a")
  message("Fitting (inhomogeneous) total: ", round(proc.time()[3] - t_fit_start, 1), " s")

  # Results table: estimate and standard error (from numerical Hessian)
  par_names <- names(fit$par)
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
    message("(Standard errors not available; Hessian inversion failed.)")
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
