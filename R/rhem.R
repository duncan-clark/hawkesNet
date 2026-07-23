#' Repeated edge-hit model likelihood
#'
#' Computes the conditional mark log-likelihood for repeated directed edge hits.
#' This is the mark-only counterpart to the RHEM mark component used inside
#' \code{\link{loglik_hawkesNet}}.
#'
#' @param params Named list with \code{RHEM_params} and optionally
#'   \code{beta_edges}.
#' @param mark_filtration Data frame with columns \code{t}, \code{i}, and
#'   \code{j}; an optional \code{weight} column is used for valued hits.
#' @param actors Actor IDs defining the directed risk set.
#' @param formula_RHS Optional ERNM valued change-statistic formula RHS.
#' @param weight_col Name of hit-weight column.
#' @param rhem_stats_cache Optional environment used by \code{PMF_mark_RHEM}.
#' @param rhem_streaming_threshold Risk-set size threshold for streaming ERNM
#'   log-probability calculations.
#' @param ... Additional arguments passed to \code{PMF_mark_RHEM}.
#' @return A list with \code{loglik}.
#' @rdname RHEM
#' @export
loglik_RHEM <- function(params,
                        mark_filtration,
                        actors = NULL,
                        formula_RHS = NULL,
                        weight_col = "weight",
                        rhem_stats_cache = NULL,
                        rhem_streaming_threshold = Inf,
                        ...) {
    hits <- .rhem_as_hits(mark_filtration) # nolint: object_usage_linter
    if (nrow(hits) == 0) {
        return(list(loglik = 0))
    }
    if (!is.null(params$beta_edges) && params$beta_edges < 0) {
        return(list(loglik = -(10 ^ 100)))
    }
    if (is.null(actors)) {
        actors <- sort(unique(c(hits$i, hits$j)))
    }
    hits <- hits[order(hits$t, seq_len(nrow(hits))), , drop = FALSE]

    log_terms <- vapply(seq_len(nrow(hits)), function(idx) {
        history <- hits[seq_len(idx), , drop = FALSE]
        pmf <- PMF_mark_RHEM( # nolint: object_usage_linter
            time = hits$t[idx],
            params = params,
            mark_filtration = history,
            actors = actors,
            formula_RHS = formula_RHS,
            weight_col = weight_col,
            rhem_stats_cache = rhem_stats_cache,
            rhem_streaming_threshold = rhem_streaming_threshold,
            ...
        )
        pmf$log_mark_density
    }, numeric(1))

    list(loglik = sum(log_terms))
}

#' @param params_init Initial parameter list for \code{\link{fit_RHEM}}.
#' @param trace Optimizer trace level.
#' @param REPORT Optimizer report frequency.
#' @param reltol Optimizer relative tolerance.
#' @param maxit Maximum optimizer iterations.
#' @param get_hessian Whether to compute a numerical Hessian.
#' @param fixed_params Character vector of parameter names to hold fixed.
#' @rdname RHEM
#' @export
fit_RHEM <- function(params_init,
                     mark_filtration,
                     actors = NULL,
                     formula_RHS = NULL,
                     trace = 0,
                     REPORT = 10,
                     reltol = 1e-8,
                     maxit = 25,
                     get_hessian = FALSE,
                     fixed_params = NULL,
                     ...) {
    params_init_old <- params_init
    if (!is.null(fixed_params)) {
        params_init[fixed_params] <- NULL
    }

    fn_wrapper <- function(par, ...) {
        params <- utils::relist(par, skeleton = params_init)
        params[fixed_params] <- params_init_old[fixed_params]
        loglik_RHEM(
            params = params,
            mark_filtration = mark_filtration,
            actors = actors,
            formula_RHS = formula_RHS,
            ...
        )$loglik
    }

    t <- proc.time()
    fit <- stats::optim(
        par = unlist(params_init),
        fn = fn_wrapper,
        method = "Nelder-Mead",
        control = list(fnscale = -1,
                       trace = trace,
                       REPORT = REPORT,
                       maxit = maxit,
                       reltol = reltol,
                       abstol = NULL),
        hessian = get_hessian,
        ...
    )
    message("RHEM fitting took ", round((proc.time() - t)[3], 2), " seconds")

    list(fit = fit,
         params_init = params_init_old,
         optimized_skeleton = params_init,
         fixed_params = fixed_params,
         formula_RHS = formula_RHS)
}

#' @param fit Fitted object from \code{\link{fit_RHEM}}, or \code{NULL} when
#'   passing fixed \code{params}.
#' @param params Fixed parameter list used when \code{fit = NULL}.
#' @param score_idx Integer indices of events to score. Defaults to all events.
#' @rdname RHEM
#' @export
score_RHEM <- function(mark_filtration,
                       fit = NULL,
                       params = NULL,
                       actors = NULL,
                       formula_RHS = NULL,
                       score_idx = NULL,
                       weight_col = "weight",
                       rhem_streaming_threshold = Inf,
                       ...) {
    hits <- .rhem_as_hits(mark_filtration) # nolint: object_usage_linter
    if (is.null(actors)) {
        actors <- sort(unique(c(hits$i, hits$j)))
    }
    if (is.null(score_idx)) {
        score_idx <- seq_len(nrow(hits))
    }
    if (is.null(params)) {
        if (is.null(fit)) {
            stop("Provide either `fit` or `params`.")
        }
        estimated <- utils::relist(fit$fit$par, skeleton = fit$optimized_skeleton)
        params <- fit$params_init
        params[names(estimated)] <- estimated
    }

    log_mark_density <- vapply(score_idx, function(idx) {
        history <- hits[seq_len(idx), , drop = FALSE]
        pmf <- PMF_mark_RHEM( # nolint: object_usage_linter
            time = hits$t[idx],
            params = params,
            mark_filtration = history,
            actors = actors,
            formula_RHS = formula_RHS %||% if (!is.null(fit)) fit$formula_RHS else NULL,
            weight_col = weight_col,
            rhem_streaming_threshold = rhem_streaming_threshold,
            ...
        )
        pmf$log_mark_density
    }, numeric(1))

    data.frame(
        event_index = score_idx,
        log_mark_density = log_mark_density,
        negative_log_mark = -log_mark_density,
        stringsAsFactors = FALSE
    )
}
