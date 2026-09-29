#' Joint marked relational Hawkes model
#'
#' Exact pilot implementation of a marked Hawkes process on directed dyads where
#' past edges, values, and times jointly affect:
#' \enumerate{
#'   \item future event intensity via mark-dependent excitation weights, and
#'   \item future edge / value marks via a shared decayed relational state.
#' }
#'
#' Factorization for an event at time \eqn{t} with mark \eqn{(i,j,v)}:
#' \deqn{\lambda(t\mid H)\, q(i,j\mid t,H)\, f(v\mid i,j,t,H)}
#' where
#' \deqn{\lambda(t\mid H)=\mu + \sum_{t_m < t} K\, w_m\, e^{-\beta (t-t_m)},}
#' \eqn{w_m=\exp(\gamma^\top \phi_m)} depends on the previous event's value and
#' dyad features, \eqn{q} is a softmax over directed dyads from decayed
#' edge/recip/activity/value statistics, and \eqn{f} is a Gaussian density for
#' \eqn{\log(1+amount)}.
#'
#' This is the **joint** target model. The package's existing RHEM + HawkesNet
#' path is the **factorized** ablation (timing independent of mark values).
#'
#' @name joint_marked_hawkes
NULL

.jmw_as_hits <- function(x) {
    if (is.null(x) || nrow(x) == 0) {
        return(data.frame(
            t = numeric(0), i = integer(0), j = integer(0),
            weight = numeric(0), amount = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    if (!inherits(x, "data.frame")) {
        stop("Hits must be a data.frame with columns t, i, j.")
    }
    need <- c("t", "i", "j")
    if (!all(need %in% names(x))) {
        stop("Hits must include columns t, i, j.")
    }
    out <- x[order(x$t, seq_len(nrow(x))), , drop = FALSE]
    if (!("weight" %in% names(out))) {
        out$weight <- if ("amount" %in% names(out)) log1p(pmax(out$amount, 0)) else 1
    }
    if (!("amount" %in% names(out))) {
        out$amount <- pmax(expm1(out$weight), 0)
    }
    out$log_amount <- log1p(pmax(out$amount, 0))
    rownames(out) <- NULL
    out
}

.jmw_risk_set <- function(actors) {
    risk <- expand.grid(i = actors, j = actors,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    risk[risk$i != risk$j, , drop = FALSE]
}

#' Default parameters for the joint marked relational Hawkes model
#'
#' @param n_events Optional seed for \code{mu} from an observed count and window.
#' @param duration Observation window length used with \code{n_events}.
#' @return Named parameter list.
#' @export
joint_marked_hawkes_default_params <- function(n_events = 100, duration = 1) {
    list(
        mu = max(n_events / max(duration, 1e-8), 1e-3),
        K = 0.2,
        beta_overall = 2,
        beta_edges = 1,
        gamma_value = 0.3,
        gamma_edge = 0.1,
        RHEM_params = c(
            intercept = 0,
            repetition = 0.5,
            reciprocity = 0.3,
            sender_activity = 0.1,
            receiver_activity = 0.1,
            edge_value = 0.4
        ),
        amount_intercept = 0,
        amount_edge_value = 0.5,
        amount_repetition = 0.2,
        amount_sd = 0.5
    )
}

#' Event-driven decayed relational state for joint marked Hawkes
#'
#' Maintains decayed directed edge weights, values, and activity using exact
#' exponential decay between events. Exact complete-risk-set calculations are
#' intended for 100-1000 actor pilots.
#'
#' @param actors Integer or character actor IDs.
#' @param beta_edges Edge / activity decay rate.
#' @export
JointMarkedState <- function(actors, beta_edges = 1) {
    actors <- sort(unique(actors))
    n <- length(actors)
    id <- stats::setNames(seq_len(n), as.character(actors))
    list(
        actors = actors,
        id = id,
        n = n,
        beta_edges = beta_edges,
        time = 0,
        edge_w = matrix(0, n, n),
        edge_v = matrix(0, n, n),
        send_act = numeric(n),
        recv_act = numeric(n)
    )
}

.jmw_decay_state <- function(state, time) {
    dt <- time - state$time
    if (dt < 0) {
        stop("JointMarkedState times must be nondecreasing.")
    }
    if (dt > 0 && state$beta_edges > 0) {
        factor <- exp(-state$beta_edges * dt)
        state$edge_w <- state$edge_w * factor
        state$edge_v <- state$edge_v * factor
        state$send_act <- state$send_act * factor
        state$recv_act <- state$recv_act * factor
    }
    state$time <- time
    state
}

.jmw_add_event <- function(state, i, j, weight = 1, amount = NULL) {
    ii <- state$id[[as.character(i)]]
    jj <- state$id[[as.character(j)]]
    if (is.null(ii) || is.null(jj) || ii == jj) {
        return(state)
    }
    w <- if (is.finite(weight)) weight else 1
    v <- if (!is.null(amount) && is.finite(amount)) log1p(max(amount, 0)) else w
    state$edge_w[ii, jj] <- state$edge_w[ii, jj] + w
    state$edge_v[ii, jj] <- state$edge_v[ii, jj] + v
    state$send_act[ii] <- state$send_act[ii] + w
    state$recv_act[jj] <- state$recv_act[jj] + w
    state
}

.jmw_dyad_stats <- function(state, i, j) {
    ii <- state$id[[as.character(i)]]
    jj <- state$id[[as.character(j)]]
    c(
        intercept = 1,
        repetition = state$edge_w[ii, jj],
        reciprocity = state$edge_w[jj, ii],
        sender_activity = state$send_act[ii],
        receiver_activity = state$recv_act[jj],
        edge_value = state$edge_v[ii, jj]
    )
}

.jmw_all_dyad_stats <- function(state) {
    actors <- state$actors
    risk <- .jmw_risk_set(actors)
    stats <- matrix(0, nrow = nrow(risk), ncol = 6)
    colnames(stats) <- c("intercept", "repetition", "reciprocity",
                         "sender_activity", "receiver_activity", "edge_value")
    for (row in seq_len(nrow(risk))) {
        ii <- state$id[[as.character(risk$i[row])]]
        jj <- state$id[[as.character(risk$j[row])]]
        stats[row, ] <- c(
            1,
            state$edge_w[ii, jj],
            state$edge_w[jj, ii],
            state$send_act[ii],
            state$recv_act[jj],
            state$edge_v[ii, jj]
        )
    }
    list(stats = stats, risk = risk)
}

.jmw_softmax <- function(eta) {
    eta <- eta - max(eta)
    p <- exp(eta)
    p / sum(p)
}

.jmw_mark_excitation_weight <- function(params, stats_prev, log_amount_prev) {
    # Mark-dependent excitation: values and edge strength scale future timing.
    gamma_value <- params$gamma_value %||% 0
    gamma_edge <- params$gamma_edge %||% 0
    rep_stat <- if (!is.null(stats_prev)) stats_prev[["repetition"]] %||% 0 else 0
    exp(gamma_value * log_amount_prev + gamma_edge * rep_stat)
}

.jmw_amount_mean <- function(params, stats) {
    (params$amount_intercept %||% 0) +
        (params$amount_edge_value %||% 0) * (stats[["edge_value"]] %||% 0) +
        (params$amount_repetition %||% 0) * (stats[["repetition"]] %||% 0)
}

.jmw_amount_logdens <- function(params, log_amount, stats) {
    sd <- max(params$amount_sd %||% 0.5, 1e-6)
    mu <- .jmw_amount_mean(params, stats)
    stats::dnorm(log_amount, mean = mu, sd = sd, log = TRUE)
}

#' Build left-continuous joint state just before each event
#'
#' @param hits Event data frame.
#' @param actors Actor IDs.
#' @param beta_edges Decay rate.
#' @return List of states of length \code{nrow(hits)} (state before event k).
#' @keywords internal
.jmw_state_path <- function(hits, actors, beta_edges) {
    hits <- .jmw_as_hits(hits)
    state <- JointMarkedState(actors, beta_edges = beta_edges)
    path <- vector("list", nrow(hits))
    for (k in seq_len(nrow(hits))) {
        state <- .jmw_decay_state(state, hits$t[k])
        path[[k]] <- state
        state <- .jmw_add_event(state, hits$i[k], hits$j[k],
                                weight = hits$weight[k], amount = hits$amount[k])
    }
    path
}

#' Log-likelihood for the joint marked relational Hawkes model (exact dyad set)
#'
#' Uses a complete directed risk set. Compensator is approximated by a
#' left-continuous piecewise-constant ground intensity between events
#' (exact for the discrete-event marked likelihood contribution at event times;
#' suitable as the pilot reference).
#'
#' @param params Parameter list from \code{\link{joint_marked_hawkes_default_params}}.
#' @param hits Event data frame with \code{t}, \code{i}, \code{j}, and amount/weight.
#' @param actors Actor universe.
#' @param time_window Numeric length-2 observation window.
#' @return List with \code{loglik} and component sums.
#' @export
loglik_joint_marked_hawkes <- function(params,
                                       hits,
                                       actors = NULL,
                                       time_window = NULL) {
    hits <- .jmw_as_hits(hits)
    if (is.null(actors)) {
        actors <- sort(unique(c(hits$i, hits$j)))
    }
    if (is.null(time_window)) {
        time_window <- c(0, max(hits$t, 1))
    }
    T1 <- time_window[1]
    T2 <- time_window[2]
    if (nrow(hits) == 0) {
        return(list(loglik = -params$mu * (T2 - T1),
                    mark_loglik = 0, amount_loglik = 0,
                    ground_loglik = -params$mu * (T2 - T1)))
    }

    theta <- params$RHEM_params
    feat <- names(theta)
    if (is.null(feat) || any(!nzchar(feat))) {
        feat <- c("intercept", "repetition", "reciprocity",
                  "sender_activity", "receiver_activity", "edge_value")
        names(theta) <- feat[seq_along(theta)]
    }

    path <- .jmw_state_path(hits, actors, params$beta_edges %||% 1)
    mark_ll <- 0
    amount_ll <- 0
    ground_ll <- 0
    mu <- params$mu %||% 1e-3
    K <- params$K %||% 0
    beta <- params$beta_overall %||% 1

    # Excitation weights for each event (mark-dependent), used for ground process
    w_evt <- numeric(nrow(hits))
    for (k in seq_len(nrow(hits))) {
        state <- path[[k]]
        all <- .jmw_all_dyad_stats(state)
        eta <- as.vector(all$stats[, feat, drop = FALSE] %*% as.numeric(theta[feat]))
        probs <- .jmw_softmax(eta)
        idx <- which(all$risk$i == hits$i[k] & all$risk$j == hits$j[k])
        if (length(idx) != 1) {
            stop("Observed dyad not in risk set.")
        }
        mark_ll <- mark_ll + log(probs[idx])
        stats_obs <- .jmw_dyad_stats(state, hits$i[k], hits$j[k])
        amount_ll <- amount_ll + .jmw_amount_logdens(params, hits$log_amount[k], stats_obs)
        w_evt[k] <- .jmw_mark_excitation_weight(params, stats_obs, hits$log_amount[k])
    }

    # Ground loglik: log intensity at events minus integral of decaying kernel
    .jmw_integral <- function(t0, t1, kernel_at_t0) {
        dt <- t1 - t0
        if (dt <= 0) return(0)
        if (beta > 1e-12) {
            mu * dt + K * kernel_at_t0 * (1 - exp(-beta * dt)) / beta
        } else {
            (mu + K * kernel_at_t0) * dt
        }
    }

    t_cursor <- T1
    kernel_cursor <- 0
    for (k in seq_len(nrow(hits))) {
        t <- hits$t[k]
        ground_ll <- ground_ll - .jmw_integral(t_cursor, t, kernel_cursor)
        if (t > t_cursor && beta > 1e-12) {
            kernel_cursor <- kernel_cursor * exp(-beta * (t - t_cursor))
        }
        lam <- mu + K * kernel_cursor
        if (!(is.finite(lam) && lam > 0)) {
            return(list(loglik = -1e10, mark_loglik = mark_ll,
                        amount_loglik = amount_ll, ground_loglik = -1e10))
        }
        ground_ll <- ground_ll + log(lam)
        kernel_cursor <- kernel_cursor + w_evt[k]
        t_cursor <- t
    }
    ground_ll <- ground_ll - .jmw_integral(t_cursor, T2, kernel_cursor)

    list(
        loglik = mark_ll + amount_ll + ground_ll,
        mark_loglik = mark_ll,
        amount_loglik = amount_ll,
        ground_loglik = ground_ll
    )
}

#' Fit joint marked relational Hawkes by numerical MLE
#'
#' @inheritParams loglik_joint_marked_hawkes
#' @param params_init Initial parameters.
#' @param fixed_params Names of parameters to hold fixed.
#' @param maxit Optimizer iterations.
#' @param trace Optim trace.
#' @return Fit object with \code{fit}, \code{params}, and component logliks.
#' @export
fit_joint_marked_hawkes <- function(params_init,
                                    hits,
                                    actors = NULL,
                                    time_window = NULL,
                                    fixed_params = NULL,
                                    maxit = 50,
                                    trace = 0) {
    hits <- .jmw_as_hits(hits)
    if (is.null(actors)) {
        actors <- sort(unique(c(hits$i, hits$j)))
    }
    if (is.null(time_window)) {
        time_window <- c(0, max(1, max(hits$t)))
    }
    params_full <- params_init
    free <- params_init
    if (!is.null(fixed_params)) {
        free[fixed_params] <- NULL
    }
    par0 <- unlist(free)
    fn <- function(par) {
        cur <- utils::relist(par, skeleton = free)
        params <- params_full
        params[names(cur)] <- cur
        # Keep nested RHEM_params names
        if (!is.null(cur$RHEM_params)) {
            params$RHEM_params <- cur$RHEM_params
            if (!is.null(names(params_full$RHEM_params))) {
                names(params$RHEM_params) <- names(params_full$RHEM_params)
            }
        }
        if (!is.null(params$mu) && params$mu <= 0) return(-1e10)
        if (!is.null(params$K) && params$K < 0) return(-1e10)
        if (!is.null(params$beta_overall) && params$beta_overall <= 0) return(-1e10)
        if (!is.null(params$beta_edges) && params$beta_edges < 0) return(-1e10)
        if (!is.null(params$amount_sd) && params$amount_sd <= 0) return(-1e10)
        loglik_joint_marked_hawkes(params, hits, actors, time_window)$loglik
    }
    fit <- stats::optim(
        par = par0,
        fn = fn,
        method = "Nelder-Mead",
        control = list(fnscale = -1, maxit = maxit, trace = trace, reltol = 1e-8)
    )
    est <- utils::relist(fit$par, skeleton = free)
    params <- params_full
    params[names(est)] <- est
    if (!is.null(est$RHEM_params) && !is.null(names(params_full$RHEM_params))) {
        names(params$RHEM_params) <- names(params_full$RHEM_params)
    }
    ll <- loglik_joint_marked_hawkes(params, hits, actors, time_window)
    list(
        fit = fit,
        params = params,
        params_init = params_init,
        fixed_params = fixed_params,
        actors = actors,
        time_window = time_window,
        loglik = ll
    )
}

#' Score events under a joint marked Hawkes fit
#'
#' Returns decomposed surprise scores: ground, mark, amount, and total.
#'
#' @param hits Event data frame.
#' @param fit Object from \code{\link{fit_joint_marked_hawkes}}, or NULL.
#' @param params Parameter list if \code{fit} is NULL.
#' @param actors Actor universe.
#' @param score_idx Indices to score (default all).
#' @return Data frame of per-event scores.
#' @export
score_joint_marked_hawkes <- function(hits,
                                      fit = NULL,
                                      params = NULL,
                                      actors = NULL,
                                      score_idx = NULL) {
    hits <- .jmw_as_hits(hits)
    if (is.null(params)) {
        if (is.null(fit)) stop("Provide fit or params.")
        params <- fit$params
    }
    if (is.null(actors)) {
        actors <- if (!is.null(fit)) fit$actors else sort(unique(c(hits$i, hits$j)))
    }
    if (is.null(score_idx)) score_idx <- seq_len(nrow(hits))

    theta <- params$RHEM_params
    feat <- names(theta)
    path <- .jmw_state_path(hits, actors, params$beta_edges %||% 1)
    kernel_sum <- 0
    rows <- vector("list", length(score_idx))
    pos <- 1L
    prev_stats <- NULL
    prev_log_amount <- 0

    for (k in seq_len(nrow(hits))) {
        state <- path[[k]]
        if (k > 1) {
            dt <- hits$t[k] - hits$t[k - 1]
            kernel_sum <- kernel_sum * exp(-(params$beta_overall %||% 1) * dt)
            kernel_sum <- kernel_sum + .jmw_mark_excitation_weight(
                params, prev_stats, prev_log_amount
            )
        }
        lam <- (params$mu %||% 1e-3) + (params$K %||% 0) * kernel_sum
        all <- .jmw_all_dyad_stats(state)
        eta <- as.vector(all$stats[, feat, drop = FALSE] %*% as.numeric(theta[feat]))
        probs <- .jmw_softmax(eta)
        idx <- which(all$risk$i == hits$i[k] & all$risk$j == hits$j[k])
        log_mark <- log(probs[idx])
        stats_obs <- .jmw_dyad_stats(state, hits$i[k], hits$j[k])
        log_amount <- .jmw_amount_logdens(params, hits$log_amount[k], stats_obs)
        log_ground <- log(lam)
        if (k %in% score_idx) {
            rows[[pos]] <- data.frame(
                event_index = k,
                t = hits$t[k],
                i = hits$i[k],
                j = hits$j[k],
                log_ground = log_ground,
                log_mark = log_mark,
                log_amount = log_amount,
                score_ground = -log_ground,
                score_mark = -log_mark,
                score_amount = -log_amount,
                score_total = -(log_ground + log_mark + log_amount),
                true_rank = rank(-probs, ties.method = "average")[idx],
                true_prob = probs[idx],
                stringsAsFactors = FALSE
            )
            pos <- pos + 1L
        }
        prev_stats <- stats_obs
        prev_log_amount <- hits$log_amount[k]
    }
    do.call(rbind, rows)
}

#' Simulate a small joint marked relational Hawkes process
#'
#' Exact thinning over the complete dyad risk set for synthetic recovery tests.
#'
#' @param params Parameter list.
#' @param time_window Observation window.
#' @param actors Actor IDs.
#' @param max_events Safety cap.
#' @return List with \code{hits} and \code{actors}.
#' @export
sim_joint_marked_hawkes <- function(params,
                                    time_window = c(0, 1),
                                    actors = 1:8,
                                    max_events = 500) {
    actors <- sort(unique(actors))
    state <- JointMarkedState(actors, beta_edges = params$beta_edges %||% 1)
    t <- time_window[1]
    T2 <- time_window[2]
    theta <- params$RHEM_params
    feat <- names(theta)
    kernel_sum <- 0
    t_last_event <- t
    prev_stats <- NULL
    prev_log_amount <- 0
    rows <- list()
    n_events <- 0L

    while (t < T2 && n_events < max_events) {
        # Upper bound: decay state to now for mark bound; use current kernel
        mu <- params$mu %||% 1e-3
        K <- params$K %||% 0
        beta <- params$beta_overall %||% 1
        # Conservative bound: current intensity with no further decay of marks
        state_now <- .jmw_decay_state(state, t)
        all <- .jmw_all_dyad_stats(state_now)
        eta <- as.vector(all$stats[, feat, drop = FALSE] %*% as.numeric(theta[feat]))
        # Bound total ground by mu + K*kernel_sum (kernel already at time t)
        lam_bar <- mu + K * max(kernel_sum, 0) + 1e-8
        if (!is.finite(lam_bar) || lam_bar <= 0) {
            break
        }
        w <- stats::rexp(1, rate = lam_bar)
        t <- t + w
        if (!is.finite(t) || t >= T2) break
        # Decay kernel and state to candidate time
        dt <- t - t_last_event
        if (n_events > 0 && is.finite(dt) && dt > 0) {
            kernel_sum <- kernel_sum * exp(-beta * dt)
        }
        if (!is.finite(kernel_sum) || kernel_sum < 0) kernel_sum <- 0
        state <- .jmw_decay_state(state, t)
        lam <- mu + K * kernel_sum
        if (!is.finite(lam) || lam <= 0) {
            t_last_event <- t
            next
        }
        if (stats::runif(1) > min(1, lam / lam_bar)) {
            t_last_event <- t
            next
        }
        all <- .jmw_all_dyad_stats(state)
        eta <- as.vector(all$stats[, feat, drop = FALSE] %*% as.numeric(theta[feat]))
        probs <- .jmw_softmax(eta)
        pick <- sample.int(nrow(all$risk), 1, prob = probs)
        i <- all$risk$i[pick]
        j <- all$risk$j[pick]
        stats_obs <- .jmw_dyad_stats(state, i, j)
        sd <- max(params$amount_sd %||% 0.5, 1e-6)
        log_amount <- stats::rnorm(1, mean = .jmw_amount_mean(params, stats_obs), sd = sd)
        log_amount <- max(log_amount, 1e-6)
        amount <- expm1(log_amount)
        weight <- log_amount
        n_events <- n_events + 1L
        rows[[n_events]] <- data.frame(
            t = t, i = i, j = j, weight = weight, amount = amount,
            stringsAsFactors = FALSE
        )
        # At event time+, include this event's excitation weight into kernel for next
        w_evt <- .jmw_mark_excitation_weight(params, stats_obs, log_amount)
        kernel_sum <- kernel_sum + w_evt
        state <- .jmw_add_event(state, i, j, weight = weight, amount = amount)
        prev_stats <- stats_obs
        prev_log_amount <- log_amount
        t_last_event <- t
    }

    hits <- if (length(rows)) do.call(rbind, rows) else
        data.frame(t = numeric(0), i = integer(0), j = integer(0),
                   weight = numeric(0), amount = numeric(0))
    list(hits = hits, actors = actors, params = params)
}

#' Benchmark exact joint marked Hawkes likelihood scaling
#'
#' Times \code{\link{loglik_joint_marked_hawkes}} across actor/event grid.
#'
#' @param actor_grid Integer vector of actor counts.
#' @param event_grid Integer vector of event counts.
#' @param reps Replications per cell.
#' @return Data frame with runtime and risk-set size.
#' @export
benchmark_joint_marked_hawkes <- function(actor_grid = c(10, 20, 40),
                                          event_grid = c(50, 100),
                                          reps = 1) {
    params <- joint_marked_hawkes_default_params()
    rows <- list()
    pos <- 1L
    for (n_act in actor_grid) {
        for (n_evt in event_grid) {
            for (r in seq_len(reps)) {
                set.seed(1000 + n_act + n_evt + r)
                actors <- seq_len(n_act)
                sim <- sim_joint_marked_hawkes(
                    params,
                    time_window = c(0, 1),
                    actors = actors,
                    max_events = n_evt
                )
                hits <- sim$hits
                if (nrow(hits) > n_evt) hits <- hits[seq_len(n_evt), , drop = FALSE]
                mem <- NA_real_
                t0 <- proc.time()[[3]]
                ll <- tryCatch(
                    loglik_joint_marked_hawkes(params, hits, actors, c(0, 1)),
                    error = function(e) NULL
                )
                elapsed <- proc.time()[[3]] - t0
                rows[[pos]] <- data.frame(
                    n_actors = n_act,
                    n_events = nrow(hits),
                    risk_dyads = n_act * (n_act - 1),
                    rep = r,
                    seconds = elapsed,
                    loglik = if (is.null(ll)) NA_real_ else ll$loglik,
                    stringsAsFactors = FALSE
                )
                pos <- pos + 1L
            }
        }
    }
    do.call(rbind, rows)
}
