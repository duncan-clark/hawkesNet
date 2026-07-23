# Generated merge: JASA main mark PMFs + RHEM / timeNet experimental backend
# See inst/RHEM_TIMENET.md

#' Mark probability mass function for the network generation process
#'
#' Calculates the probability mass function (PMF) for marks (network structures) at a given time, supporting multiple network growth models (Barabasi-Albert, Change Statistic Hawkes, and BA-bipartite). 
#'
#' @references
#' Barabasi, A.-L. & Albert, R. (1999). Emergence of scaling in random networks. *Science*, 286, 509–512. \doi{10.1126/science.286.5439.509}
#' @param time Numeric. The time at which to evaluate the PMF.
#' @param params Named list. Model parameter values required for the chosen \code{type}:
#'   - For \code{type = "BA"}: \code{beta_edges} (numeric).
#'   - For \code{type = "CS"}: \code{beta_edges} (numeric), \code{node_lambda} (numeric), \code{CS_params} (numeric vector of coefficients for \code{formula_RHS}).
#'   - For \code{type = "BA-bip"}: similar to BA, plus any relevant bipartite parameters.
#' @param mark_filtration Network or compatible object. The network history (filtration) up to the current time.
#' @param type Character. Model type: one of \code{"BA"}, \code{"CS"}, or \code{"BA-bip"}. Default: \code{"BA"}.
#' @param mark Network or NULL. The current mark/network structure. If \code{NULL}, derived from \code{mark_filtration}.
#' @param generate_mark Logical. If \code{TRUE}, generates a new mark/sample; otherwise computes density for the supplied mark. Default: \code{FALSE}.
#' @param generate_density Logical. If \code{TRUE}, computes the density for the provided mark. Default: \code{TRUE}.
#' @param grad Logical. If \code{TRUE}, also computes gradients of the mark density. Default: \code{FALSE}.
#' @param new_edge_hash Hash or NULL. Optional hashed edge list for fast lookup. Default: \code{NULL}.
#' @param truncation Integer or NULL. For \code{"CS"} models, controls which edges are considered (e.g., \code{1} means only new-to-old). Default: \code{NULL}.
#' @param formula_RHS Character or formula. For \code{"CS"} models, specifies the right-hand-side for change statistics calculation. Default: \code{NULL}.
#' @param mark_decay Character or NULL. How edge decay is modeled, e.g. \code{"node_entrance"}, \code{"activity"}, etc. Default: \code{NULL}.
#' @param model Object or NULL. Preconstructed model object (for efficiency); if \code{NULL}, will be built internally. Default: \code{NULL}.
#' @param max_node_time Numeric. The last time at which a node can enter the network (CS model). Default: \code{10}.
#' @param ... Additional arguments, passed to model-specific PMF functions.
#' @return Named list containing:
#'   \item{mark_density}{Numeric. Density of the provided or generated mark.}
#'   \item{log_mark_density}{Numeric. Log-density of the mark.}
#'   \item{edge_probs}{Numeric vector. Probabilities for each possible edge.}
#'   \item{mark_grad}{Numeric vector. Gradient of the mark density (if \code{grad = TRUE}).}
#'   \item{decay_grad}{Numeric vector. Gradient with respect to decay (if \code{grad = TRUE}).}
#'   \item{mark_sample}{Network. The sampled mark/network object.}
#'   \item{mark_sample_density}{Numeric. Density of the sampled mark.}
#'   \item{log_mark_sample_density}{Numeric. Log-density of the sampled mark.}
#'
#' @details Computes the mark PMF, \eqn{q(m\vert t,\mathcal{H}_{t})} (see \code{\link{cond_intensity}}).
#' Currently three options: \code{type = "BA"}, \code{type = "CS"}, and \code{type = "BA-bip"}. 
#'
#' For \code{type = "BA"} the Barabasi Albert (BA) preferential attachment model is used where the mark distribution is defined as
#' \deqn{
#'   q(m \mid t, \mathcal{H}_t) =
#'   \prod_{i=1}^{N_{t-}} \left(p_i^{BA}\right)^{e_i} \cdot
#'   \left(1 - p_i^{BA}\right)^{1 - e_i}.
#' }
#' Here, the attachment probability, \eqn{p_i^{BA}}, is defined as
#' \deqn{
#'   p_i^{BA} = \frac{\delta_i}{\sum_{k=1}^{N} \delta_k}
#' }
#' where \eqn{\delta_{i}^{t} = \exp(\tau \cdot (t - t_i)) \cdot d_{i}^{t}}, 
#' and \eqn{d_{i}^{t}} is the sna::degree of node \eqn{i} just before time \eqn{t}.
#' 
#' For \code{type = "CS"} the change statistic (CS) model is used where the mark distribution is defined (similar to above) as
#' \deqn{
#'   q(m \mid t, \mathcal{H}_t) =
#'   \prod_{i=1}^{N_{t-}} \left(p_i^{CS}\right)^{e_i} \cdot
#'   \left(1 - p_i^{CS}\right)^{1 - e_i}.
#' }
#' where the attachment probability, \eqn{p_i^{CS}}, is defined as
#' \deqn{
#' p_i^{CS} = \left(\nu + \exp(\tau \cdot(t - t_i))\right)\cdot\frac{1}{1 + \exp(-\theta^{\top} \cdot C_{i,N_t})}
#' }
#' @examples
#' \dontrun{
#' if(interactive()){
#'  data(net, package = "hawkesNet")
#' time <- get_times(net)$times
#' ## BA
#' params_ba <-  list( beta_edges = 0.1)
#' mark_filtration <-  filtration_to_net(net,10)
#' pmf_ba <- PMF_mark(time[10],  params_ba, mark_filtration)
#' ## CS
#' devtools::install_github("duncan-clark/ernm", ref = "R_change_stats")
#' require(ernm)
#' params_cs <-  list(beta_edges = 0.1,node_lambda = 1,CS_params =  c(-10,0,0,0))
#' pmf_cs <- PMF_mark(time = time[10],  params = params_cs,
#' mark_filtration = mark_filtration, type = "CS",  truncation = 1,
#' formula_RHS = "edges + triangles + star(c(2,3))",
#' max_node_time = 1)
#'  }
#' }
#' @seealso
#' \code{\link{cond_intensity}}
#' \code{\link{PMF_mark_BA}}
#' \code{\link{PMF_mark_CS}}
#' \code{\link[network]{network}}, \code{\link[network]{add.vertices}}
#' \code{\link[ernm]{as.BinaryNet}}
#' @rdname PMF_mark
#' @export
PMF_mark <- function(time,
                     params,
                     mark_filtration,
                     type = c("BA", "CS", "RHEM"),
                     mark = NULL,
                     generate_mark = FALSE,
                     generate_density = TRUE,
                     grad = FALSE,
                     new_edge_hash = NULL,
                     truncation = NULL,
                     formula_RHS = NULL,
                     mark_decay = 'node_entrance',
                     model = NULL,
                     max_node_time = 10,
                     ...){
    type <- type[1]
    if (!(type %in% c("BA", "CS", "BA-bip", "RHEM"))) {
        stop("type can only be one of `BA` for Barabasi-Albert, `CS` for change statistic Hawkes, `BA-bip` for bipartite Barabasi-Albert, or `RHEM` for repeated edge-hit models.")
    }
    if(type == "BA"){
        pmf <- PMF_mark_BA(time, params, mark_filtration, mark,
                           generate_mark, generate_density, grad,
                           new_edge_hash, truncation, ...)
    }else{
        if(type == "BA-bip"){
            pmf <- PMF_mark_BA_bipartite(time, params,  mark_filtration,
                                      mark, generate_mark, generate_density,
                                      grad, new_edge_hash, truncation, ...)
        }else{
            if(type == "CS"){
                pmf <- PMF_mark_CS(time, params, mark_filtration,
                                   mark, generate_mark, generate_density,
                                   grad,  new_edge_hash, truncation, formula_RHS,
                                   mark_decay, model, max_node_time, ...)
            }else{
                if(type == "RHEM"){
                    pmf <- PMF_mark_RHEM(time, params, mark_filtration,
                                         mark, generate_mark, generate_density,
                                         grad, formula_RHS = formula_RHS, ...)
                }
            }
        }
    }
    return(pmf)
}


#' Repeated edge-hit mark probability mass function
#'
#' Computes a softmax mark PMF for recurrent directed edge hits. The filtration is
#' a data frame with columns \code{t}, \code{i}, and \code{j}; an optional
#' \code{weight} column contributes to decayed edge states. The mark distribution
#' follows the change-statistic paradigm: each candidate dyad receives statistics
#' computed from the left-continuous hit history and is scored by
#' \code{RHEM_params}.
#'
#' @inheritParams PMF_mark
#' @param actors Optional vector of actor IDs defining the directed risk set.
#' @param hit_features Character vector of supported features. Defaults to
#'   \code{names(params$RHEM_params)} when named, otherwise
#'   \code{c("intercept", "repetition", "reciprocity")}.
#' @param weight_col Name of the optional hit-weight column. Default:
#'   \code{"weight"}.
#' @return Named list with the same shape as \code{\link{PMF_mark}}.
#' @rdname PMF_mark_RHEM
#' @export
PMF_mark_RHEM <- function(time,
                          params,
                          mark_filtration,
                          mark = NULL,
                          generate_mark = FALSE,
                          generate_density = TRUE,
                          grad = FALSE,
                          formula_RHS = NULL,
                          actors = NULL,
                          hit_features = NULL,
                          weight_col = "weight",
                          hit_delta = 1,
                          rhem_stats_cache = NULL,
                          rhem_streaming_threshold = Inf,
                          backend = c("auto", "simple", "ernm", "timeNet"),
                          timeNet_state = NULL,
                          ...){
    hits <- .rhem_as_hits(mark_filtration)
    mark_hits <- if (is.null(mark)) hits else .rhem_as_hits(mark)
    current <- .rhem_current_mark(mark_hits, time)
    history <- hits[hits$t < time, , drop = FALSE]

    if (is.null(actors)) {
        actors <- sort(unique(c(hits$i, hits$j, mark_hits$i, mark_hits$j)))
    }
    actors <- sort(unique(actors))
    if (length(actors) < 2) {
        stop("PMF_mark_RHEM requires at least two actors in `actors` or the hit history.")
    }

    backend <- .rhem_select_backend(backend, formula_RHS)
    if (backend == "timeNet") {
        pmf <- .rhem_timeNet_pmf(time = time,
                                 params = params,
                                 history = history,
                                 current = current,
                                 actors = actors,
                                 formula_RHS = formula_RHS,
                                 generate_mark = generate_mark,
                                 generate_density = generate_density,
                                 grad = grad,
                                 weight_col = weight_col,
                                 hit_delta = hit_delta,
                                 rhem_streaming_threshold = rhem_streaming_threshold,
                                 rhem_stats_cache = rhem_stats_cache,
                                 timeNet_state = timeNet_state)
        return(pmf)
    }

    edge_state <- .rhem_edge_state(history, time, actors, params$beta_edges %||% 0, weight_col)

    risk <- .rhem_risk_set(actors)
    if (!is.null(formula_RHS)) {
        stats <- .rhem_ernm_change_stats(risk, edge_state, actors, formula_RHS,
                                         hit_delta, rhem_stats_cache)
        theta <- .rhem_theta(params, colnames(stats))
    } else {
        features <- .rhem_features(params, hit_features)
        theta <- .rhem_theta(params, features)
        stats <- .rhem_candidate_stats(risk, edge_state, features)
    }
    eta <- as.vector(stats %*% theta)
    probs <- .rhem_softmax(eta)

    out <- list(mark_density = 1,
                log_mark_density = 0,
                edge_probs = probs,
                mark_grad = NULL,
                decay_grad = 0)

    if (generate_density && !is.null(current)) {
        idx <- which(risk$i == current$i[1] & risk$j == current$j[1])
        if (length(idx) != 1) {
            stop("Observed repeated-hit mark is not in the RHEM risk set.")
        }
        log_density <- log(probs[idx])
        out$mark_density <- exp(log_density)
        out$log_mark_density <- log_density
        if (grad) {
            expected_stats <- colSums(stats * probs)
            out$mark_grad <- stats[idx, ] - expected_stats
        }
    }

    if (generate_mark) {
        idx <- sample.int(nrow(risk), size = 1, prob = probs)
        sample_hit <- data.frame(t = time,
                                 i = risk$i[idx],
                                 j = risk$j[idx])
        sample_hit[[weight_col]] <- 1
        mark_sample <- rbind(hits, sample_hit)
        out$mark_sample <- mark_sample
        out$mark_sample_density <- probs[idx]
        out$log_mark_sample_density <- log(probs[idx])
    } else {
        out$mark_sample <- NULL
        out$mark_sample_density <- NULL
        out$log_mark_sample_density <- NULL
    }

    out
}

.rhem_as_hits <- function(x) {
    if (is.null(x)) {
        return(data.frame(t = numeric(0), i = integer(0), j = integer(0)))
    }
    if (!inherits(x, "data.frame")) {
        stop("RHEM mark filtrations must be data frames with columns `t`, `i`, and `j`.")
    }
    if (!all(c("t", "i", "j") %in% names(x))) {
        stop("RHEM mark filtrations must include columns `t`, `i`, and `j`.")
    }
    x <- x[order(x$t), , drop = FALSE]
    rownames(x) <- NULL
    x
}

.rhem_current_mark <- function(hits, time) {
    current <- hits[hits$t == time, , drop = FALSE]
    if (nrow(current) == 0) {
        return(NULL)
    }
    if (nrow(current) > 1) {
        stop("PMF_mark_RHEM expects at most one observed hit at each event time.")
    }
    current
}

.rhem_features <- function(params, hit_features) {
    if (!is.null(hit_features)) {
        return(hit_features)
    }
    if (!is.null(names(params$RHEM_params)) && all(nzchar(names(params$RHEM_params)))) {
        return(names(params$RHEM_params))
    }
    c("intercept", "repetition", "reciprocity")
}

.rhem_theta <- function(params, features) {
    theta <- params$RHEM_params
    if (is.null(theta)) {
        stop("params$RHEM_params is required for PMF_mark_RHEM.")
    }
    if (!is.null(names(theta)) && all(features %in% names(theta))) {
        theta <- theta[features]
    }
    if (length(theta) != length(features)) {
        stop("Length of params$RHEM_params must match the selected hit features.")
    }
    as.numeric(theta)
}

.rhem_risk_set <- function(actors) {
    risk <- expand.grid(i = actors, j = actors,
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
    risk <- risk[risk$i != risk$j, , drop = FALSE]
    rownames(risk) <- NULL
    risk
}

.rhem_edge_state <- function(history, time, actors, beta_edges, weight_col) {
    state <- matrix(0, nrow = length(actors), ncol = length(actors),
                    dimnames = list(as.character(actors), as.character(actors)))
    if (nrow(history) == 0) {
        return(state)
    }
    weights <- if (weight_col %in% names(history)) history[[weight_col]] else rep(1, nrow(history))
    weights[is.na(weights)] <- 1
    decays <- exp(-beta_edges * (time - history$t))
    vals <- weights * decays
    idx_i <- match(history$i, actors)
    idx_j <- match(history$j, actors)
    keep <- !is.na(idx_i) & !is.na(idx_j) & idx_i != idx_j
    for (k in which(keep)) {
        state[idx_i[k], idx_j[k]] <- state[idx_i[k], idx_j[k]] + vals[k]
    }
    state
}

.rhem_candidate_stats <- function(risk, edge_state, features) {
    stats <- matrix(0, nrow = nrow(risk), ncol = length(features))
    colnames(stats) <- features
    for (feature in features) {
        stats[, feature] <- switch(feature,
            intercept = 1,
            repetition = edge_state[cbind(as.character(risk$i), as.character(risk$j))],
            reciprocity = edge_state[cbind(as.character(risk$j), as.character(risk$i))],
            sender_activity = rowSums(edge_state)[as.character(risk$i)],
            receiver_activity = colSums(edge_state)[as.character(risk$j)],
            stop("Unsupported RHEM hit feature: ", feature)
        )
    }
    stats
}

.rhem_softmax <- function(eta) {
    eta <- eta - max(eta)
    probs <- exp(eta)
    probs / sum(probs)
}

.rhem_select_backend <- function(backend = c("auto", "simple", "ernm", "timeNet"), formula_RHS) {
    backend <- match.arg(backend)
    if (backend != "auto") {
        return(backend)
    }
    if (is.null(formula_RHS)) {
        return("simple")
    }
    if (.rhem_has_timeNet_terms(formula_RHS) || .rhem_has_legacy_valued_terms(formula_RHS)) {
        return("timeNet")
    }
    "ernm"
}

.rhem_formula_rhs <- function(formula_RHS) {
    if (inherits(formula_RHS, "formula")) {
        return(as.character(formula_RHS)[length(formula_RHS)])
    }
    as.character(formula_RHS)
}

.rhem_formula_terms <- function(formula_RHS) {
    rhs <- .rhem_formula_rhs(formula_RHS)
    attr(stats::terms(stats::as.formula(paste("~", rhs))), "term.labels")
}

.rhem_has_timeNet_terms <- function(formula_RHS) {
    if (is.null(formula_RHS)) {
        return(FALSE)
    }
    supported <- c("decayedEdgeValue", "decayedRecipValue",
                   "decayedSenderActivity", "decayedReceiverActivity",
                   "decayedTransitiveValue", "decayedCycleValue",
                   "decayedCommonSourceValue", "decayedCommonTargetValue",
                   "edgeAge", "timeSinceLastEdge")
    rhs <- .rhem_formula_rhs(formula_RHS)
    any(vapply(supported, function(term) grepl(paste0("\\b", term, "\\b"), rhs), logical(1)))
}

.rhem_legacy_valued_map <- function() {
    c(edgeValue = "decayedEdgeValue",
      recipValue = "decayedRecipValue",
      senderValueActivity = "decayedSenderActivity",
      receiverValueActivity = "decayedReceiverActivity",
      transitiveValue = "decayedTransitiveValue",
      cycleValue = "decayedCycleValue",
      commonSourceValue = "decayedCommonSourceValue",
      commonTargetValue = "decayedCommonTargetValue")
}

.rhem_has_legacy_valued_terms <- function(formula_RHS) {
    if (is.null(formula_RHS)) {
        return(FALSE)
    }
    terms <- .rhem_formula_terms(formula_RHS)
    length(terms) > 0 && all(terms %in% names(.rhem_legacy_valued_map()))
}

.rhem_timeNet_term_spec <- function(formula_RHS, params) {
    if (!requireNamespace("timeNet", quietly = TRUE)) {
        stop("The `timeNet` package is required for temporal RHEM change stats.")
    }
    if (.rhem_has_legacy_valued_terms(formula_RHS)) {
        terms <- .rhem_formula_terms(formula_RHS)
        mapped <- unname(.rhem_legacy_valued_map()[terms])
        out <- data.frame(term = mapped,
                          beta = rep(params$beta_edges %||% 0, length(terms)),
                          weight_col = rep("weight", length(terms)),
                          param_name = terms,
                          stringsAsFactors = FALSE)
        return(out)
    }
    out <- timeNet::parse_time_formula(formula_RHS)
    out$param_name <- out$term
    out
}

.rhem_timeNet_state <- function(history, actors, weight_col, timeNet_state = NULL) {
    if (!requireNamespace("timeNet", quietly = TRUE)) {
        stop("The `timeNet` package is required for temporal RHEM change stats.")
    }
    if (!is.null(timeNet_state)) {
        return(timeNet_state)
    }
    net <- timeNet::temporal_directed_net(length(actors))
    if (nrow(history) > 0) {
        tails <- match(history$i, actors)
        heads <- match(history$j, actors)
        weights <- if (weight_col %in% names(history)) history[[weight_col]] else rep(1, nrow(history))
        keep <- !is.na(tails) & !is.na(heads) & tails != heads
        if (any(keep)) {
            net$addEvents(history$t[keep], tails[keep], heads[keep], weights[keep])
        }
    }
    net
}

.rhem_timeNet_cache_key <- function(time, history, actors, term_spec, hit_delta) {
    if (nrow(history) == 0) {
        history_fingerprint <- c(0, 0, 0, 0, 0)
    } else {
        history_fingerprint <- c(
            nrow(history),
            sum(history$t),
            sum(history$t * history$t),
            sum(as.numeric(history$i) * seq_len(nrow(history))),
            sum(as.numeric(history$j) * seq_len(nrow(history)))
        )
    }
    paste(c("timeNet",
            "time", format(time, digits = 17),
            "terms", paste(term_spec$term, collapse = ","),
            "params", paste(term_spec$param_name %||% term_spec$term, collapse = ","),
            "beta", paste(format(term_spec$beta, digits = 17), collapse = ","),
            "hit_delta", format(hit_delta, digits = 17),
            "actors", paste(actors, collapse = ","),
            "history", paste(format(history_fingerprint, digits = 17, scientific = TRUE),
                             collapse = ",")),
          collapse = "\r")
}

.rhem_timeNet_change_stats <- function(temporal_net,
                                       time,
                                       risk,
                                       actors,
                                       term_spec,
                                       hit_delta,
                                       history,
                                       stats_cache = NULL) {
    cache_key <- NULL
    if (!is.null(stats_cache)) {
        if (!is.environment(stats_cache)) {
            stop("`rhem_stats_cache` must be an environment when supplied.")
        }
        cache_key <- .rhem_timeNet_cache_key(time, history, actors, term_spec, hit_delta)
        if (exists(cache_key, envir = stats_cache, inherits = FALSE)) {
            return(get(cache_key, envir = stats_cache, inherits = FALSE))
        }
    }
    stats <- timeNet::time_change_stats(
        temporal_net,
        time,
        match(risk$i, actors),
        match(risk$j, actors),
        term_spec,
        hit_delta
    )
    stats <- as.matrix(stats)
    if (!is.null(stats_cache)) {
        assign(cache_key, stats, envir = stats_cache)
    }
    stats
}

.rhem_timeNet_pmf <- function(time,
                              params,
                              history,
                              current,
                              actors,
                              formula_RHS,
                              generate_mark,
                              generate_density,
                              grad,
                              weight_col,
                              hit_delta,
                              rhem_streaming_threshold,
                              rhem_stats_cache = NULL,
                              timeNet_state = NULL) {
    term_spec <- .rhem_timeNet_term_spec(formula_RHS, params)
    param_names <- term_spec$param_name %||% term_spec$term
    theta <- .rhem_theta(params, param_names)
    temporal_net <- .rhem_timeNet_state(history, actors, weight_col, timeNet_state)
    risk_size <- length(actors) * (length(actors) - 1)

    if (generate_density && !grad && !generate_mark && !is.null(current) &&
        is.finite(rhem_streaming_threshold) && risk_size > rhem_streaming_threshold) {
        tail_idx <- match(current$i[1], actors)
        head_idx <- match(current$j[1], actors)
        log_prob <- timeNet::time_log_prob(
            temporal_net, time, tail_idx, head_idx, theta, term_spec, hit_delta
        )
        return(list(mark_density = exp(log_prob[["log_probability"]]),
                    log_mark_density = log_prob[["log_probability"]],
                    edge_probs = NULL,
                    mark_grad = NULL,
                    decay_grad = 0,
                    mark_sample = NULL,
                    mark_sample_density = NULL,
                    log_mark_sample_density = NULL))
    }

    risk <- .rhem_risk_set(actors)
    stats <- .rhem_timeNet_change_stats(temporal_net, time, risk, actors, term_spec,
                                        hit_delta, history, rhem_stats_cache)
    colnames(stats) <- param_names
    eta <- as.vector(stats %*% theta)
    probs <- .rhem_softmax(eta)

    out <- list(mark_density = 1,
                log_mark_density = 0,
                edge_probs = probs,
                mark_grad = NULL,
                decay_grad = 0)
    if (generate_density && !is.null(current)) {
        idx <- which(risk$i == current$i[1] & risk$j == current$j[1])
        if (length(idx) != 1) {
            stop("Observed repeated-hit mark is not in the RHEM risk set.")
        }
        log_density <- log(probs[idx])
        out$mark_density <- exp(log_density)
        out$log_mark_density <- log_density
        if (grad) {
            expected_stats <- colSums(stats * probs)
            out$mark_grad <- stats[idx, ] - expected_stats
        }
    }
    if (generate_mark) {
        idx <- sample.int(nrow(risk), size = 1, prob = probs)
        sample_hit <- data.frame(t = time, i = risk$i[idx], j = risk$j[idx])
        sample_hit[[weight_col]] <- 1
        out$mark_sample <- rbind(history, sample_hit)
        out$mark_sample_density <- probs[idx]
        out$log_mark_sample_density <- log(probs[idx])
    } else {
        out$mark_sample <- NULL
        out$mark_sample_density <- NULL
        out$log_mark_sample_density <- NULL
    }
    out
}

.rhem_fast_valued_terms <- function(formula_RHS) {
    supported <- c("edgeValue", "recipValue", "senderValueActivity",
                   "receiverValueActivity", "transitiveValue", "cycleValue",
                   "commonSourceValue", "commonTargetValue")
    terms <- .rhem_formula_terms(formula_RHS)
    if (length(terms) == 0 || !all(terms %in% supported)) {
        return(NULL)
    }
    terms
}

.rhem_ernm_cache_key <- function(edge_state, actors, formula_RHS, hit_delta) {
    edge_vec <- as.vector(edge_state)
    weights <- seq_along(edge_vec)
    paste(c("formula", as.character(formula_RHS),
            "hit_delta", format(hit_delta, digits = 17),
            "actors", as.character(actors),
            "dim", dim(edge_state),
            "nnz", sum(edge_vec != 0),
            "sum", format(sum(edge_vec), digits = 17, scientific = TRUE),
            "sumsq", format(sum(edge_vec * edge_vec), digits = 17, scientific = TRUE),
            "weighted", format(sum(edge_vec * weights), digits = 17, scientific = TRUE),
            "weighted_squares", format(sum(edge_vec * weights * weights), digits = 17, scientific = TRUE)),
          collapse = "\r")
}

.rhem_ernm_change_stats <- function(risk, edge_state, actors, formula_RHS, hit_delta,
                                    stats_cache = NULL) {
    if (!requireNamespace("ernm", quietly = TRUE)) {
        stop("The `ernm` package is required for formula-based RHEM change stats.")
    }
    cache_key <- NULL
    if (!is.null(stats_cache)) {
        if (!is.environment(stats_cache)) {
            stop("`rhem_stats_cache` must be an environment when supplied.")
        }
        cache_key <- .rhem_ernm_cache_key(edge_state, actors, formula_RHS, hit_delta)
        if (exists(cache_key, envir = stats_cache, inherits = FALSE)) {
            return(get(cache_key, envir = stats_cache, inherits = FALSE))
        }
    }

    state_net <- .rhem_state_network(edge_state, actors)
    model <- ernm::createCppModel(stats::as.formula(paste("state_net ~", formula_RHS)))
    tails <- match(risk$i, actors)
    heads <- match(risk$j, actors)
    stats <- model$computeChangeStats(tails, heads)
    stats <- as.matrix(stats)
    if (is.null(colnames(stats))) {
        stat_names <- names(model$statistics())
        colnames(stats) <- if (!is.null(stat_names) && length(stat_names) == ncol(stats)) {
            stat_names
        } else {
            paste0("stat", seq_len(ncol(stats)))
        }
    }
    if (!is.null(stats_cache)) {
        assign(cache_key, stats, envir = stats_cache)
    }
    stats
}

.rhem_state_network <- function(edge_state, actors) {
    n <- length(actors)
    net <- network::network(matrix(0, nrow = n, ncol = n), directed = TRUE)
    network::set.vertex.attribute(net, "vertex.names", as.character(actors))
    nz <- which(edge_state > 0, arr.ind = TRUE)
    if (nrow(nz) > 0) {
        network::add.edges(net, tail = nz[, 1], head = nz[, 2])
        network::set.edge.attribute(net, "value", edge_state[nz])
    }
    net
}


#' Internal function to prepare for mark PMFs
#' @inheritParams PMF_mark
#' @noRd
mark_setup <- function(mark, mark_filtration, time){
    if(is.null(mark)){
        mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
    }
    last_net <- filtration_to_net(mark_filtration, time, equals = FALSE)
    new_net <- last_net
    if(last_net %n% 'n' != 0){
        new_nodes <- mark %n% 'n'
        old_nodes <- last_net %n% 'n'
        network::add.vertices(new_net, nv = new_nodes)
        set.vertex.attribute(new_net, "time", c(get.vertex.attribute(last_net, "time"), rep(time, new_nodes)))
    } else {
        last_net <- NULL
        new_net <- network(matrix(0, 1, 1), directed = FALSE)
        set.vertex.attribute(new_net, "time", time)
        old_nodes <- 0
        new_nodes <- 1
    }
    ## get the possible edges for the given truncation
    
    # ====================
    # THIS IS A BUG WHEN old_nodes = new_nodes!!!!!
    # ====================
    poss_tails <- seq.int((old_nodes + 1), new_nodes)
    
    poss_tails <- poss_tails[poss_tails>0]
    poss_heads <- 1:old_nodes
    poss_heads <- poss_heads[poss_heads>0]
    poss_edges <- expand.grid(poss_tails,poss_heads)
    poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
    tails <- poss_edges[,1]
    heads <- poss_edges[,2]
    ## only consider edges that were not already in the old network
    if(!is.null(last_net)){
        in_old_net <- sapply(1:length(heads),function(i){
            length(get.edgeIDs(last_net, heads[i],tails[i])) !=0
        })
        tails <- tails[!in_old_net]
        heads <- heads[!in_old_net]
    }
    
    return(list(mark = mark, new_net = new_net, last_net = last_net,
                poss_tails = poss_tails, poss_heads = poss_heads,
                poss_edges = poss_edges,
                tails = tails, heads = heads, new_nodes = new_nodes,
                old_nodes = old_nodes))
}
#' Internal function to prepare for bipartite mark PMFs
#' @inheritParams PMF_mark
#' @noRd
mark_setup_bipartite <- function(mark = NULL, mark_filtration, time){
    if(is.null(mark)){
        mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
    }
    last_net <- filtration_to_net(mark_filtration, time, equals = FALSE)
    new_net <- last_net
    if(is.null(last_net) || (last_net %n% "n") == 0){
        last_net <- NULL
        new_net <- network::network(matrix(0, 1, 1), directed = FALSE, bipartite = 0)
        set.vertex.attribute(new_net, "time", time)
        set.vertex.attribute(new_net, "role", "event")  
        old_nodes <- 0
        new_nodes <- 1
    } else {
        new_nodes <- mark %n% "n"
        old_nodes <- last_net %n% "n"
        network::add.vertices(new_net, nv = new_nodes - old_nodes)
        set.vertex.attribute(new_net, "time",
                             c(get.vertex.attribute(last_net, "time"),
                               rep(time, new_nodes - old_nodes)))
        roles <- get.vertex.attribute(last_net, "role")
        if(is.null(roles)) roles <- rep("perp", old_nodes)
        new_roles <- c(roles, rep("event", new_nodes - old_nodes))
        set.vertex.attribute(new_net, "role", new_roles)
    }
    if(is.null(last_net)){
        last_net <- network::network(matrix(0, 0, 0), directed = FALSE, bipartite = 0)
        set.vertex.attribute(last_net, "time", numeric(0))
        set.vertex.attribute(last_net, "role", character(0))
    }
    perp_nodes <- which(get.vertex.attribute(new_net, "role") == "perp")
    event_nodes <- which(get.vertex.attribute(new_net, "role") == "event")
    poss_tails <- seq.int(from = old_nodes + 1, to = new_nodes)
    poss_tails <- poss_tails[poss_tails > 0]
    poss_heads <- perp_nodes
    poss_heads <- poss_heads[poss_heads > 0]
    poss_edges <- expand.grid(poss_tails, poss_heads)
    colnames(poss_edges) <- c("tail", "head")
    tails <- poss_edges[, "tail"]
    heads <- poss_edges[, "head"]
    if(!is.null(last_net) && length(heads) > 0){
        in_old_net <- sapply(1:length(heads), function(i){
            length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
        })
        tails <- tails[!in_old_net]
        heads <- heads[!in_old_net]
    }
    return(list(mark = mark,
                new_net = new_net,
                last_net = last_net,
                poss_tails = poss_tails,
                poss_heads = poss_heads,
                poss_edges = poss_edges,
                tails = tails,
                heads = heads,
                new_nodes = new_nodes,
                old_nodes = old_nodes))
}

#' Mark PMF for Barabási–Albert-style (degree-weighted) attachment
#'
#' Probability mass function for the mark at each event: one new node and K edges from that node to existing nodes,
#' where K ~ Poisson(m). So \code{m} is the expected number of edges added per time step. Edges are sampled
#' without replacement with probabilities proportional to (degree * time decay); the likelihood uses the
#' approximation that the probability of the edge set is the product of those degree-based Bernoulli probabilities.
#'
#' @param time Current event time.
#' @param params List with \code{beta_edges}, \code{m} (expected edges per event; default 1), and optionally
#'   \code{beta_overall}, \code{K}, \code{mu}, \code{vertex_categorical}, \code{vertex_categorical_levels}.
#' @param mark_filtration Observed network up to \code{time}.
#' @param mark Optional network state at \code{time}; if \code{NULL}, derived from \code{mark_filtration}.
#' @param generate_mark If \code{TRUE}, sample K ~ Poisson(m) and then K distinct edges (default \code{FALSE}).
#' @param generate_density If \code{TRUE} (default), compute the log-density of the mark.
#' @param new_edge_hash Optional hash of existing edges for fast lookup.
#' @param truncation Optional integer cap on the number of edges per event.
#' @param mark_decay Character string controlling how temporal weights decay.
#'   One of \code{"node_entrance"} (default) or \code{"activity"}.
#' @param ... Additional arguments (currently unused).
#' @details
#' The \code{params} list may optionally include:
#' \describe{
#'   \item{vertex_categorical}{Named list of multinomial proportions per vertex attribute.}
#'   \item{vertex_categorical_levels}{Named list of level names per attribute; last level is reference.}
#' }
#' @return List with \code{log_mark_density}, \code{log_density_func}, and optionally sampled mark / probabilities.
#' @examples
#' \donttest{
#' params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
#' net <- network::network.initialize(5, directed = FALSE)
#' network::set.vertex.attribute(net, "time", seq(0.1, 0.5, length.out = 5))
#' # Compute mark density for a new node at time 0.6
#' pmf <- PMF_mark_BA(0.6, params, net)
#' pmf$log_mark_density
#' }
#' @seealso \code{\link[network]{network}}, \code{\link[network]{add.vertices}}
#' @rdname PMF_mark_BA
#' @export
PMF_mark_BA <- function(time,
                        params,
                        mark_filtration,
                        mark = NULL,
                        generate_mark = FALSE,
                        generate_density = TRUE,
                        new_edge_hash = NULL,
                        truncation = NULL,
                        mark_decay = 'node_entrance',
                        ...){
  # Expected number of edges per event (Poisson rate)
  m_val <- if (!is.null(params$m) && is.numeric(params$m) && length(params$m) == 1L && is.finite(params$m) && params$m > 0) params$m else 1

  if(is.null(mark)){
    mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
  }
  last_net <- filtration_to_net(mark_filtration, time, equals = FALSE)
  new_net <- last_net

  if(last_net %n% 'n' != 0){
    new_nodes <- mark %n% 'n'
    old_nodes <- last_net %n% 'n'
    add.vertices(new_net, nv = new_nodes)
    set.vertex.attribute(new_net, "time", c(get.vertex.attribute(last_net, "time"), rep(time, new_nodes)))
  } else {
    last_net <- NULL
    new_net <- network(matrix(0, 1, 1), directed = FALSE)
    set.vertex.attribute(new_net, "time", time)
    old_nodes <- 0
    new_nodes <- 1
  }

  # get the possible edges: new nodes -> old nodes, with optional truncation
  if (!is.null(truncation) && !is.null(last_net) && old_nodes > truncation) {
    # Truncate which old nodes are considered as attachment targets
    if (mark_decay == "activity") {
      activity_times <- get_latest_times(last_net)
      # Pick the truncation most recently active old nodes
      ord <- order(activity_times[seq_len(old_nodes)], seq_len(old_nodes), decreasing = TRUE)
      eligible_heads <- sort(ord[seq_len(min(truncation, old_nodes))])
    } else {
      # node_entrance: take the most recently entered old nodes
      eligible_heads <- seq(max(1L, old_nodes - truncation + 1L), old_nodes)
    }
  } else {
    eligible_heads <- seq_len(old_nodes)
  }
  poss_tails <- ((old_nodes+1) : new_nodes)
  poss_tails <- poss_tails[poss_tails>0]
  poss_heads <- eligible_heads[eligible_heads > 0]
  poss_edges <- expand.grid(poss_tails,poss_heads)
  poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
  tails <- poss_edges[,1]
  heads <- poss_edges[,2]
  
  # only consider edges that were not already in the old network
  if(!is.null(last_net)){
    in_old_net <- sapply(seq_along(heads), function(i) {
      length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
    })
    tails <- tails[!in_old_net]
    heads <- heads[!in_old_net]
  }

  if(!is.null(last_net) && (last_net %n% 'n' > 2)){
    times <- get.vertex.attribute(last_net, "time")
    node_degrees <- degree(last_net)
    degs <- node_degrees * exp(-params$beta_edges * (time - times))
    degs[is.na(degs) | is.nan(degs)] <- 0
    total_deg <- sum(degs)
    if (total_deg <= 0 || !is.finite(total_deg)) {
      probs <- rep(1 / length(heads), length(heads))
    } else {
      probs <- degs[heads] / total_deg
    }
    probs[is.na(probs) | is.nan(probs)] <- 0
    probs[probs < 0] <- 0
    probs[probs > 1] <- 1
    if (length(probs) > 0 && all(probs == 0)) probs[] <- 1 / length(probs)
  } else {
    node_degrees <- NULL
    probs <- rep(1, length(heads))
    times <- numeric(0)  # not used when node_degrees is NULL; avoids missing 'times' in closure env
  }
  
  if(!is.null(mark) && length(probs) !=0 && (last_net %n% 'n' > 2)){
    if(is.null(new_edge_hash)){
      in_mark <- sapply(seq_along(heads), function(i) {
        length(get.edgeIDs(mark, heads[i], tails[i])) != 0
      })
    } else {
      in_mark <- has_edge(heads, tails, new_edge_hash)
    }
    K_obs <- sum(in_mark, na.rm = TRUE)
    log_poisson <- dpois(K_obs, m_val, log = TRUE)
    p_in <- pmax(probs[in_mark], .Machine$double.eps)
    p_out <- pmax(1 - probs[!in_mark], .Machine$double.eps)
    log_mark_density <- log_poisson + sum(log(p_in), na.rm = TRUE) + sum(log(p_out), na.rm = TRUE)
    observed_categorical <- list()
    if (!is.null(mark) && (new_nodes - old_nodes) > 0) {
      vcat_res <- log_categorical_density(params, mark, old_nodes, new_nodes, eps = 1e-10)
      if (!is.null(vcat_res$log_dens) && is.finite(vcat_res$log_dens)) log_mark_density <- log_mark_density + vcat_res$log_dens
      observed_categorical <- vcat_res$observed
    }
    mark_density <- exp(log_mark_density)
  } else {
    in_mark <- rep(1, length(heads))
    K_obs <- if (length(heads) > 0) sum(in_mark, na.rm = TRUE) else 0
    log_mark_density <- dpois(K_obs, m_val, log = TRUE)
    observed_categorical <- list()
    if (!is.null(mark) && (new_nodes - old_nodes) > 0) {
      vcat_res <- log_categorical_density(params, mark, old_nodes, new_nodes, eps = 1e-10)
      if (!is.null(vcat_res$log_dens) && is.finite(vcat_res$log_dens)) log_mark_density <- log_mark_density + vcat_res$log_dens
      observed_categorical <- vcat_res$observed
    }
    mark_density <- exp(log_mark_density)
  }
  
  # Pre-compute level names for vertex categorical (for closure)
  level_names_by_attr <- list()
  if (length(observed_categorical) > 0) {
    for (an in names(observed_categorical)) {
      level_names_by_attr[[an]] <- vertex_categorical_level_names(params, an, mark)
    }
  }

  # 1. Define the lightweight log-density function for BA (includes Poisson(m) for K_obs edges + vertex categorical)
  log_density_func_light <- function(params) {
    m_p <- if (!is.null(params$m) && is.numeric(params$m) && length(params$m) == 1L && is.finite(params$m) && params$m > 0) params$m else 1
    log_poisson <- dpois(K_obs, m_p, log = TRUE)
    if (!is.finite(log_poisson)) return(-1e10)
    node_dens <- 0
    if (!is.null(node_degrees)) {
      degs <- node_degrees * exp(-params$beta_edges * (time - times))
      degs[is.na(degs) | is.nan(degs)] <- 0
      total_deg <- sum(degs)
      if (!is.finite(total_deg) || total_deg <= 0) return(-1e10)
      probs <- degs[heads] / total_deg
      probs[is.na(probs) | is.nan(probs)] <- 0
      probs[probs < 0] <- 0
      probs[probs > 1] <- 1
      if (length(probs) > 0 && all(probs == 0)) probs[] <- 1 / length(probs)
      p_in <- pmax(probs[in_mark], .Machine$double.eps)
      p_out <- pmax(1 - probs[!in_mark], .Machine$double.eps)
      node_dens <- sum(log(p_in), na.rm = TRUE) + sum(log(p_out), na.rm = TRUE)
    }
    vcat <- params$vertex_categorical
    if (!is.null(vcat) && is.list(vcat) && length(observed_categorical) > 0) {
      eps_cl <- 1e-10
      for (attr_name in names(observed_categorical)) {
        if (attr_name %in% names(vcat)) {
          levs_attr <- level_names_by_attr[[attr_name]]
          if (is.null(levs_attr) && !is.null(params$vertex_categorical_levels) && attr_name %in% names(params$vertex_categorical_levels))
            levs_attr <- params$vertex_categorical_levels[[attr_name]]
          if (!is.null(levs_attr)) {
            p_attr <- expand_vertex_categorical_probs(vcat[[attr_name]], levs_attr, eps = eps_cl)
            if (!is.null(p_attr)) {
              obs <- observed_categorical[[attr_name]]
              idx <- match(obs, names(p_attr))
              idx[is.na(idx)] <- match("unknown", names(p_attr))
              idx[is.na(idx)] <- 1L
              node_dens <- node_dens + sum(log(pmax(p_attr[idx], eps_cl)))
            }
          }
        }
      }
    }
    return(log_poisson + node_dens)
  }

  # 2. Capture the environment (include K_obs for Poisson term + vertex categorical)
  environment(log_density_func_light) <- list2env(
    list(
      heads = heads,
      time = time,
      times = times,
      node_degrees = node_degrees,
      in_mark = in_mark,
      K_obs = K_obs,
      dpois = stats::dpois,
      observed_categorical = observed_categorical,
      level_names_by_attr = level_names_by_attr,
      expand_vertex_categorical_probs = expand_vertex_categorical_probs
    ),
    parent = baseenv()
  )
  
  # 3. Define the wrapper
  density_func_light <- function(params) {
    exp(log_density_func_light(params))
  }
  
  # 4. Link wrapper to the log function
  # CRITICAL: You must include 'log_density_func_light' in the environment!
  environment(density_func_light) <- list2env(
    list(log_density_func_light = log_density_func_light), 
    parent = baseenv()
  )

  if(generate_mark){
    last_net <- mark
    times <- get.vertex.attribute(last_net, "time")
    if (!is.null(last_net) && (last_net %n% 'n') > 1) {
      mark_sample <- network.copy(last_net)
      old_nodes <- last_net %n% 'n'
      new_nodes <- 1
      mark_sample <- add.vertices(mark_sample, new_nodes)
      set.vertex.attribute(mark_sample, "time", c((last_net %v% 'time'), rep(time, new_nodes)))
      mark_sample <- sample_vertex_attrs(params, last_net, mark_sample, old_nodes, new_nodes, eps = 1e-10)
      new_node_idx <- old_nodes + 1

      if (!is.null(truncation) && old_nodes > truncation) {
        if (mark_decay == "activity") {
          activity_times <- get_latest_times(last_net)
          ord <- order(activity_times[seq_len(old_nodes)], seq_len(old_nodes), decreasing = TRUE)
          eligible_heads <- sort(ord[seq_len(min(truncation, old_nodes))])
        } else {
          eligible_heads <- seq(max(1L, old_nodes - truncation + 1L), old_nodes)
        }
      } else {
        eligible_heads <- seq_len(old_nodes)
      }
      degs <- degree(last_net) * exp(-params$beta_edges * (time - times))
      degs[is.na(degs) | is.nan(degs)] <- 0
      total_deg <- sum(degs)
      if (total_deg <= 0 || !is.finite(total_deg)) {
        probs <- rep(1 / length(eligible_heads), length(eligible_heads))
      } else {
        probs <- degs[eligible_heads] / total_deg
      }
      probs[is.na(probs) | is.nan(probs)] <- 0
      probs[probs < 0] <- 0
      probs[probs > 1] <- 1
      if (length(probs) > 0 && all(probs == 0)) probs[] <- 1 / length(probs)
      # Ensure all probs > 0 so sample(..., replace = FALSE, prob = probs) never fails with "too few positive probabilities"
      eps_p <- max(.Machine$double.eps, 1e-10)
      probs[probs <= 0 | !is.finite(probs)] <- eps_p
      probs <- probs / sum(probs)

      K <- rpois(1, m_val)
      # Join to at most all eligible targets (one node added per event; edges capped by available targets)
      K <- min(K, length(eligible_heads))
      if (K > 0 && length(eligible_heads) > 0) {
        sampled_idx <- sample(length(eligible_heads), size = K, replace = FALSE, prob = probs)
        heads_to_add <- eligible_heads[sampled_idx]
        tails_to_add <- rep(new_node_idx, K)
        add.edges(mark_sample, heads_to_add, tails_to_add)
        p_chosen <- pmax(probs[sampled_idx], .Machine$double.eps)
        log_mark_sample_density <- dpois(K, m_val, log = TRUE) + sum(log(p_chosen), na.rm = TRUE)
      } else {
        log_mark_sample_density <- dpois(K, m_val, log = TRUE)
      }
      mark_sample_density <- exp(log_mark_sample_density)
    } else {
      if (is.null(last_net) || (last_net %n% 'n') < 2) {
        if (is.null(last_net)) {
          mark_sample <- network(matrix(1), directed = FALSE)
          set.vertex.attribute(mark_sample, "time", time)
        } else {
          mark_sample <- network.copy(last_net)
        }
        times <- mark_sample %v% 'time'
        mark_sample <- add.vertices(mark_sample, 1)
        if (mark_sample %n% 'n' == 2) {
          add.edges(mark_sample, 2, 1)
        }
        set.vertex.attribute(mark_sample, "time", c(times, time))
        mark_sample <- sample_vertex_attrs(params, last_net, mark_sample, mark_sample %n% 'n' - 1L, 1L, eps = 1e-10)
        log_mark_sample_density <- dpois(1, m_val, log = TRUE)
        mark_sample_density <- exp(log_mark_sample_density)
      }
    }
  } else {
    mark_sample <- new_net
    log_mark_sample_density <- 0
    mark_sample_density <- 1
  }
  
  return(list(
    mark_density = mark_density,
    log_mark_density = log_mark_density,
    density_func = density_func_light,
    log_density_func = log_density_func_light,
    edge_probs = probs,
    # mark sample:
    mark_sample = mark_sample,
    mark_sample_density = exp(log_mark_sample_density),
    log_mark_sample_density = log_mark_sample_density
  ))
}

#' Mark PMF for change statistic (ERGM-style) attachment
#'
#' Probability mass function for the mark using an ERNM/ERGM-style model with change statistics and optional truncation.
#'
#' @param time Current event time.
#' @param params List including \code{node_lambda}, \code{CS_params}, and optionally \code{beta_edges}, \code{K}, etc.
#' @param mark_filtration Observed network up to \code{time}.
#' @param mark Optional network at \code{time}; if \code{NULL}, derived from \code{mark_filtration}.
#' @param generate_mark If \code{TRUE}, sample a new edge (default \code{FALSE}).
#' @param new_edge_hash Optional hash of existing edges for fast lookup.
#' @param formula_RHS Character RHS of the ERNM formula (e.g. \code{"edges + triangles() + star(c(2,3))"}).
#' @param truncation Truncation window: 1 = only new-to-old edges; k = edges from k steps before new nodes.
#' @return List with \code{log_mark_density}, \code{log_density_func}, and optionally sampled edge / probabilities.
#' @details
#' The \code{params} list may optionally include \code{vertex_categorical},
#' a named list of multinomial proportions per vertex attribute.
#' @seealso \code{\link[network]{network}}, \code{\link[network]{add.vertices}}, \code{\link[ernm]{as.BinaryNet}}

#' Sanitize edge probabilities: replace non-finite, clamp to \code{[eps, 1-eps]}, fallback to uniform
#' @param probs Numeric vector of probabilities.
#' @param eps Small positive floor/ceiling value (default 1e-10).
#' @param context String appended to warning messages for debugging context.
#' @return Sanitized probability vector.
#' @noRd
sanitize_probs <- function(probs, eps = 1e-10, context = "") {
  if (any(!is.finite(probs))) {
    warning("PMF_mark_CS", context, ": non-finite edge probs; replacing with 0 before clamp.")
    probs[!is.finite(probs)] <- 0
  }
  probs <- pmin(pmax(probs, eps), 1 - eps)
  if (length(probs) > 0L && all(probs <= eps)) {
    warning("PMF_mark_CS", context, ": all edge probs effectively zero; using uniform.")
    probs[] <- 1 / length(probs)
  }
  probs
}

#' Compute log multinomial density for observed discrete vertex attributes on new nodes
#' @param params Parameter list (must contain vertex_categorical and vertex_categorical_levels).
#' @param mark Network at current time (to read observed attributes).
#' @param old_nodes Number of nodes before this event.
#' @param new_nodes Total number of nodes including new arrivals.
#' @param eps Small positive value for log floor (default 1e-10).
#' @return Scalar log-density contribution.
#' @noRd
log_categorical_density <- function(params, mark, old_nodes, new_nodes, eps = 1e-10) {
  ld <- 0
  observed <- list()
  vcat <- params$vertex_categorical
  if (is.null(vcat) || !is.list(vcat) || (new_nodes - old_nodes) <= 0) {
    return(list(log_dens = ld, observed = observed))
  }
  for (attr_name in names(vcat)) {
    if (attr_name %in% list.vertex.attributes(mark)) {
      obs_vals <- (mark %v% attr_name)[(old_nodes + 1):new_nodes]
      observed[[attr_name]] <- obs_vals
      level_names <- vertex_categorical_level_names(params, attr_name, mark)
      p <- expand_vertex_categorical_probs(vcat[[attr_name]], level_names, eps = eps)
      if (!is.null(p)) {
        idx <- match(obs_vals, names(p))
        idx[is.na(idx)] <- match("unknown", names(p))
        idx[is.na(idx)] <- 1L
        ld <- ld + sum(log(pmax(p[idx], eps)))
      }
    }
  }
  list(log_dens = ld, observed = observed)
}

#' Sample discrete vertex attributes for new nodes and set them on the network
#' @param params Parameter list (must contain vertex_categorical and vertex_categorical_levels).
#' @param last_net Network before this event (to get existing attribute values).
#' @param mark_sample Network to set attributes on.
#' @param old_nodes Number of nodes before this event.
#' @param new_nodes Number of new nodes added.
#' @param eps Small positive value (default 1e-10).
#' @return The updated mark_sample network (invisible).
#' @noRd
sample_vertex_attrs <- function(params, last_net, mark_sample, old_nodes, new_nodes, eps = 1e-10) {
  vcat <- params$vertex_categorical
  if (is.null(vcat) || !is.list(vcat)) return(mark_sample)
  for (attr_name in names(vcat)) {
    levs <- vertex_categorical_level_names(params, attr_name, mark_sample)
    p <- expand_vertex_categorical_probs(vcat[[attr_name]], levs, eps = eps)
    if (is.null(p)) next
    levs <- names(p)
    if (any(!is.finite(p)) || sum(p) <= 0) { p <- rep(1 / length(levs), length(levs)); names(p) <- levs }
    existing <- if (is.null(last_net) || !(attr_name %in% list.vertex.attributes(last_net))) rep(levs[1L], old_nodes) else last_net %v% attr_name
    if (new_nodes > 0) {
      sampled <- sample(levs, size = new_nodes, replace = TRUE, prob = p)
      set.vertex.attribute(mark_sample, attr_name, c(existing, sampled))
    } else {
      set.vertex.attribute(mark_sample, attr_name, existing)
    }
  }
  mark_sample
}

#' Ensure all required vertex attributes are present on a network (fill missing with defaults)
#' @param params Parameter list (must contain vertex_categorical).
#' @param net Network to check/fix.
#' @return The updated network (invisible).
#' @noRd
ensure_vertex_attrs <- function(params, net) {
  nv <- network.size(net)
  if (nv == 0) return(net)
  
  # 1. Handle attributes explicitly defined in params$vertex_categorical
  vcat <- params$vertex_categorical
  if (!is.null(vcat) && is.list(vcat)) {
    for (attr_name in names(vcat)) {
      levs <- if (!is.null(params$vertex_categorical_levels) && attr_name %in% names(params$vertex_categorical_levels)) {
        params$vertex_categorical_levels[[attr_name]]
      } else { c("unknown") }
      if (is.null(levs) || length(levs) == 0) levs <- c("unknown")
      
      if (!attr_name %in% list.vertex.attributes(net)) {
        set.vertex.attribute(net, attr_name, rep(levs[1L], nv))
      } else {
        attr_vals <- net %v% attr_name
        if (length(attr_vals) < nv || any(is.na(attr_vals)) || any(attr_vals == "")) {
          if (length(attr_vals) < nv) attr_vals <- c(attr_vals, rep(levs[1L], nv - length(attr_vals)))
          attr_vals[is.na(attr_vals) | attr_vals == ""] <- levs[1L]
          set.vertex.attribute(net, attr_name, attr_vals)
        }
      }
    }
  }
  
  # 2. Safety: Sanitize ALL vertex attributes to prevent C++ segfaults in as.BinaryNet.
  # Rcpp/ERNM expects attributes to be purely numeric or character, and NO NAs.
  # Non-ASCII (e.g. accented names) or very long strings can trigger segfaults.
  all_attrs <- list.vertex.attributes(net)
  all_attrs <- setdiff(all_attrs, "na")
  
  for (a in all_attrs) {
    vals <- get.vertex.attribute(net, a)
    if (length(vals) < nv) {
      default_val <- if (is.numeric(vals)) 0 else "unknown"
      vals <- c(vals, rep(default_val, nv - length(vals)))
    }
    if (is.factor(vals)) vals <- as.character(vals)
    if (any(is.na(vals))) {
      if (is.numeric(vals)) {
        vals[is.na(vals)] <- 0
      } else {
        vals <- as.character(vals)
        vals[is.na(vals)] <- "unknown"
      }
    }
    if (!is.numeric(vals)) {
      vals <- sanitize_vertex_attr_for_binarynet(vals, default_val = "unknown", max_len = 200L)
    }
    set.vertex.attribute(net, a, vals)
  }
  
  net
}

#' Vertex attributes required by an ERNM formula (nodeMatch/nodeMix).
#' @param formula_RHS Character RHS of formula.
#' @return Character vector of attribute names, or character(0) if none.
#' @noRd
formula_vertex_attrs_PMF <- function(formula_RHS) {
  if (is.null(formula_RHS) || !nzchar(trimws(formula_RHS))) return(character(0))
  m <- gregexpr("node(?:Match|Mix)\\s*\\(\\s*['\"]([^'\"]+)['\"]", formula_RHS, perl = TRUE)[[1]]
  if (m[1] == -1) return(character(0))
  s <- attr(m, "capture.start")
  l <- attr(m, "capture.length")
  unique(substring(formula_RHS, s[, 1], s[, 1] + l[, 1] - 1))
}

#' Drop unneeded vertex attributes before ERNM as.BinaryNet.
#'
#' ERNM's BinaryNet may attempt to register *all* character vertex attributes
#' as discrete variables (addDiscreteVar). High-cardinality attributes such as
#' publication dates, titles, author names, etc. can cause crashes/segfaults.
#' Keep only attributes required by the ERNM formula (nodeMatch/nodeMix) and
#' those explicitly modeled via params$vertex_categorical.
#'
#' @param net Network to modify.
#' @param formula_RHS ERNM formula RHS (character).
#' @param params Parameter list (optional; used for vertex_categorical names).
#' @return Modified net.
#' @noRd
strip_vertex_attrs_for_ernm <- function(net, formula_RHS, params = NULL) {
  if (is.null(net) || network.size(net) == 0L) return(net)
  keep <- formula_vertex_attrs_PMF(formula_RHS)
  if (!is.null(params) && !is.null(params$vertex_categorical) && is.list(params$vertex_categorical)) {
    keep <- unique(c(keep, names(params$vertex_categorical)))
  }
  # CRITICAL: Always preserve "time" and "vertex.names" — these are structural

  # attributes used by filtration_to_net, get_times, and the mark density
  # (dpois for node count). Stripping "time" breaks the likelihood entirely.
  keep <- unique(c(keep, "time", "vertex.names"))
  # Always drop the problematic placeholder attr if present.
  all_attrs <- setdiff(list.vertex.attributes(net), "na")
  drop <- setdiff(all_attrs, keep)
  if (length(drop) > 0) {
    for (a in drop) {
      tryCatch(delete.vertex.attribute(net, a), error = function(e) NULL)
    }
  }
  if ("na" %in% list.vertex.attributes(net)) {
    tryCatch(delete.vertex.attribute(net, "na"), error = function(e) NULL)
  }
  net
}

#' @rdname PMF_mark_CS
#' @examples
#' normalize_vertex_categorical_probs(c(0.4, 0.4, 0.2))
#' @export
normalize_vertex_categorical_probs <- function(probs, eps = 1e-10) {
  if (is.null(probs) || length(probs) == 0) return(NULL)
  nms <- names(probs)
  probs <- as.numeric(probs)
  probs[!is.finite(probs)] <- eps
  probs <- pmax(probs, eps)
  s <- sum(probs)
  if (!is.finite(s) || s <= 0) return(NULL)
  probs <- probs / s
  if (!is.null(nms)) names(probs) <- nms
  probs
}

#' Expand n-1 vertex categorical probabilities to full n probabilities
#'
#' Given n-1 probabilities for the first n-1 levels and a vector of all n level
#' names (last is the reference), returns a named vector of length n where the
#' last level gets probability \code{1 - sum(p_n1)}.
#'
#' @param p_n1 Named numeric vector of length n-1 (probabilities for non-reference levels).
#' @param level_names Character vector of all n level names; the last element is the reference level.
#' @param eps Small positive value used as floor for any probability (default 1e-10).
#' @return Named numeric vector of length n summing to 1, or NULL if inputs are invalid.
#' @examples
#' expand_vertex_categorical_probs(c(male=0.4, female=0.4), c("male", "female", "unknown"))
#' @export
expand_vertex_categorical_probs <- function(p_n1, level_names, eps = 1e-10) {
  if (is.null(p_n1) || length(p_n1) == 0) return(NULL)
  if (is.null(level_names) || length(level_names) < 2L) return(NULL)
  n <- length(level_names)
  if (length(p_n1) != n - 1L) return(NULL)
  p_ref <- 1 - sum(p_n1)
  full <- c(as.numeric(p_n1), p_ref)
  full <- pmax(full, eps)
  names(full) <- level_names
  full
}

#' Get level names for a vertex categorical attribute
#'
#' Retrieves the level names from \code{params$vertex_categorical_levels[[attr_name]]}.
#' Falls back to the unique values observed on the network \code{mark} if available.
#'
#' @param params Parameter list (must contain \code{vertex_categorical_levels}).
#' @param attr_name Name of the vertex attribute.
#' @param mark Optional network object to fall back on for observed levels.
#' @return Character vector of level names, or NULL.
#' @noRd
vertex_categorical_level_names <- function(params, attr_name, mark = NULL) {
  if (!is.null(params$vertex_categorical_levels) && attr_name %in% names(params$vertex_categorical_levels)) {
    return(params$vertex_categorical_levels[[attr_name]])
  }
  if (!is.null(mark) && attr_name %in% list.vertex.attributes(mark)) {
    return(sort(unique(mark %v% attr_name)))
  }
  NULL
}

#' Expected parameter names for PMF_mark_BA
#'
#' @return List with \code{required} (character vector: \code{beta_edges}, \code{m}) and
#'   \code{optional} (character vector: \code{vertex_categorical}, \code{vertex_categorical_levels}).
#' @examples
#' expected_params_PMF_mark_BA()
#' @export
expected_params_PMF_mark_BA <- function() {
  list(required = c("beta_edges", "m"), optional = c("vertex_categorical", "vertex_categorical_levels"))
}

#' Expected parameter structure for PMF_mark_CS
#'
#' Returns required parameter names and the expected length of \code{CS_params}
#' (number of change statistics from the ERNM formula).
#'
#' @param mark_filtration Observed network (filtration).
#' @param formula_RHS Character RHS of the ERNM formula (e.g. \code{"edges + triangles"}).
#' @param ... Ignored.
#' @return List with \code{required} and \code{CS_params_length} (NA if cannot be computed).
#' @examples
#' \donttest{
#' net <- network::network.initialize(5, directed = FALSE)
#' network::set.vertex.attribute(net, "time", seq(0.1, 0.5, length.out = 5))
#' expected_params_PMF_mark_CS(net, "edges + triangles")
#' }
#' @export
expected_params_PMF_mark_CS <- function(mark_filtration, formula_RHS, ...) {
  required <- c("node_lambda", "CS_params", "beta_edges")
  CS_params_length <- NA_integer_
  if (is.null(formula_RHS) || is.null(mark_filtration)) {
    return(list(required = required, CS_params_length = CS_params_length))
  }
  times <- get_times(mark_filtration)$times
  if (length(times) == 0) return(list(required = required, CS_params_length = CS_params_length))
  # use the last net since stuff ight not be added til the end
  net <- filtration_to_net(mark_filtration, times[length(times)], equals = TRUE)
  nv <- network.size(net)
  if (!is.finite(nv) || is.na(nv)) nv <- 0
  if (nv < 4) {
    add.vertices(net, 4 - nv)
    t0 <- if (nv > 0) (net %v% "time")[1] else times[1]
    set.vertex.attribute(net, "time", c(net %v% "time", rep(t0, 4 - nv)))
  }
  if ("na" %in% list.vertex.attributes(net)) {
    delete.vertex.attribute(net, "na")
  }
  # If formula uses nodeMatch('gender'), ensure network has 'gender' so createCppModel/setNetwork succeed
  if (grepl("nodeMatch\\s*\\(\\s*['\"]gender['\"]", formula_RHS) && nv > 0L) {
    if (!"gender" %in% list.vertex.attributes(net)) {
      set.vertex.attribute(net, "gender", rep("unknown", nv))
    } else {
      attr_vals <- get.vertex.attribute(net, "gender")
      if (length(attr_vals) < nv || any(is.na(attr_vals)) || any(attr_vals == "")) {
        attr_vals <- if (length(attr_vals) < nv) c(attr_vals, rep("unknown", nv - length(attr_vals))) else attr_vals
        attr_vals[is.na(attr_vals) | attr_vals == ""] <- "unknown"
        set.vertex.attribute(net, "gender", attr_vals)
      }
    }
  }
  CS_params_names <- NULL
  tryCatch({
    net_safe <- sanitize_net_for_binarynet(network::network.copy(net))
    net_safe <- strip_vertex_attrs_for_ernm(net_safe, formula_RHS, params = NULL)
    model <- createCppModel(as.formula(paste("net_safe ~ ", formula_RHS)))
    model$setNetwork(as.BinaryNet(net_safe))
    model$calculate()
    stats <- model$statistics()
    CS_params_length <- length(stats)
    CS_params_names <- names(stats)
  }, error = function(e) NULL)
  list(required = required,
       optional = c("m", "vertex_categorical", "vertex_categorical_levels"),
       CS_params_length = CS_params_length, CS_params_names = CS_params_names)
}

#' Validate parameters for the given mark PMF
#'
#' Checks that \code{params} contains required names and (for CS) that
#' \code{CS_params} length matches the number of change statistics.
#'
#' @param params List of parameters passed to the mark PMF.
#' @param PMF_mark Mark PMF function (e.g. \code{PMF_mark_BA} or \code{PMF_mark_CS}).
#' @param mark_filtration Observed network; required for CS to validate \code{CS_params} length.
#' @param ... Passed through (e.g. \code{formula_RHS} for CS).
#' @return Invisible \code{TRUE}, or an error is thrown.
#' @examples
#' \donttest{
#' params <- list(beta_edges = 0.5, m = 1)
#' validate_params_for_PMF(params, PMF_mark_BA)
#' }
#' @export
validate_params_for_PMF <- function(params, PMF_mark, mark_filtration = NULL, ...) {
  if (identical(PMF_mark, PMF_mark_BA)) {
    exp_ba <- expected_params_PMF_mark_BA()
    missing <- setdiff(exp_ba$required, names(params))
    if (length(missing) > 0) {
      stop("PMF_mark_BA requires the following parameters: ", paste(missing, collapse = ", "))
    }
    message("Mark PMF parameters validated (BA): required parameters present.")
    return(invisible(TRUE))
  }
  if (identical(PMF_mark, PMF_mark_CS)) {
    formula_RHS <- list(...)$formula_RHS
    exp_cs <- expected_params_PMF_mark_CS(mark_filtration, formula_RHS)
    missing <- setdiff(exp_cs$required, names(params))
    if (length(missing) > 0) {
      stop("PMF_mark_CS requires the following parameters: ", paste(missing, collapse = ", "))
    }
    if (!is.na(exp_cs$CS_params_length)) {
      if (length(params$CS_params) != exp_cs$CS_params_length) {
        stop("PMF_mark_CS: length(CS_params) must be ", exp_cs$CS_params_length,
             " (number of change statistics from formula_RHS), got ", length(params$CS_params))
      }
      message("Mark PMF parameters validated (CS): required parameters present; length(CS_params) = ", exp_cs$CS_params_length, " matches formula.")
    } else {
      message("Mark PMF parameters validated (CS): required parameters present.")
    }
    if (!is.null(params$vertex_categorical)) {
      if (!is.list(params$vertex_categorical)) {
        stop("PMF_mark_CS: vertex_categorical must be a list of named numeric vectors (n-1 per attribute)")
      }
      if (is.null(params$vertex_categorical_levels) || !is.list(params$vertex_categorical_levels)) {
        stop("PMF_mark_CS: when vertex_categorical is set, vertex_categorical_levels must be a list of level name vectors (e.g. list(gender = c('female', 'male', 'unknown'))); last level is reference")
      }
      for (attr_name in names(params$vertex_categorical)) {
        p <- params$vertex_categorical[[attr_name]]
        if (!is.numeric(p) || length(p) == 0) {
          stop("PMF_mark_CS: vertex_categorical$", attr_name, " must be a non-empty numeric vector (n-1 parameters)")
        }
        if (is.null(names(p))) {
          stop("PMF_mark_CS: vertex_categorical$", attr_name, " must have names (level labels for non-reference levels)")
        }
        levs <- params$vertex_categorical_levels[[attr_name]]
        if (is.null(levs) || length(levs) < 2L) {
          stop("PMF_mark_CS: vertex_categorical_levels$", attr_name, " must be a character vector of length >= 2 (last is reference)")
        }
        if (length(p) != length(levs) - 1L) {
          stop("PMF_mark_CS: vertex_categorical$", attr_name, " must have length ", length(levs) - 1L, " (n-1 for ", length(levs), " levels), got ", length(p))
        }
        if (!all(names(p) %in% levs)) {
          stop("PMF_mark_CS: names of vertex_categorical$", attr_name, " must be in vertex_categorical_levels$", attr_name)
        }
        if (any(!is.finite(p)) || any(p < 0) || sum(p) >= 1) {
          stop("PMF_mark_CS: vertex_categorical$", attr_name, " must be non-negative, finite, and sum < 1 (reference level gets 1 - sum)")
        }
      }
    }
    return(invisible(TRUE))
  }
  invisible(TRUE)
}

#' Compute candidate edges for truncation
#'
#' When \code{mark_decay = "node_entrance"} (default), selects the most recently
#' *entered* nodes (by index). When \code{mark_decay = "activity"}, selects
#' the most recently *active* nodes (latest edge or entry time).
#'
#' @param net Network to compute candidates from.
#' @param new_nodes Total number of nodes in the current mark.
#' @param old_nodes Number of nodes before this event.
#' @param truncation Maximum number of nodes to consider.
#' @param mark_decay Either \code{"node_entrance"} or \code{"activity"}.
#' @param growth_only Logical; if \code{TRUE}, only allow edges from new nodes to old nodes.
#' @return List with \code{tails} and \code{heads} integer vectors.
#' @examples
#' \donttest{
#' net <- network::network.initialize(5, directed = FALSE)
#' # Get candidates for a new node (index 6) with truncation 3
#' get_truncated_candidates(net, 6, 5, 3, "node_entrance")
#' }
#' @export
get_truncated_candidates <- function(net, new_nodes, old_nodes, truncation, mark_decay, growth_only = FALSE) {
  # Add a tiny wait to ensure it's not a race condition in PSOCK
  # (though export should handle it)
  n <- max(new_nodes, old_nodes)
  if (n == 0) return(list(tails = integer(0), heads = integer(0)))

  if (mark_decay == "activity" && !is.null(net) && (net %n% 'n') > 0) {
    # Select the truncation most recently active nodes
    activity_times <- get_latest_times(net)
    # Include any new nodes (they get time = current event time, already set on net)
    # Rank nodes by activity time (most recent first); break ties by index (higher = newer)
    node_ids <- seq_len(n)
    ord <- order(activity_times[seq_len(n)], node_ids, decreasing = TRUE)
    active_nodes <- sort(ord[seq_len(min(truncation, n))])
    # All pairs among active nodes
    if (length(active_nodes) < 2) return(list(tails = integer(0), heads = integer(0)))
    poss_edges <- expand.grid(active_nodes, active_nodes)
    poss_edges <- poss_edges[poss_edges[, 1] > poss_edges[, 2], , drop = FALSE]
    tails <- poss_edges[, 1]
    heads <- poss_edges[, 2]
  } else {
    # Default: node_entrance -- use index-based window (original behavior)
    # Ensure we consider edges between new nodes and recent nodes, 
    # and among new nodes themselves.
    # poss_tails: nodes that can be the 'tail' (higher index)
    # poss_heads: nodes that can be the 'head' (lower index)
    
    # We want to consider all pairs where at least one node is "new" (index > old_nodes)
    # OR both nodes are within the truncation window of the most recent nodes.
    
    # Let's simplify: consider all pairs within the truncation window of the current total size.
    window_start <- max(1L, n - truncation + 1L)
    active_nodes <- window_start:n
    
    if (length(active_nodes) < 2) return(list(tails = integer(0), heads = integer(0)))
    poss_edges <- expand.grid(active_nodes, active_nodes)
    poss_edges <- poss_edges[poss_edges[, 1] > poss_edges[, 2], , drop = FALSE]
    tails <- poss_edges[, 1]
    heads <- poss_edges[, 2]
  }
  
  # Growth-only constraint: at least one node in the pair must be "new" (index > old_nodes)
  if (growth_only) {
    if (new_nodes <= old_nodes) {
      # No new nodes added; no new edges allowed
      return(list(tails = integer(0), heads = integer(0)))
    }
    is_new_edge <- (tails > old_nodes) | (heads > old_nodes)
    tails <- tails[is_new_edge]
    heads <- heads[is_new_edge]
  }
  
  list(tails = tails, heads = heads)
}

#' Mark probability mass function using Change Statistics (CS/ERNM)
#'
#' \code{...} can include \code{return_combined_inputs = TRUE}: then the return list
#' gets \code{combined_inputs} (change_stats, in_mark, diffs, etc.) for this event so the
#' caller can build one combined intensity closure (saves closure envs; same data, no memory blow-up).
#'
#' @param time Current event time.
#' @param params List with \code{node_lambda}, \code{CS_params}, \code{beta_edges},
#'   and optionally \code{vertex_categorical}, \code{vertex_categorical_levels}.
#' @param mark_filtration Observed network (filtration) up to \code{time}.
#' @param mark Optional network state at \code{time}; if \code{NULL}, derived from \code{mark_filtration}.
#' @param generate_mark If \code{TRUE}, generate a random mark sample.
#' @param generate_density If \code{TRUE} (default), compute the log-density of the mark.
#' @param new_edge_hash Optional hash of existing edges for fast lookup.
#' @param formula_RHS Character RHS of the ERNM formula (e.g. \code{"edges + triangles"}).
#' @param truncation Maximum number of nodes to consider for edge candidates (default 1).
#' @param mark_decay Character string controlling how temporal weights decay.
#'   One of \code{"node_entrance"} (default) or \code{"activity"}.
#' @param growth_only Logical; if \code{TRUE}, edges only form when a node enters the network.
#'   (At least one node in the pair must be a new entrant). Default \code{FALSE}.
#' @param model Optional pre-built ERNM model object to reuse.
#' @param max_node_time Optional maximum node time for temporal truncation.
#' @param ... Additional arguments; pass \code{return_combined_inputs = TRUE} to
#'   include change-statistic inputs in the return list.
#' @param probs Named numeric vector of vertex categorical probabilities (used by helpers).
#' @param eps Small positive value for probability clamping (default 1e-10).
#' @return List with \code{log_mark_density}, \code{log_density_func}, and optionally sampled mark / probabilities.
#' @examples
#' \donttest{
#' params <- list(node_lambda = 0.5, CS_params = c(-5, 0.5), beta_edges = 0.5)
#' net <- network::network.initialize(5, directed = FALSE)
#' network::set.vertex.attribute(net, "time", seq(0.1, 0.5, length.out = 5))
#' # Compute mark density for a new node at time 0.6 with triangles
#' pmf <- PMF_mark_CS(0.6, params, net, formula_RHS = "edges + triangles", truncation = 30)
#' pmf$log_mark_density
#' }
#' @rdname PMF_mark_CS
#' @export
PMF_mark_CS <- function(time,
                        params,
                        mark_filtration,
                        mark = NULL,
                        generate_mark = FALSE,
                        generate_density = TRUE,
                        new_edge_hash = NULL,
                        formula_RHS,
                        truncation = 1,
                        mark_decay = 'node_entrance',
                        growth_only = FALSE,
                        model = NULL,
                        max_node_time = NULL,
                        ...
){
  if (growth_only && !is.null(formula_RHS)) {
    zero_terms <- c("triangles", "triangle", "gwesp", "gwdsp", "esp", "dsp",
                     "ttriple", "ctriple", "kstar")
    found <- zero_terms[sapply(zero_terms, function(t) grepl(t, formula_RHS, ignore.case = TRUE))]
    if (length(found) > 0) {
      warn_key <- paste0("PMF_mark_CS_growth_only_zero_terms_", paste(found, collapse = "_"))
      if (is.null(getOption(warn_key))) {
        warning(
          "growth_only = TRUE: formula contains terms that will be structurally zero ",
          "when edges are evaluated independently from a degree-0 new node: ",
          paste(found, collapse = ", "), ". ",
          "These change statistics will always be 0 and their parameters unidentified. ",
          "Consider using only degree-based terms (e.g. edges, gwdegree) and ",
          "node-level covariates (e.g. nodeMatch, nodeMix, nodeCov).",
          call. = FALSE
        )
        options(setNames(list(TRUE), warn_key))
      }
    }
  }
  eps <- 1e-10  # used for probability clamping and safe log (CS safety)
  # Optional m parameter: when present, number of edges per event ~ Poisson(m)
  # and CS probabilities act as sampling weights (consistent with BA model).
  use_m <- !is.null(params$m) && is.numeric(params$m) && length(params$m) == 1L && is.finite(params$m) && params$m > 0
  m_val <- if (use_m) params$m else NULL

  if(is.null(mark)){
    mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
  }
  last_net <- filtration_to_net(mark_filtration, time, equals = FALSE)
  new_net <- last_net
  if(is.null(max_node_time)){
    max_node_time <- Inf
  }
  
  if(last_net %n% 'n' != 0){
    new_nodes <- mark %n% 'n'
    old_nodes <- last_net %n% 'n'
    add.vertices(new_net,new_nodes-old_nodes)
    set.vertex.attribute(new_net,"time",c(last_net %v% 'time',rep(time,new_nodes-old_nodes)))
  }else{
    last_net <- NULL
    new_net <- network(matrix(1),directed = FALSE)
    delete.vertex.attribute(new_net,'na')
    set.vertex.attribute(new_net,"time",time)
    old_nodes <- 0
    new_nodes <- 1
  }
  # get the possible edges for the given truncation:
  cands <- get_truncated_candidates(new_net, new_nodes, old_nodes, truncation, mark_decay, growth_only = growth_only)
  tails <- cands$tails
  heads <- cands$heads

  # only consider edges that were not already in the old network
  if(!is.null(last_net) & length(heads) != 0){
    in_old_net <- sapply(seq_along(heads), function(i) {
      length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
    })
    tails <- tails[!in_old_net]
    heads <- heads[!in_old_net]
  }

  if(!is.null(last_net) & generate_density){
    if(last_net %n% 'n' > 0){
      # If new net has less than 4 nodes add some (ERNM C++ safety)
      if(new_net %n% 'n' < 4){
        old_new_net <- new_net
        new_net <- add.vertices(new_net,4 - (new_net %n% 'n'))
      }else{
        old_new_net <- new_net
      }
      
      if(max(tails)>new_net %n% 'n'){
        stop("accidently adding a edge into the network that doesn't have that node yet")
      }
      
    # Copy discrete vertex attributes from mark so change stats (e.g. nodeMix) are correct
    if (!is.null(params$vertex_categorical) && is.list(params$vertex_categorical)) {
      nv <- network.size(new_net)
      for (attr_name in names(params$vertex_categorical)) {
        if (attr_name %in% list.vertex.attributes(mark)) {
          g <- mark %v% attr_name
          nm <- length(g)
          if (nv <= nm) {
            set.vertex.attribute(new_net, attr_name, g[seq_len(nv)])
          } else {
            levs <- vertex_categorical_level_names(params, attr_name, mark)
            if (is.null(levs) || length(levs) == 0) levs <- "unknown"
            # Sanitize: ensure no NAs in the attribute vector
            g_padded <- c(g, rep(levs[1L], nv - nm))
            g_padded[is.na(g_padded)] <- levs[1L]
            set.vertex.attribute(new_net, attr_name, g_padded)
          }
        }
      }
    }
      
      # delete NAs to prevent C++ using them
      delete.vertex.attribute(new_net,'na')
      if(is.null(model)){
        model <- createCppModel(as.formula(paste("new_net ~ ",formula_RHS)))
      }else{
        # CRITICAL: Strip high-cardinality metadata vertex attrs before BinaryNet conversion.
        new_net <- strip_vertex_attrs_for_ernm(new_net, formula_RHS, params = params)
        sanitize_net_for_binarynet(new_net)
        model$setNetwork(as.BinaryNet(new_net))
      }
      new_net <- old_new_net
      model$calculate()
      stat <- model$statistics()
      change_stats <- model$computeChangeStats(tails, heads)
      n_cs <- length(params$CS_params)
      if (ncol(change_stats) != n_cs) {
        change_stats <- matrix(0, nrow = NROW(change_stats), ncol = n_cs)
      }
      
      eta <- as.vector(change_stats %*% params$CS_params)
      probs <- sanitize_probs(stats::plogis(eta), eps, " (density, post-logistic)")
      # use either node times or last node activity:
      if(mark_decay == 'activity'){
        node_times <- get_latest_times(new_net)
        # For activity decay, use the MOST RECENTLY active of the two endpoints
        # so that an edge is viable as long as at least one endpoint is active.
        # This breaks the cold-start feedback loop where stale nodes could never
        # receive new edges.
        diffs <- time - pmax(node_times[heads], node_times[tails])
      } else {
        # For node_entrance: decay from the head (lower-index / older) node's
        # entrance time, matching the paper's specification: exp(-τ(t - t_i)).
        node_times <- new_net %v% 'time'
        diffs <- time - node_times[heads]
      }
      # --- Safety: sanitize diffs/factor before multiplying probs ---
      if (any(!is.finite(diffs))) {
        warning("PMF_mark_CS: non-finite time diffs in density path; replacing with 0.")
        diffs[!is.finite(diffs)] <- 0
      }
      factor <- exp(-params$beta_edges*(diffs))
      if (any(!is.finite(factor))) {
        warning("PMF_mark_CS: non-finite decay factor in density path; replacing with 1.")
        factor[!is.finite(factor)] <- 1
      }
      probs <- sanitize_probs(probs * factor, eps, " (density, post-decay)")

      if(length(probs)==0){
        in_mark <- logical(0)
        probs <- NULL
      }
      
    }else{
      change_stats <- matrix(0, nrow = 0, ncol = length(params$CS_params))
      in_mark <- logical(0)
      times <- last_net %v% 'time'
      probs <- c(1)
    }
  }else{
    change_stats <- matrix(0, nrow = 0, ncol = length(params$CS_params))
    in_mark <- logical(0)
    probs <- NULL
  }
  
  if(!is.null(mark) & !is.null(probs)){
    if(is.null(new_edge_hash)){
      in_mark <- sapply(seq_along(heads), function(i) {
        length(get.edgeIDs(mark, heads[i], tails[i])) != 0
      })
    }else{
      in_mark <- has_edge(heads,tails,new_edge_hash)
    }
    
    if(length(probs)==1){
      log_mark_density <- 0
      mark_density <-1
    }else{
      if(time >max_node_time){
        node_dens <- 0
      }else{
        dval <- dpois(new_nodes-old_nodes, params$node_lambda)
        if (!is.finite(dval) || dval <= 0) {
          warning("PMF_mark_CS: degenerate node count density (dpois=0 or non-finite); using large negative log-density.")
          node_dens <- -1e10
        } else {
          node_dens <- log(dval)
        }
      }
      # Log multinomial contribution for discrete vertex attributes on new nodes
      cat_result <- log_categorical_density(params, mark, old_nodes, new_nodes, eps)
      node_dens <- node_dens + cat_result$log_dens
      observed_categorical <- cat_result$observed
      # --- Safety: safe log with clamped probs ---
      p_in <- pmax(probs[in_mark], eps, na.rm = TRUE)
      p_out <- pmax(1 - probs[!in_mark], eps, na.rm = TRUE)
      if (any(!is.finite(p_in)) || any(!is.finite(p_out))) {
        warning("PMF_mark_CS: non-finite probs in log_mark_density; using epsilon for log.")
      }
      log_mark_density <- sum(log(p_in), na.rm = TRUE) + sum(log(p_out), na.rm = TRUE) + node_dens
      # When m is specified, add Poisson(K_obs; m) term for number of edges
      # (consistent with BA model density: Bernoulli product + Poisson count).
      if (use_m) {
        K_obs <- sum(in_mark, na.rm = TRUE)
        log_mark_density <- log_mark_density + dpois(K_obs, m_val, log = TRUE)
      }
      mark_density <- exp(log_mark_density)
    }
  }else{
    mark_density <- 1
    log_mark_density <- 0
  }
  
  # Ensure node_dens is initialized (may not be set in degenerate/no-mark paths)
  if (!exists("node_dens", inherits = FALSE)) node_dens <- 0
  # K_obs for Poisson(m) term in closure (only meaningful when use_m = TRUE)
  if (!exists("K_obs", inherits = FALSE)) K_obs <- sum(in_mark, na.rm = TRUE)
  
  # log_density_func_light: full version with decay, vertex_categorical, safety clamping.
  # (Environment is rebound to a minimal env after definition.)
  # Use plogis from env (not stats::) so closure is self-contained when serialized to parallel workers.
  log_density_func_light <- function(params) {
    # Everything it needs will come from its environment:
    # change_stats, in_mark, diffs, new_nodes, old_nodes, time, max_node_time, degenerate_edges
    
    if (degenerate_edges) {
      return(0)
    }
    
    eta    <- as.vector(change_stats %*% params$CS_params)
    p_base <- plogis(eta)
    
    # same decay factor as direct
    p <- p_base * exp(-params$beta_edges * diffs)
    
    # 2. SAFETY CLAMP
    # Ensure p is never exactly 0 or 1. 
    # This prevents log(0) and log(1-1) errors.
    epsilon <- 1e-10
    p[p > (1 - epsilon)] <- 1 - epsilon
    p[p < epsilon] <- epsilon
    
    if (anyNA(p)) return(NA_real_)
    log_edge_part <- sum(log(p[in_mark])) + sum(log1p(-p[!in_mark]))
    # Poisson(K_obs; m) term when m-parameter is active
    if (use_m) {
      m_p <- if (!is.null(params$m) && is.numeric(params$m) && length(params$m) == 1L && is.finite(params$m) && params$m > 0) params$m else 1
      log_edge_part <- log_edge_part + dpois(K_obs, m_p, log = TRUE)
    }
    # --- Safety: handle dpois=0 or non-finite in closure (suggestion 8) ---
    node_dens <- if (!is.null(max_node_time) && time > max_node_time) {
      0
    } else {
      dval <- dpois(new_nodes - old_nodes, params$node_lambda)
      if (!is.finite(dval) || dval <= 0) -1e10 else log(dval)
    }
    vcat <- params$vertex_categorical
    if (!is.null(vcat) && is.list(vcat) && length(observed_categorical) > 0) {
      eps_cl <- 1e-10
      for (attr_name in names(observed_categorical)) {
        if (attr_name %in% names(vcat)) {
          levs_attr <- level_names_by_attr[[attr_name]]
          if (is.null(levs_attr) && !is.null(params$vertex_categorical_levels) && attr_name %in% names(params$vertex_categorical_levels))
            levs_attr <- params$vertex_categorical_levels[[attr_name]]
          if (!is.null(levs_attr)) {
            p_attr <- expand_vertex_categorical_probs(vcat[[attr_name]], levs_attr, eps = eps_cl)
            if (!is.null(p_attr)) {
              obs <- observed_categorical[[attr_name]]
              idx <- match(obs, names(p_attr))
              idx[is.na(idx)] <- match("unknown", names(p_attr))
              idx[is.na(idx)] <- 1L
              node_dens <- node_dens + sum(log(pmax(p_attr[idx], eps_cl)))
            }
          }
        }
      }
    }
    log_edge_part + node_dens
  }
  
  # Decide if you're in the same degenerate branch as the direct computation
  degenerate_edges <- is.null(probs) || length(probs) == 1L
  
  # Sanitize diffs for closure (may not exist in degenerate paths)
  if (!exists("diffs", inherits = FALSE)) diffs <- numeric(0)
  diffs_for_closure <- diffs
  if (length(diffs_for_closure) > 0L && any(!is.finite(diffs_for_closure))) {
    warning("PMF_mark_CS: non-finite diffs passed to log_density_func closure; replacing with 0.")
    diffs_for_closure[!is.finite(diffs_for_closure)] <- 0
  }
  
  # Pre-compute level names per attribute for the closure
  if (!exists("observed_categorical", inherits = FALSE)) observed_categorical <- list()
  obs_cat <- observed_categorical
  level_names_by_attr <- list()
  for (an in names(obs_cat)) level_names_by_attr[[an]] <- vertex_categorical_level_names(params, an, mark)
  
  # Now *force* a tiny environment (no local needed).
  # Include expand_vertex_categorical_probs so the closure finds it when parent is baseenv().
  environment(log_density_func_light) <- list2env(
    list(
      change_stats                    = change_stats,
      in_mark                         = in_mark,
      diffs                           = diffs_for_closure,
      new_nodes                       = new_nodes,
      old_nodes                       = old_nodes,
      time                            = time,
      max_node_time                   = max_node_time,
      degenerate_edges                = degenerate_edges,
      use_m                           = use_m,
      K_obs                           = K_obs,
      observed_categorical            = obs_cat,
      level_names_by_attr             = level_names_by_attr,
      expand_vertex_categorical_probs = expand_vertex_categorical_probs,
      dpois                           = stats::dpois,
      plogis                          = stats::plogis
    ),
    parent = baseenv()
  )
  
  density_func_light <- function(params) exp(log_density_func_light(params))
  environment(density_func_light) <- list2env(
    list(log_density_func_light = log_density_func_light),
    parent = baseenv()
  )

  # Optional: return ingredients for combined intensity (one closure over all events).
  # Caller collects these from each event and builds a single closure; saves N-1 closure envs.
  return_combined_inputs <- list(...)$return_combined_inputs
  combined_inputs <- if (isTRUE(return_combined_inputs)) {
    list(
      change_stats = change_stats,
      in_mark = in_mark,
      diffs = diffs_for_closure,
      new_nodes = new_nodes,
      old_nodes = old_nodes,
      time = time,
      max_node_time = max_node_time,
      degenerate_edges = degenerate_edges,
      use_m = use_m,
      K_obs = K_obs,
      observed_categorical = obs_cat,
      level_names_by_attr = level_names_by_attr
    )
  } else NULL

  # =============
  # generate mark
  # =============
  if(generate_mark){
    # use latest mark as baseline:
    last_net <- mark
    # use function sample a new mark
    # since we have poisson number of nodes added  we need to redo the probabilities
    if(!is.null(last_net) && (last_net %n% 'n') >= 1){
      mark_sample <- last_net
      old_nodes <- last_net %n% 'n'
      
      if(mark_sample %n% 'n' < 4){
        # ERNM model requires >= 4 nodes to compute change statistics;
        # add just enough nodes to reach 4.  Once we have >= 4 nodes the
        # Poisson(node_lambda) model takes over for subsequent events.
        new_nodes <- 4L - (mark_sample %n% 'n')
      }else{
        if(time > max_node_time){
          new_nodes <- 0
        }else{
          lam <- params$node_lambda
          if (!is.finite(lam) || lam < 0) lam <- 0
          new_nodes <- rpois(1, lam)
        }
      }
      new_nodes <- as.integer(round(new_nodes))
      if (!is.finite(new_nodes) || new_nodes < 0) new_nodes <- 0L
      # Guarantee at least 4 total nodes for ERNM (safety net)
      min_needed <- 4L - (mark_sample %n% 'n')
      if (min_needed > 0L && new_nodes < min_needed) new_nodes <- min_needed
      mark_sample <- add.vertices(mark_sample, new_nodes)
      # if(mark_sample %n% 'n' > 4){
      #   browser()
      # }
      set.vertex.attribute(mark_sample,"time",c((last_net %v% 'time'),rep(time,new_nodes)))
      # Set discrete vertex attributes for new nodes (sample from vertex_categorical; n-1 params)
      mark_sample <- sample_vertex_attrs(params, last_net, mark_sample, old_nodes, new_nodes, eps)
      new_size <- mark_sample %n% 'n'
      
      # get new poss edges (truncation based on mark_decay)
      cands <- get_truncated_candidates(mark_sample, new_size, old_nodes, truncation, mark_decay, growth_only = growth_only)
      tails <- cands$tails
      heads <- cands$heads

      # only consider edges that are not in the old net
      if(!is.null(last_net)){
        in_old_net <- sapply(seq_along(heads), function(i) {
          length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
        })
        tails <- tails[!in_old_net]
        heads <- heads[!in_old_net]
      }

      delete.vertex.attribute(mark_sample,'na')
      
      # Ensure ALL nodes have required vertex attributes before createCppModel
      mark_sample <- ensure_vertex_attrs(params, mark_sample)
      if ("na" %in% list.vertex.attributes(mark_sample)) {
        delete.vertex.attribute(mark_sample, "na")
      }
      # CRITICAL: Drop unneeded high-cardinality metadata vertex attributes (e.g. publication_date, title)
      # before any ERNM BinaryNet conversion. ERNM may register *all* character vertex attrs as discrete
      # vars (addDiscreteVar), which can segfault for very high-cardinality attributes.
      mark_sample <- strip_vertex_attrs_for_ernm(mark_sample, formula_RHS, params = params)
      
      # Reuse ERNM model per formula (avoids createCppModel every event; big speedup for nodeMatch)
      cache <- get0(".ernm_model_cache", envir = asNamespace("hawkesNet"), inherits = FALSE)
      if (is.null(cache)) {
        cache <- new.env()
        cache[[".cache_pid"]] <- Sys.getpid()
        assign(".ernm_model_cache", cache, envir = asNamespace("hawkesNet"))
      }
      # If cache exists but lacks pid marker (older installs / sourced code), initialize it.
      if (is.null(cache[[".cache_pid"]]) || length(cache[[".cache_pid"]]) != 1L || !is.finite(cache[[".cache_pid"]])) {
        cache[[".cache_pid"]] <- Sys.getpid()
      }
      # Invalidate cache in forked children: C++ pointers become invalid after fork.
      if (Sys.getpid() != cache[[".cache_pid"]]) {
        rm(list = setdiff(names(cache), ".cache_pid"), envir = cache)
        cache[[".cache_pid"]] <- Sys.getpid()
      }
      key <- formula_RHS
      # Create or re-create the model.  Use mark_sample (not an empty g0) so
      # that ERNM can find vertex attributes such as 'gender' for nodeMatch.
      # Also re-create if the cached C++ pointer is stale (e.g. after fork).
      need_create <- is.null(cache[[key]])
      if (!need_create) {
        need_create <- tryCatch({
          sanitize_net_for_binarynet(mark_sample)
          cache[[key]]$setNetwork(as.BinaryNet(mark_sample))
          cache[[key]]$calculate()
          FALSE
        }, error = function(e) TRUE)
      }
      if (need_create) {
        ms_ref <- mark_sample
        cache[[key]] <- createCppModel(as.formula(paste("ms_ref ~ ", formula_RHS)))
        sanitize_net_for_binarynet(mark_sample)
        cache[[key]]$setNetwork(as.BinaryNet(mark_sample))
      }
      model <- cache[[key]]
      sanitize_net_for_binarynet(mark_sample)
      model$setNetwork(as.BinaryNet(mark_sample))
      model$calculate()
      change_stats <- model$computeChangeStats(tails, heads)
      eta <- as.vector(change_stats %*% params$CS_params)
      probs <- sanitize_probs(stats::plogis(eta), eps, " (generate_mark, post-logistic)")

      # reset to when we did not add more edges
      # logistic regression on change stats:
      dot_list <- list(...)
      stop_on_full_network <- if ("stop_on_full_network" %in% names(dot_list)) dot_list$stop_on_full_network else TRUE
      if (length(change_stats) == 0) {
        if (stop_on_full_network) {
          stop("these parameters result in full networks - you probably don't want this")
        }
        warning("PMF_mark_CS (generate_mark): no candidate edges (full network); returning mark with no new edges (stop_on_full_network = FALSE).")
        if (new_nodes == 0) {
          nv_before <- network.size(mark_sample)
          mark_sample <- add.vertices(mark_sample, 1)
          set.vertex.attribute(mark_sample, "time", c(mark_sample %v% "time", time))
          mark_sample <- sample_vertex_attrs(params, mark_sample, mark_sample, nv_before, 1L, eps)
        }
        mark_sample_density <- 1
        log_mark_sample_density <- 0
      } else {
      if(mark_decay == 'activity'){
        node_times <- get_latest_times(mark_sample)
        # For activity decay: most recently active of the two endpoints.
        diffs <- time - pmax(node_times[heads], node_times[tails])
      } else {
        # For node_entrance: decay from head node's entrance time (paper spec).
        node_times <- mark_sample %v% 'time'
        diffs <- time - node_times[heads]
      }
      # --- Safety: sanitize diffs/factor in generate_mark (suggestion 2 & 9) ---
      if (any(!is.finite(diffs))) {
        warning("PMF_mark_CS (generate_mark): non-finite time diffs; replacing with 0.")
        diffs[!is.finite(diffs)] <- 0
      }
      factor <- exp(-params$beta_edges*(diffs))
      if (any(!is.finite(factor))) {
        warning("PMF_mark_CS (generate_mark): non-finite decay factor; replacing with 1.")
        factor[!is.finite(factor)] <- 1
      }
      probs <- sanitize_probs(factor * probs, eps, " (generate_mark, post-decay)")
      
      verbose_mark <- if ("verbose" %in% names(dot_list)) dot_list$verbose else FALSE
      if (verbose_mark) {
        n_iso <- sum(sna::degree(mark_sample, gmode = "graph") == 0)
        median_diff <- median(diffs)
        cat(sprintf("    [Mark] Cands: %d | E[edges]: %.2f | max_p: %.4e | mean_p: %.4e | med_diff: %.4f | isolates: %d/%d | new_nodes: %d | use_m: %s\n", 
                    length(probs), sum(probs), max(probs), mean(probs), median_diff,
                    n_iso, mark_sample %n% "n", new_nodes, use_m))
      }

      if (use_m) {
        # --- m-parameter mode: sample K ~ Poisson(m) edges using probs as weights ---
        # Consistent with BA model: Poisson count + weighted sampling without replacement.
        K_gen <- rpois(1, m_val)
        K_gen <- min(K_gen, length(probs))
        if (K_gen > 0 && length(probs) > 0) {
          w <- probs / sum(probs)
          w[!is.finite(w)] <- eps
          w[w <= 0] <- eps
          w <- w / sum(w)
          sampled_idx <- sample(length(probs), size = K_gen, replace = FALSE, prob = w)
          add.edges(mark_sample, heads[sampled_idx], tails[sampled_idx])
          set.edge.attribute(mark_sample, "time", c(mark_sample %e% 'time', rep(time, K_gen)))
          p_chosen <- pmax(probs[sampled_idx], eps)
          log_edge_gen <- dpois(K_gen, m_val, log = TRUE) + sum(log(p_chosen), na.rm = TRUE)
        } else {
          log_edge_gen <- dpois(K_gen, m_val, log = TRUE)
        }
        dpois_val <- dpois(new_nodes - old_nodes, params$node_lambda)
        if (!is.finite(dpois_val) || dpois_val <= 0) dpois_val <- 1e-300
        log_multinomial_sample <- log_categorical_density(params, mark_sample, old_nodes, new_size, eps)$log_dens
        log_mark_sample_density <- log_edge_gen + log(dpois_val) + log_multinomial_sample
        mark_sample_density <- exp(log_mark_sample_density)
      } else {
        # --- Original Bernoulli mode: independent draw per candidate ---
        add <- runif(length(probs)) < probs
        if (any(is.na(add))) {
          warning("PMF_mark_CS (generate_mark): NA in edge add vector; treating as FALSE.")
          add[is.na(add)] <- FALSE
        }
        n_added_edges <- sum(add)
        if (n_added_edges > 0) {
          add.edges(mark_sample, heads[add], tails[add])
          set.edge.attribute(mark_sample, "time", c(mark_sample %e% 'time', rep(time, n_added_edges)))
        }
        dpois_val <- dpois(new_nodes - old_nodes, params$node_lambda)
        if (!is.finite(dpois_val) || dpois_val <= 0) {
          warning("PMF_mark_CS (generate_mark): degenerate dpois for node count; using small positive value for density.")
          dpois_val <- 1e-300
        }
        p_add <- pmax(probs[add], eps, na.rm = TRUE)
        p_not <- pmax(1 - probs[!add], eps, na.rm = TRUE)
        if (any(!is.finite(p_add)) || any(!is.finite(p_not))) {
          warning("PMF_mark_CS (generate_mark): non-finite probs in log_mark_sample_density; using epsilon for log.")
        }
        log_multinomial_sample <- log_categorical_density(params, mark_sample, old_nodes, new_size, eps)$log_dens
        mark_sample_density <- prod(p_add) * prod(p_not) * dpois_val * exp(log_multinomial_sample)
        log_mark_sample_density <- sum(log(p_add), na.rm = TRUE) +
                                   sum(log(p_not), na.rm = TRUE) +
                                   log(dpois_val) + log_multinomial_sample
      }
      }
      }else{
        if(is.null(last_net)){
          mark_sample <- network(matrix(1),directed = FALSE)
          set.vertex.attribute(mark_sample,"time",time)
        }else{
          mark_sample <- last_net
        }
        times <- mark_sample %v% 'time'
        mark_sample <- add.vertices(mark_sample,1)
        set.vertex.attribute(mark_sample,
                             "time",
                             c(times,time))
        mark_sample <- sample_vertex_attrs(params, if (!is.null(last_net)) last_net else mark_sample, mark_sample, length(times), 1L, eps)
        mark_sample_density <- 1
        log_mark_sample_density <- 0
      }
    }else{
      mark_sample <- new_net
      mark_sample_density <- 1
      log_mark_sample_density <- 0
    }

  out <- list(
    # density of provided marks
    mark_density = mark_density,
    log_mark_density = log_mark_density,
    density_func = density_func_light,
    log_density_func = log_density_func_light,
    edge_probs = probs,
    # mark_sample
    mark_sample = mark_sample,
    mark_sample_density = mark_sample_density,
    log_mark_sample_density = log_mark_sample_density
  )
  if (!is.null(combined_inputs)) out$combined_inputs <- combined_inputs
  return(out)
}
