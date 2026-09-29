#' Specify topology feedback in the ground Hawkes intensity
#'
#' M0 gives every complete network update unit excitation weight. M1 freezes
#' each update's structural feature immediately after that whole update; M2
#' reevaluates it in the current graph. The kernel convention is
#' `K * alpha * exp(-beta_overall * age)`, so K is an amplitude, not a
#' branching ratio. Updates sharing a timestamp form one complete mark.
#'
#' With `feature = "degree"`, let U contain all vertices added by the update
#' and all endpoints of its added edges. If d is their mean current degree,
#' the feature is `d / (scale + d)` and the default link is
#' `alpha = exp(feedback_gamma * feature)`.
#' With `feature = "triangles"`, the feature is the fraction of the update's
#' added edges currently belonging to at least one triangle (zero for a
#' node-only update). Its default link is
#' `alpha = 2 * plogis(feedback_gamma * feature)`.
#' Either link can be used with either feature. Both links give alpha = 1 at
#' `feedback_gamma = 0`; the parameter defaults to zero if absent from params.
#'
#' M1 and M2 currently support undirected simple network growth, with one
#' vertex/edge `time` attribute recording each addition. They do not implement
#' deletions, revisions, directed, hypergraph or multiple-edge feedback.
#' NA vertex/edge times describe a static initial graph and contribute no
#' event. Every initial edge must connect initial vertices. All other times
#' must be finite. Original vertex indices are retained throughout feature
#' construction. M2 caches features after each update and therefore uses
#' quadratic memory in the number of updates.
#'
#' @param model One of `"M0"`, `"M1"`, `"M2"`.
#' @param feature Structural feature, `"degree"` or `"triangles"`.
#' @param link Positive amplitude link, `"exp"` or `"logistic"`; NULL chooses
#'   the default associated with the feature.
#' @param scale Positive degree saturation scale; ignored for triangles.
#' @return A `hawkes_timing` specification.
#' @examples
#' hawkes_timing("M1", "degree", scale = 5)
#' hawkes_timing("M2", "triangles")
#' @export
hawkes_timing <- function(model = c("M0", "M1", "M2"),
                          feature = c("degree", "triangles"),
                          link = NULL, scale = 1) {
  model <- match.arg(model)
  feature <- match.arg(feature)
  if (is.null(link)) link <- if (feature == "degree") "exp" else "logistic"
  link <- match.arg(link, c("exp", "logistic"))
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0)
    stop("scale must be one finite positive number")
  structure(list(model = model, feature = feature, link = link, scale = scale),
            class = "hawkes_timing")
}

.network_timing_validate <- function(timing) {
  if (!inherits(timing, "hawkes_timing") ||
      !identical(names(timing), c("model", "feature", "link", "scale")))
    stop("timing must be a specification from hawkes_timing()")
  check <- hawkes_timing(timing$model, timing$feature, timing$link, timing$scale)
  if (!identical(check, timing)) stop("invalid hawkes_timing specification")
  invisible(timing)
}

.network_timing_gamma <- function(params) {
  gamma <- params$feedback_gamma
  if (is.null(gamma)) gamma <- 0
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma))
    stop("feedback_gamma must be one finite number")
  gamma
}

.network_timing_params <- function(params) {
  if (!is.list(params)) stop("params must be a named parameter list")
  for (name in c("mu", "K", "beta_overall")) {
    value <- params[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || (name == "beta_overall" && value == 0))
      stop(name, " must be one finite ",
           if (name == "beta_overall") "positive" else "nonnegative", " number")
  }
  .network_timing_gamma(params)
  invisible(params)
}

# Uniform alpha bound on the feature range [0, 1]. This is a pathwise
# envelope, not a contraction/stationarity test for the complete marked law.
.network_timing_bound <- function(params, timing) {
  .network_timing_validate(timing)
  gamma <- .network_timing_gamma(params)
  if (timing$model == "M0") return(1)
  bound <- if (timing$link == "exp") exp(max(0, gamma)) else
    2 * stats::plogis(max(0, gamma))
  if (!is.finite(bound)) stop("feedback_gamma produces a non-finite alpha bound")
  bound
}

.network_timing_alpha <- function(features, params, timing) {
  if (timing$model == "M0") return(rep(1, length(features)))
  gamma <- .network_timing_gamma(params)
  alpha <- if (timing$link == "exp") exp(gamma * features) else
    2 * stats::plogis(gamma * features)
  if (any(!is.finite(alpha))) stop("feedback_gamma produces non-finite excitation")
  alpha
}

.network_timing_times <- function(mark_filtration) {
  if (is.null(mark_filtration)) return(numeric(0))
  if (inherits(mark_filtration, "data.frame")) {
    if (!"t" %in% names(mark_filtration)) stop("event data must contain a t column")
    times <- mark_filtration$t
    if (!is.numeric(times) || any(!is.finite(times)))
      stop("event times must be finite numeric values")
    return(sort(unique(times)))
  }
  if (!network::is.network(mark_filtration))
    stop("mark_filtration must be a network or event data.frame")
  times <- c(network::get.vertex.attribute(mark_filtration, "time"),
             network::get.edge.attribute(mark_filtration, "time"))
  # network returns logical NA when all initial times are unspecified.
  if (is.logical(times) && all(is.na(times))) times <- as.numeric(times)
  if (length(times) && (!is.numeric(times) || any(!is.finite(times[!is.na(times)])) ||
                        any(is.nan(times))))
    stop("event times must be finite numeric values (NA denotes initial structure)")
  as.numeric(sort(unique(times[!is.na(times)])))
}

# The cache is independent of mu, K, beta and feedback_gamma. For M2,
# features_after[[j]] describes marks 1:j in the graph just after event j.
# No vertex deletion/reindexing is used when reconstructing earlier graphs.
.prepare_network_timing <- function(mark_filtration, timing = hawkes_timing()) {
  .network_timing_validate(timing)
  if (timing$model != "M0" && !is.null(mark_filtration)) {
    if (!network::is.network(mark_filtration))
      stop("M1/M2 topology feedback requires an undirected simple growth network; data.frame marks are unsupported")
    unsupported <- vapply(c("directed", "multiple", "hyper", "loops"), function(a)
      isTRUE(network::get.network.attribute(mark_filtration, a)), logical(1))
    if (any(unsupported))
      stop("M1/M2 topology feedback requires an undirected simple growth network (no directed, multiple, hyper or loop edges)")
  }
  times <- .network_timing_times(mark_filtration)
  n <- length(times)
  cache <- structure(list(timing = timing, times = times,
                          birth_features = numeric(n), features_after = NULL),
                     class = "hawkes_timing_cache")
  if (timing$model == "M0" || is.null(mark_filtration)) return(cache)

  nv <- network::network.size(mark_filtration)
  node_times <- network::get.vertex.attribute(mark_filtration, "time")
  if (length(node_times) != nv) stop("each vertex must have a time attribute (NA for initial vertices)")
  ne <- network::network.edgecount(mark_filtration)
  if (ne > 0L) {
    # as.edgelist attaches the attribute to its corresponding endpoints before
    # sorting rows; get.edge.attribute alone need not have the same row order.
    edge_rows <- network::as.edgelist(mark_filtration, attrname = "time")
    if (ncol(edge_rows) != 3L || nrow(edge_rows) != ne)
      stop("each edge must have a time attribute (NA for initial edges)")
    edge_ends <- edge_rows[, 1:2, drop = FALSE]
    edge_times <- edge_rows[, 3]
    if (!is.numeric(edge_times) || any(is.nan(edge_times)) ||
        any(!is.finite(edge_times[!is.na(edge_times)])))
      stop("edge event times must be finite numeric values or NA for initial edges")
    keys <- paste(pmin(edge_ends[, 1], edge_ends[, 2]),
                  pmax(edge_ends[, 1], edge_ends[, 2]), sep = ":")
    if (any(edge_ends[, 1] == edge_ends[, 2]) || anyDuplicated(keys))
      stop("M1/M2 topology feedback requires an undirected simple growth network")
    for (column in 1:2) {
      birth <- node_times[edge_ends[, column]]
      if (any(is.na(edge_times) & !is.na(birth)) ||
          any(!is.na(birth) & !is.na(edge_times) & birth > edge_times))
        stop("an edge cannot precede the birth of either endpoint")
    }
  } else {
    edge_ends <- matrix(integer(0), ncol = 2L)
    edge_times <- numeric(0)
    keys <- character(0)
  }
  neighbors <- vector("list", nv)
  degrees <- numeric(nv)
  triangle_edge <- rep(FALSE, ne)
  present_edges <- new.env(hash = TRUE, parent = emptyenv())
  edge_key <- function(u, v) paste(min(u, v), max(u, v), sep = ":")
  insert_edge <- function(eid) {
    u <- edge_ends[eid, 1]
    v <- edge_ends[eid, 2]
    common <- intersect(neighbors[[u]], neighbors[[v]])
    if (length(common)) {
      triangle_edge[eid] <<- TRUE
      for (w in common) {
        triangle_edge[present_edges[[edge_key(u, w)]]] <<- TRUE
        triangle_edge[present_edges[[edge_key(v, w)]]] <<- TRUE
      }
    }
    neighbors[[u]] <<- c(neighbors[[u]], v)
    neighbors[[v]] <<- c(neighbors[[v]], u)
    degrees[c(u, v)] <<- degrees[c(u, v)] + 1
    present_edges[[keys[eid]]] <- eid
  }
  for (eid in which(is.na(edge_times))) insert_edge(eid)
  mark_vertices <- vector("list", n)
  mark_edges <- vector("list", n)
  feature_now <- function(i) {
    if (timing$feature == "degree") {
      d <- if (length(mark_vertices[[i]])) mean(degrees[mark_vertices[[i]]]) else 0
      d / (timing$scale + d)
    } else {
      ids <- mark_edges[[i]]
      if (length(ids)) mean(triangle_edge[ids]) else 0
    }
  }
  if (timing$model == "M2") cache$features_after <- vector("list", n)
  for (j in seq_len(n)) {
    mark_edges[[j]] <- which(!is.na(edge_times) & edge_times == times[j])
    new_vertices <- which(!is.na(node_times) & node_times == times[j])
    mark_vertices[[j]] <- unique(c(new_vertices,
                                   as.vector(edge_ends[mark_edges[[j]], , drop = FALSE])))
    for (eid in mark_edges[[j]]) insert_edge(eid)
    cache$birth_features[j] <- feature_now(j)
    if (timing$model == "M2")
      cache$features_after[[j]] <- vapply(seq_len(j), feature_now, numeric(1))
  }
  cache
}

.network_timing_cache <- function(mark_filtration, timing, timing_cache) {
  .network_timing_validate(timing)
  if (is.null(timing_cache)) return(.prepare_network_timing(mark_filtration, timing))
  if (!inherits(timing_cache, "hawkes_timing_cache") ||
      !identical(timing_cache$timing, timing))
    stop("timing_cache must be prepared for the same timing specification")
  timing_cache
}

.network_timing_features <- function(timing_cache, t) {
  j <- sum(timing_cache$times < t)
  if (j == 0L) return(numeric(0))
  if (timing_cache$timing$model == "M2") return(timing_cache$features_after[[j]])
  timing_cache$birth_features[seq_len(j)]
}

#' Evaluate the ground Hawkes intensity with topology feedback
#'
#' Uses only updates strictly before each query time. In M2, features of all
#' such updates are evaluated in the graph at that left limit. The baseline
#' here is the constant params$mu; an externally evaluated background can be
#' substituted by adding `mu_at_t - params$mu` to the result.
#'
#' @param t Finite numeric query times.
#' @param params Named list containing nonnegative mu and K, positive
#'   beta_overall, and optional finite feedback_gamma.
#' @param mark_filtration Network history (or event data.frame for M0).
#' @param timing Specification from [hawkes_timing()].
#' @param timing_cache Optional internal cache prepared from this exact
#'   history and timing specification. The caller must rebuild it if the
#'   history changes; only its timing configuration is checked on reuse.
#' @return Numeric ground intensities, one per query time.
#' @export
hawkes_ground_intensity <- function(t, params, mark_filtration,
                                    timing = hawkes_timing(), timing_cache = NULL) {
  .network_timing_params(params)
  if (!is.numeric(t) || any(!is.finite(t))) stop("t must contain finite numeric query times")
  cache <- .network_timing_cache(mark_filtration, timing, timing_cache)
  vapply(t, function(query) {
    features <- .network_timing_features(cache, query)
    prior <- cache$times[seq_along(features)]
    params$mu + params$K * sum(.network_timing_alpha(features, params, timing) *
                               exp(-params$beta_overall * (query - prior)))
  }, numeric(1), USE.NAMES = FALSE)
}

#' Integrate the ground Hawkes intensity with topology feedback
#'
#' The exponential excitation integral is exact, including any history before
#' the window. M2 splits the window at network updates and uses the features
#' in force on each graph-constant interval. Future marks do not contribute.
#'
#' @param params,mark_filtration,timing,timing_cache See [hawkes_ground_intensity()].
#' @param time_window Two finite numbers giving the observation start and end.
#' @param integral_bg Optional finite nonnegative integral of an externally
#'   specified background over this window; replaces mu times window length.
#' @return The scalar integrated ground intensity over the window.
#' @export
hawkes_ground_compensator <- function(params, time_window, mark_filtration,
                                      timing = hawkes_timing(), timing_cache = NULL,
                                      integral_bg = NULL) {
  .network_timing_params(params)
  if (!is.numeric(time_window) || length(time_window) != 2L ||
      any(!is.finite(time_window)) || time_window[2] < time_window[1])
    stop("time_window must contain finite start and end with end >= start")
  if (!is.null(integral_bg) && (!is.numeric(integral_bg) || length(integral_bg) != 1L ||
      !is.finite(integral_bg) || integral_bg < 0))
    stop("integral_bg must be one finite nonnegative number")
  cache <- .network_timing_cache(mark_filtration, timing, timing_cache)
  start <- time_window[1]
  end <- time_window[2]
  baseline <- if (is.null(integral_bg)) params$mu * (end - start) else integral_bg
  if (end == start || params$K == 0 || !length(cache$times)) return(baseline)
  beta <- params$beta_overall
  if (timing$model != "M2") {
    ids <- which(cache$times < end)
    from <- pmax(start, cache$times[ids])
    weights <- .network_timing_alpha(cache$birth_features[ids], params, timing)
    excitation <- sum(weights * exp(-beta * (from - cache$times[ids])) *
                        (-expm1(-beta * (end - from))) / beta)
  } else {
    cuts <- c(start, cache$times[cache$times > start & cache$times < end], end)
    excitation <- 0
    for (k in seq_len(length(cuts) - 1L)) {
      from <- cuts[k]
      to <- cuts[k + 1L]
      # Integration on (from, to): updates at the left endpoint have occurred.
      j <- sum(cache$times <= from)
      if (j == 0L) next
      ids <- seq_len(j)
      weights <- .network_timing_alpha(cache$features_after[[j]], params, timing)
      excitation <- excitation + sum(weights * exp(-beta * (from - cache$times[ids]))) *
        (-expm1(-beta * (to - from))) / beta
    }
  }
  baseline + params$K * excitation
}
