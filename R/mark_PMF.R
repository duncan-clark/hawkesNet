#' Mark PMF for Barabási–Albert-style (degree-weighted) attachment
#'
#' Probability mass function for the mark at each event: one new node and K edges from that node to existing nodes,
#' where K ~ Poisson(m). So \code{m} is the expected number of edges added per time step. Edges are sampled
#' without replacement with probabilities proportional to (degree * time decay); the likelihood uses the
#' approximation that the probability of the edge set is the product of those degree-based Bernoulli probabilities.
#'
#' @param time Current event time.
#' @param params List with \code{beta_edges}, \code{m} (expected edges per event; default 1), and optionally \code{beta_overall}, \code{K}, \code{mu}.
#' @param mark_filtration Observed network up to \code{time}.
#' @param mark Optional network state at \code{time}; if \code{NULL}, derived from \code{mark_filtration}.
#' @param generate_mark If \code{TRUE}, sample K ~ Poisson(m) and then K distinct edges (default \code{FALSE}).
#' @param new_edge_hash Optional hash of existing edges for fast lookup.
#' @return List with \code{log_mark_density}, \code{log_density_func}, and optionally sampled mark / probabilities.
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
    network::add.vertices(new_net, nv = new_nodes)
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
    log_poisson <- stats::dpois(K_obs, m_val, log = TRUE)
    p_in <- pmax(probs[in_mark], .Machine$double.eps)
    p_out <- pmax(1 - probs[!in_mark], .Machine$double.eps)
    log_mark_density <- log_poisson + sum(log(p_in), na.rm = TRUE) + sum(log(p_out), na.rm = TRUE)
    mark_density <- exp(log_mark_density)
  } else {
    in_mark <- rep(1, length(heads))
    K_obs <- if (length(heads) > 0) sum(in_mark, na.rm = TRUE) else 0
    log_mark_density <- stats::dpois(K_obs, m_val, log = TRUE)
    mark_density <- exp(log_mark_density)
  }
  
  # 1. Define the lightweight log-density function for BA (includes Poisson(m) for K_obs edges)
  log_density_func_light <- function(params) {
    m_p <- if (!is.null(params$m) && is.numeric(params$m) && length(params$m) == 1L && is.finite(params$m) && params$m > 0) params$m else 1
    log_poisson <- stats::dpois(K_obs, m_p, log = TRUE)
    if (!is.finite(log_poisson)) return(-1e10)
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
      edge_part <- sum(log(p_in), na.rm = TRUE) + sum(log(p_out), na.rm = TRUE)
      return(log_poisson + edge_part)
    }
    return(log_poisson)
  }

  # 2. Capture the environment (include K_obs for Poisson term)
  environment(log_density_func_light) <- list2env(
    list(
      heads = heads,
      time = time,
      times = times,
      node_degrees = node_degrees,
      in_mark = in_mark,
      K_obs = K_obs
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
      mark_sample <- last_net
      old_nodes <- last_net %n% 'n'
      new_nodes <- 1
      mark_sample <- network::add.vertices(mark_sample, new_nodes)
      set.vertex.attribute(mark_sample, "time", c((last_net %v% 'time'), rep(time, new_nodes)))
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

      K <- stats::rpois(1, m_val)
      # Join to at most all eligible targets (one node added per event; edges capped by available targets)
      K <- min(K, length(eligible_heads))
      if (K > 0 && length(eligible_heads) > 0) {
        sampled_idx <- sample(length(eligible_heads), size = K, replace = FALSE, prob = probs)
        heads_to_add <- eligible_heads[sampled_idx]
        tails_to_add <- rep(new_node_idx, K)
        network::add.edges(mark_sample, heads_to_add, tails_to_add)
        p_chosen <- pmax(probs[sampled_idx], .Machine$double.eps)
        log_mark_sample_density <- stats::dpois(K, m_val, log = TRUE) + sum(log(p_chosen), na.rm = TRUE)
      } else {
        log_mark_sample_density <- stats::dpois(K, m_val, log = TRUE)
      }
      mark_sample_density <- exp(log_mark_sample_density)
    } else {
      if (is.null(last_net) || (last_net %n% 'n') < 2) {
        if (is.null(last_net)) {
          mark_sample <- network::network(matrix(1), directed = FALSE)
          set.vertex.attribute(mark_sample, "time", time)
        } else {
          mark_sample <- last_net
        }
        times <- mark_sample %v% 'time'
        mark_sample <- network::add.vertices(mark_sample, 1)
        if (mark_sample %n% 'n' == 2) {
          network::add.edges(mark_sample, 2, 1)
        }
        set.vertex.attribute(mark_sample, "time", c(times, time))
        log_mark_sample_density <- stats::dpois(1, m_val, log = TRUE)
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
#' @param vertex_categorical Optional named list in \code{params}: for each discrete vertex attribute, a named
#'   numeric vector of multinomial proportions (e.g. \code{list(gender = c(male=0.4, female=0.4, unknown=0.2))}).
#'   Used for likelihood (observed new-node attribute values) and simulation (sampling new-node attributes).
#'   Values are normalized to sum to 1 per attribute; any positive values are valid.
#' @seealso \code{\link[network]{network}}, \code{\link[network]{add.vertices}}, \code{\link[ernm]{as.BinaryNet}}
#' @rdname PMF_mark_CS
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
  if (!is.null(mark) && attr_name %in% network::list.vertex.attributes(mark)) {
    return(sort(unique(mark %v% attr_name)))
  }
  NULL
}

#' Expected parameter names for PMF_mark_BA
#'
#' @return List with \code{required} (character vector: \code{beta_edges}, \code{m}).
#' @export
expected_params_PMF_mark_BA <- function() {
  list(required = c("beta_edges", "m"))
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
  nv <- network::network.size(net)
  if (!is.finite(nv) || is.na(nv)) nv <- 0
  if (nv < 4) {
    network::add.vertices(net, 4 - nv)
    t0 <- if (nv > 0) (net %v% "time")[1] else times[1]
    network::set.vertex.attribute(net, "time", c(net %v% "time", rep(t0, 4 - nv)))
  }
  if ("na" %in% network::list.vertex.attributes(net)) {
    delete.vertex.attribute(net, "na")
  }
  # If formula uses nodeMatch('gender'), ensure network has 'gender' so createCppModel/setNetwork succeed
  if (grepl("nodeMatch\\s*\\(\\s*['\"]gender['\"]", formula_RHS) && nv > 0L) {
    if (!"gender" %in% network::list.vertex.attributes(net)) {
      network::set.vertex.attribute(net, "gender", rep("unknown", nv))
    } else {
      attr_vals <- network::get.vertex.attribute(net, "gender")
      if (length(attr_vals) < nv || any(is.na(attr_vals)) || any(attr_vals == "")) {
        attr_vals <- if (length(attr_vals) < nv) c(attr_vals, rep("unknown", nv - length(attr_vals))) else attr_vals
        attr_vals[is.na(attr_vals) | attr_vals == ""] <- "unknown"
        network::set.vertex.attribute(net, "gender", attr_vals)
      }
    }
  }
  CS_params_names <- NULL
  tryCatch({
    model <- createCppModel(as.formula(paste("net ~ ", formula_RHS)))
    model$setNetwork(as.BinaryNet(net))
    model$calculate()
    stats <- model$statistics()
    CS_params_length <- length(stats)
    CS_params_names <- names(stats)
  }, error = function(e) NULL)
  list(required = required, CS_params_length = CS_params_length, CS_params_names = CS_params_names)
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
#' @return List with \code{tails} and \code{heads} integer vectors.
#' @noRd
get_truncated_candidates <- function(net, new_nodes, old_nodes, truncation, mark_decay) {
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
    poss_tails <- (old_nodes - truncation):(new_nodes)
    poss_tails <- poss_tails[poss_tails > 0]
    poss_heads <- (new_nodes - truncation - 1):new_nodes
    poss_heads <- poss_heads[poss_heads > 0]
    poss_edges <- expand.grid(poss_tails, poss_heads)
    poss_edges <- poss_edges[poss_edges[, 1] > poss_edges[, 2], , drop = FALSE]
    tails <- poss_edges[, 1]
    heads <- poss_edges[, 2]
  }
  list(tails = tails, heads = heads)
}

#' Mark probability mass function using Change Statistics (CS/ERNM)
#'
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
                        model = NULL,
                        max_node_time = NULL,
                        ...
){
  eps <- 1e-10  # used for probability clamping and safe log (CS safety)
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
    network::add.vertices(new_net,new_nodes-old_nodes)
    set.vertex.attribute(new_net,"time",c(last_net %v% 'time',rep(time,new_nodes)))
  }else{
    last_net <- NULL
    new_net <- network::network(matrix(1),directed = F)
    delete.vertex.attribute(new_net,'na')
    set.vertex.attribute(new_net,"time",time)
    old_nodes <- 0
    new_nodes <- 1
  }
  # get the possible edges for the given truncation:
  cands <- get_truncated_candidates(new_net, new_nodes, old_nodes, truncation, mark_decay)
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
      # if new net has less than 4 nodes add some:
      if(new_net %n% 'n' < 4){
        old_new_net <- new_net
        new_net <- network::add.vertices(new_net,4 - (new_net %n% 'n'))
      }else{
        old_new_net <- new_net
      }
      
      if(max(tails)>new_net %n% 'n'){
        stop("accidently adding a edge into the network that doesn't have that node yet")
      }
      
      # Copy discrete vertex attributes from mark so change stats (e.g. nodeMix) are correct
      vcat <- params$vertex_categorical
      if (!is.null(vcat) && is.list(vcat)) {
        nv <- network::network.size(new_net)
        for (attr_name in names(vcat)) {
          if (attr_name %in% network::list.vertex.attributes(mark)) {
            g <- mark %v% attr_name
            nm <- length(g)
            if (nv <= nm) {
              network::set.vertex.attribute(new_net, attr_name, g[seq_len(nv)])
            } else {
              levs <- vertex_categorical_level_names(params, attr_name, mark)
              if (is.null(levs) || length(levs) == 0) levs <- "unknown"
              network::set.vertex.attribute(new_net, attr_name, c(g, rep(levs[1L], nv - nm)))
            }
          }
        }
      }
      
      # delete NAs to prevent C++ using them
      delete.vertex.attribute(new_net,'na')
      if(is.null(model)){
        model <- createCppModel(as.formula(paste("new_net ~ ",formula_RHS)))
      }else{
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
      probs <- plogis(eta)
      # --- Safety: sanitize probs after logistic (suggestions 1 & 9) ---
      if (any(!is.finite(probs))) {
        warning("PMF_mark_CS: NA/NaN/Inf in edge probs after logistic; replacing with 0 before clamp.")
        probs[!is.finite(probs)] <- 0
      }
      probs <- pmin(pmax(probs, eps), 1 - eps)
      if (length(probs) > 0L && all(probs <= eps)) {
        warning("PMF_mark_CS: all edge probs effectively zero after logistic; using uniform probs.")
        probs[] <- 1 / length(probs)
      }
      # use either node times or last node activity:
      if(mark_decay == 'activity'){
        node_times <- get_latest_times(new_net)
      }
      if(mark_decay == 'node_entrance'){
        node_times <- new_net %v% 'time'
      }
      diffs <- time - node_times[heads]
      # --- Safety: sanitize diffs/factor before multiplying probs (suggestion 2 & 9) ---
      if (any(!is.finite(diffs))) {
        warning("PMF_mark_CS: non-finite time diffs in density path; replacing with 0.")
        diffs[!is.finite(diffs)] <- 0
      }
      factor <- exp(-params$beta_edges*(diffs))
      if (any(!is.finite(factor))) {
        warning("PMF_mark_CS: non-finite decay factor in density path; replacing with 1.")
        factor[!is.finite(factor)] <- 1
      }
      probs <- probs * factor
      # --- Safety: sanitize probs after factor (suggestion 3 & 9) ---
      if (any(!is.finite(probs))) {
        warning("PMF_mark_CS: NA/NaN/Inf in edge probs after decay factor; replacing with 0 before clamp.")
        probs[!is.finite(probs)] <- 0
      }
      probs <- pmin(pmax(probs, eps), 1 - eps)
      if (length(probs) > 0L && all(probs <= eps)) {
        warning("PMF_mark_CS: all edge probs effectively zero after decay; using uniform probs.")
        probs[] <- 1 / length(probs)
      }

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
        dval <- stats::dpois(new_nodes-old_nodes, params$node_lambda)
        if (!is.finite(dval) || dval <= 0) {
          warning("PMF_mark_CS: degenerate node count density (dpois=0 or non-finite); using large negative log-density.")
          node_dens <- -1e10
        } else {
          node_dens <- log(dval)
        }
      }
      # Log multinomial contribution for discrete vertex attributes on new nodes
      observed_categorical <- list()
      vcat <- params$vertex_categorical
      if (!is.null(vcat) && is.list(vcat) && (new_nodes - old_nodes) > 0) {
        for (attr_name in names(vcat)) {
          if (attr_name %in% network::list.vertex.attributes(mark)) {
            observed_categorical[[attr_name]] <- (mark %v% attr_name)[(old_nodes + 1):new_nodes]
            level_names <- vertex_categorical_level_names(params, attr_name, mark)
            p <- expand_vertex_categorical_probs(vcat[[attr_name]], level_names, eps = eps)
            if (!is.null(p)) {
              idx <- match(observed_categorical[[attr_name]], names(p))
              idx[is.na(idx)] <- match("unknown", names(p))
              idx[is.na(idx)] <- 1L
              node_dens <- node_dens + sum(log(pmax(p[idx], eps)))
            }
          }
        }
      }
      # --- Safety: safe log with clamped probs (suggestion 4) ---
      p_in <- pmax(probs[in_mark], eps, na.rm = TRUE)
      p_out <- pmax(1 - probs[!in_mark], eps, na.rm = TRUE)
      if (any(!is.finite(p_in)) || any(!is.finite(p_out))) {
        warning("PMF_mark_CS: non-finite probs in log_mark_density; using epsilon for log.")
      }
      log_mark_density <- sum(log(p_in), na.rm = TRUE) + sum(log(p_out), na.rm = TRUE) + node_dens
      mark_density <- exp(log_mark_density)
    }
  }else{
    mark_density <- 1
    mark_density_normalized <- NULL
    log_mark_density <- 0
  }
  
  # if node_dens doesn't exist set it to 0
  if(!exists("node_dens")){
    node_dens <- 0
  }
  
  log_density_func_light <- function(params) {
    eta <- change_stats %*% params$CS_params
    p   <- stats::plogis(eta)
    log_edge_part <- sum(log(p[in_mark])) + sum(log1p(-p[!in_mark]))
    log_edge_part + node_dens
  }

  
  environment(log_density_func_light) <- list2env(
    list(
      change_stats = change_stats,
      in_mark      = in_mark,
      node_dens    = node_dens,
      plogis = stats::plogis
    ),
    parent = baseenv()
  )
  
  # define the function (will be rebound to a minimal env right after)
  log_density_func_light <- function(params) {
    # Everything it needs will come from its environment:
    # change_stats, in_mark, diffs, new_nodes, old_nodes, time, max_node_time, degenerate_edges
    
    if (degenerate_edges) {
      return(0)
    }
    
    eta    <- as.vector(change_stats %*% params$CS_params)
    p_base <- stats::plogis(eta)
    
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
    # --- Safety: handle dpois=0 or non-finite in closure (suggestion 8) ---
    node_dens <- if (!is.null(max_node_time) && time > max_node_time) {
      0
    } else {
      dval <- stats::dpois(new_nodes - old_nodes, params$node_lambda)
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
  
  # --- Safety: sanitize diffs for closure (suggestion 10) ---
  diffs_for_closure <- if (exists("diffs", inherits = FALSE)) {
    d <- diffs
    if (any(!is.finite(d))) {
      warning("PMF_mark_CS: non-finite diffs passed to log_density_func closure; replacing with 0.")
      d[!is.finite(d)] <- 0
    }
    d
  } else numeric(0)
  
  # Pre-compute level names per attribute for the closure
  obs_cat <- if (exists("observed_categorical", inherits = FALSE)) observed_categorical else list()
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
      observed_categorical            = obs_cat,
      level_names_by_attr             = level_names_by_attr,
      expand_vertex_categorical_probs = expand_vertex_categorical_probs
    ),
    parent = baseenv()
  )
  
  density_func_light <- function(params) exp(log_density_func_light(params))
  environment(density_func_light) <- list2env(
    list(log_density_func_light = log_density_func_light),
    parent = baseenv()
  )

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
        old_new_net <- mark_sample
        new_nodes <- 4 - (mark_sample %n% 'n')
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
      mark_sample <- network::add.vertices(mark_sample, new_nodes)
      # if(mark_sample %n% 'n' > 4){
      #   browser()
      # }
      set.vertex.attribute(mark_sample,"time",c((last_net %v% 'time'),rep(time,new_nodes)))
      # Set discrete vertex attributes for new nodes (sample from vertex_categorical; n-1 params)
      vcat <- params$vertex_categorical
      if (!is.null(vcat) && is.list(vcat)) {
        for (attr_name in names(vcat)) {
          levs <- vertex_categorical_level_names(params, attr_name, mark_sample)
          p <- expand_vertex_categorical_probs(vcat[[attr_name]], levs, eps = eps)
          if (is.null(p)) next
          levs <- names(p)
          if (any(!is.finite(p)) || sum(p) <= 0) p <- rep(1 / length(levs), length(levs)); names(p) <- levs
          existing <- if (attr_name %in% network::list.vertex.attributes(last_net)) last_net %v% attr_name else rep(levs[1L], old_nodes)
          if (new_nodes > 0) {
            sampled <- sample(levs, size = new_nodes, replace = TRUE, prob = p)
            network::set.vertex.attribute(mark_sample, attr_name, c(existing, sampled))
          } else {
            network::set.vertex.attribute(mark_sample, attr_name, existing)
          }
        }
      }
      new_size <- mark_sample %n% 'n'
      
      # get new poss edges (truncation based on mark_decay)
      cands <- get_truncated_candidates(mark_sample, new_size, old_nodes, truncation, mark_decay)
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
      
      # CRITICAL: Ensure ALL nodes have required vertex attributes before createCppModel
      vcat <- params$vertex_categorical
      if (!is.null(vcat) && is.list(vcat)) {
        nv <- network::network.size(mark_sample)
        for (attr_name in names(vcat)) {
          if (!attr_name %in% network::list.vertex.attributes(mark_sample)) {
            levs <- if (!is.null(params$vertex_categorical_levels) && attr_name %in% names(params$vertex_categorical_levels)) {
              params$vertex_categorical_levels[[attr_name]]
            } else {
              c("unknown")
            }
            if (is.null(levs) || length(levs) == 0) levs <- c("unknown")
            network::set.vertex.attribute(mark_sample, attr_name, rep(levs[1L], nv))
          } else {
            attr_vals <- mark_sample %v% attr_name
            if (length(attr_vals) < nv || any(is.na(attr_vals)) || any(attr_vals == "")) {
              levs <- if (!is.null(params$vertex_categorical_levels) && attr_name %in% names(params$vertex_categorical_levels)) {
                params$vertex_categorical_levels[[attr_name]]
              } else {
                unique_vals <- unique(attr_vals[!is.na(attr_vals) & attr_vals != ""])
                if (length(unique_vals) > 0) sort(unique_vals) else c("unknown")
              }
              if (is.null(levs) || length(levs) == 0) levs <- c("unknown")
              if (length(attr_vals) < nv) {
                attr_vals <- c(attr_vals, rep(levs[1L], nv - length(attr_vals)))
              }
              attr_vals[is.na(attr_vals) | attr_vals == ""] <- levs[1L]
              network::set.vertex.attribute(mark_sample, attr_name, attr_vals)
            }
          }
        }
      }
      if ("na" %in% network::list.vertex.attributes(mark_sample)) {
        delete.vertex.attribute(mark_sample, "na")
      }
      
      # Reuse ERNM model per formula (avoids createCppModel every event; big speedup for nodeMatch)
      cache <- get0(".ernm_model_cache", envir = asNamespace("hawkesNet"), inherits = FALSE)
      if (is.null(cache)) {
        cache <- new.env()
        assign(".ernm_model_cache", cache, envir = asNamespace("hawkesNet"))
      }
      key <- formula_RHS
      if (is.null(cache[[key]])) {
        g0 <- network::network.initialize(0L, directed = FALSE)
        cache[[key]] <- createCppModel(as.formula(paste("g0 ~ ", formula_RHS)))
        cache[[key]]$setNetwork(as.BinaryNet(g0))
      }
      model <- cache[[key]]
      model$setNetwork(as.BinaryNet(mark_sample))
      model$calculate()
      change_stats <- model$computeChangeStats(tails, heads)
      eta <- as.vector(change_stats %*% params$CS_params)
      probs <- plogis(eta)
      # --- Safety: sanitize probs after logistic in generate_mark (suggestions 1 & 9) ---
      if (any(!is.finite(probs))) {
        warning("PMF_mark_CS (generate_mark): NA/NaN/Inf in edge probs after logistic; replacing with 0 before clamp.")
        probs[!is.finite(probs)] <- 0
      }
      probs <- pmin(pmax(probs, eps), 1 - eps)
      if (length(probs) > 0L && all(probs <= eps)) {
        warning("PMF_mark_CS (generate_mark): all edge probs effectively zero after logistic; using uniform probs.")
        probs[] <- 1 / length(probs)
      }

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
          mark_sample <- network::add.vertices(mark_sample, 1)
          network::set.vertex.attribute(mark_sample, "time", c(mark_sample %v% "time", time))
          vcat <- params$vertex_categorical
          if (!is.null(vcat) && is.list(vcat)) {
            for (attr_name in names(vcat)) {
              levs <- vertex_categorical_level_names(params, attr_name, mark_sample)
              p <- expand_vertex_categorical_probs(vcat[[attr_name]], levs, eps = eps)
              if (is.null(p)) next
              levs <- names(p)
              if (any(!is.finite(p)) || sum(p) <= 0) p <- rep(1 / length(levs), length(levs))
              names(p) <- levs
              nv <- network::network.size(mark_sample)
              existing <- if (attr_name %in% network::list.vertex.attributes(mark_sample)) (mark_sample %v% attr_name)[seq_len(nv - 1)] else rep(levs[1L], nv - 1)
              sampled_one <- sample(levs, size = 1L, replace = TRUE, prob = p)
              network::set.vertex.attribute(mark_sample, attr_name, c(existing, sampled_one))
            }
          }
        }
        mark_sample_density <- 1
        log_mark_sample_density <- 0
      } else {
      if(mark_decay == 'activity'){
        node_times <- get_latest_times(mark_sample)
      }
      if(mark_decay == 'node_entrance'){
        node_times <- mark_sample %v% 'time'
      }
      diffs <- sapply(seq_along(tails), function(i) {
        node_times[tails[i]] - node_times[heads[i]]
      })
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
      probs <- factor * probs
      # --- Safety: sanitize probs after decay in generate_mark (suggestion 3 & 9) ---
      if (any(!is.finite(probs))) {
        warning("PMF_mark_CS (generate_mark): NA/NaN/Inf in edge probs after decay; replacing with 0 before clamp.")
        probs[!is.finite(probs)] <- 0
      }
      probs <- pmin(pmax(probs, eps), 1 - eps)
      if (length(probs) > 0L && all(probs <= eps)) {
        warning("PMF_mark_CS (generate_mark): all edge probs effectively zero after decay; using uniform probs.")
        probs[] <- 1 / length(probs)
      }

      add <- runif(length(probs)) < probs
      # --- Safety: no NA in add before add.edges (suggestion 5 & 9) ---
      if (any(is.na(add))) {
        warning("PMF_mark_CS (generate_mark): NA in edge add vector; treating as FALSE (do not add edge).")
        add[is.na(add)] <- FALSE
      }
      add.edges(mark_sample,
                heads[add],
                tails[add]
      )
      set.edge.attribute(mark_sample,"time",c(mark_sample %e% 'time',rep(time,sum(add))))
      # --- Safety: safe log and dpois for sample density (suggestion 6 & 8) ---
      dpois_val <- stats::dpois(new_nodes-old_nodes, params$node_lambda)
      if (!is.finite(dpois_val) || dpois_val <= 0) {
        warning("PMF_mark_CS (generate_mark): degenerate dpois for node count; using small positive value for density.")
        dpois_val <- 1e-300
      }
      p_add <- pmax(probs[add], eps, na.rm = TRUE)
      p_not <- pmax(1 - probs[!add], eps, na.rm = TRUE)
      if (any(!is.finite(p_add)) || any(!is.finite(p_not))) {
        warning("PMF_mark_CS (generate_mark): non-finite probs in log_mark_sample_density; using epsilon for log.")
      }
      log_multinomial_sample <- 0
      vcat <- params$vertex_categorical
      if (!is.null(vcat) && is.list(vcat) && (new_size - old_nodes) > 0) {
        for (attr_name in names(vcat)) {
          if (attr_name %in% network::list.vertex.attributes(mark_sample)) {
            levs <- vertex_categorical_level_names(params, attr_name, mark_sample)
            p <- expand_vertex_categorical_probs(vcat[[attr_name]], levs, eps = eps)
            if (!is.null(p)) {
              obs <- (mark_sample %v% attr_name)[(old_nodes + 1):new_size]
              idx <- match(obs, names(p))
              idx[is.na(idx)] <- match("unknown", names(p))
              idx[is.na(idx)] <- 1L
              log_multinomial_sample <- log_multinomial_sample + sum(log(pmax(p[idx], eps)))
            }
          }
        }
      }
      mark_sample_density <- prod(p_add) * prod(p_not) * dpois_val * exp(log_multinomial_sample)
      log_mark_sample_density <- sum(log(p_add), na.rm = TRUE) +
                                 sum(log(p_not), na.rm = TRUE) +
                                 log(dpois_val) + log_multinomial_sample
      }
      }else{
        if(is.null(last_net)){
          mark_sample <- network::network(matrix(1),directed = F)
          set.vertex.attribute(mark_sample,"time",time)
        }else{
          mark_sample <- last_net
        }
        times <- mark_sample %v% 'time'
        mark_sample <- network::add.vertices(mark_sample,1)
        set.vertex.attribute(mark_sample,
                             "time",
                             c(times,time))
        vcat <- params$vertex_categorical
        if (!is.null(vcat) && is.list(vcat)) {
          for (attr_name in names(vcat)) {
            levs <- vertex_categorical_level_names(params, attr_name, mark)
            p <- expand_vertex_categorical_probs(vcat[[attr_name]], levs, eps = eps)
            if (is.null(p)) next
            levs <- names(p)
            if (any(!is.finite(p)) || sum(p) <= 0) p <- rep(1 / length(levs), length(levs)); names(p) <- levs
            existing <- if (attr_name %in% network::list.vertex.attributes(mark_sample)) (mark_sample %v% attr_name)[seq_len(length(times))] else rep(levs[1L], length(times))
            sampled_one <- sample(levs, size = 1L, replace = TRUE, prob = p)
            network::set.vertex.attribute(mark_sample, attr_name, c(existing, sampled_one))
          }
        }
        mark_sample_density <- 1
        log_mark_sample_density <- 0
      }
    }else{
      mark_sample <- new_net
      mark_sample_density <- 1
      log_mark_sample_density <- 0
    }

  return(list(
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
  ))
}
