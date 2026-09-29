# Exact whole-update CS probabilities. ERNM toggles below are calculations,
# never an ordering of the simultaneous edges in a scientific mark.

.cs_change_model <- function(net, formula_RHS) {
  baseline <- network::network.copy(net)
  if ("na" %in% network::list.vertex.attributes(baseline)) {
    network::delete.vertex.attribute(baseline, "na")
  }
  model <- ernm::createCppModel(stats::as.formula(
    paste("baseline ~", formula_RHS), env = environment()))
  model$calculate()
  model
}

#' Calculate the joint change statistics of simultaneous edge additions
#'
#' Uses sequential ERNM toggles to evaluate the exact change caused by the whole
#' unordered edge set. The calculation does not define an edge-order likelihood.
#' @param net Baseline network, already containing any nodes added by the event.
#' @param edges Two-column matrix of distinct, absent edges to add.
#' @param formula_RHS ERNM statistic formula right-hand side.
#' @return List with total change, sequential increments, and endpoint statistics.
#' The input network is not modified. Changes due solely to adding vertices are
#' outside this edge calculation.
#' @export
multi_edge_change_stats <- function(net, edges, formula_RHS) {
  edges <- as.matrix(edges)
  if (length(edges) == 0L) edges <- matrix(integer(0), ncol = 2L)
  if (ncol(edges) != 2L || anyNA(edges) || any(!is.finite(edges)) ||
      any(edges != as.integer(edges)) || any(edges < 1L) ||
      any(edges > network::network.size(net)) ||
      any(edges[, 1L] == edges[, 2L])) {
    stop("edges must contain valid, distinct, non-loop vertex pairs")
  }
  storage.mode(edges) <- "integer"
  canonical <- if (network::is.directed(net)) edges else
    cbind(pmin(edges[, 1L], edges[, 2L]), pmax(edges[, 1L], edges[, 2L]))
  if (anyDuplicated(data.frame(canonical))) stop("edges contains duplicate dyads")
  model <- .cs_change_model(net, formula_RHS)
  binary <- model$getNetwork()
  if (nrow(edges) && any(binary$getDyads(edges[, 1L], edges[, 2L]))) {
    stop("multi_edge_change_stats accepts edge additions only")
  }
  before <- model$statistics()
  increments <- matrix(0, nrow(edges), length(before),
                       dimnames = list(NULL, names(before)))
  previous <- before
  for (j in seq_len(nrow(edges))) {
    model$dyadUpdate(edges[j, 1L], edges[j, 2L])
    binary$setDyads(edges[j, 1L], edges[j, 2L], TRUE)
    current <- model$statistics()
    increments[j, ] <- current - previous
    previous <- current
  }
  list(total = previous - before, increments = increments,
       before = before, after = previous)
}

.cs_state_table <- function(net, tails, heads, formula_RHS, max_candidates) {
  d <- length(tails)
  if (length(max_candidates) != 1L || !is.finite(max_candidates) ||
      max_candidates < 0L || max_candidates > 20L ||
      max_candidates != as.integer(max_candidates)) {
    stop("max_candidates must be an integer between 0 and 20")
  }
  if (d > max_candidates) {
    stop(sprintf(paste0("Exact whole-update CS requires 2^D subset states; D=%d ",
                        "exceeds max_candidates=%d. Reduce the candidate window. ",
                        "No independent-edge approximation is substituted."),
                 d, max_candidates))
  }
  model <- .cs_change_model(net, formula_RHS)
  binary <- model$getNetwork()
  before <- model$statistics()
  nstates <- as.integer(2^d)
  states <- matrix(FALSE, nstates, d)
  changes <- matrix(0, nstates, length(before),
                    dimnames = list(NULL, names(before)))
  if (d) {
    ids <- seq.int(0L, nstates - 1L)
    for (j in seq_len(d)) states[, j] <- bitwAnd(ids, bitwShiftL(1L, j - 1L)) != 0L
    # Gray-code traversal changes one dyad at a time. Store each total under its
    # unordered subset bit mask, independent of the traversal used to obtain it.
    previous <- 0L
    for (step in seq_len(nstates - 1L)) {
      mask <- bitwXor(step, bitwShiftR(step, 1L))
      changed <- bitwXor(mask, previous)
      j <- which(as.logical(intToBits(changed)[seq_len(d)]))
      model$dyadUpdate(tails[j], heads[j])
      binary$setDyads(tails[j], heads[j], states[mask + 1L, j])
      changes[mask + 1L, ] <- model$statistics() - before
      previous <- mask
    }
  }
  list(states = states, changes = changes, tails = tails, heads = heads,
       baseline_stats = before)
}

.cs_logsumexp <- function(x) {
  peak <- max(c(-Inf, x))
  if (!is.finite(peak)) return(peak)
  peak + log(sum(exp(x - peak)))
}

# Sum positive update masses directly: subtracting a rounded empty probability
# from one loses valid support when edge additions are very unlikely.
.cs_nonempty_log <- function(lambda, zero_log_probs) {
  .cs_logsumexp(c(log(-expm1(-lambda)),
                 -lambda + .cs_logsumexp(zero_log_probs[-1L])))
}

.cs_table_distribution <- function(tab, params) {
  theta <- as.numeric(params$CS_params)
  if (length(theta) != ncol(tab$changes) || any(!is.finite(theta))) {
    stop("CS_params must match the finite whole-update statistic vector")
  }
  m <- params$m
  tau <- params$beta_edges
  if (length(m) != 1L || !is.finite(m) || m < 0 ||
      length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("m and beta_edges must be finite nonnegative scalars")
  }
  d <- ncol(tab$states)
  if (!d) return(list(log_probs = 0, probs = 1, reference_probs = numeric(0), log_normalizer = 0))
  log_weights <- -tau * tab$ages
  p <- exp(log_weights - max(log_weights))
  p <- p / sum(p)
  x <- m * p
  log_selected <- x + log(-expm1(-x))
  zero <- !is.finite(log_selected)
  log_selected[zero] <- 0
  log_ref <- -m + as.vector(tab$states %*% log_selected)
  if (any(zero)) log_ref[rowSums(tab$states[, zero, drop = FALSE]) > 0L] <- -Inf
  log_score <- log_ref + as.vector(tab$changes %*% theta)
  peak <- max(log_score)
  log_z <- peak + log(sum(exp(log_score - peak)))
  log_probs <- log_score - log_z
  list(log_probs = log_probs, probs = exp(log_probs),
       reference_probs = -expm1(-x), log_normalizer = log_z)
}

.cs_joint_context <- function(old, proposed, time, params, formula_RHS,
                              truncation, mark_decay, growth_only, max_candidates) {
  old_n <- network::network.size(old)
  n <- network::network.size(proposed)
  baseline <- network::network.copy(old)
  if (n > old_n) {
    baseline <- network::add.vertices(baseline, n - old_n)
    network::set.vertex.attribute(baseline, "time", c(
      network::get.vertex.attribute(old, "time"), rep(time, n - old_n)))
    for (attr in names(params$vertex_categorical)) {
      values <- network::get.vertex.attribute(proposed, attr)
      network::set.vertex.attribute(baseline, attr, values)
    }
  }
  baseline <- strip_vertex_attrs_for_ernm(baseline, formula_RHS, params)
  cands <- get_truncated_candidates(baseline, n, old_n, truncation,
                                    mark_decay, growth_only)
  existing <- vapply(seq_along(cands$tails), function(i) {
    length(network::get.edgeIDs(baseline, cands$tails[i], cands$heads[i])) > 0L
  }, logical(1))
  tails <- as.integer(cands$tails[!existing])
  heads <- as.integer(cands$heads[!existing])
  tab <- .cs_state_table(baseline, tails, heads, formula_RHS, max_candidates)
  # Conditional on the births, their statistic contribution is constant and
  # cancels from the edge normalizer. Retain it for the complete-mark report.
  tab$birth_change <- if (n > old_n) {
    before_births <- strip_vertex_attrs_for_ernm(old, formula_RHS, params)
    tab$baseline_stats - .cs_change_model(before_births, formula_RHS)$statistics()
  } else rep(0, length(tab$baseline_stats))
  times <- if (identical(mark_decay, "activity")) get_latest_times(baseline) else
    network::get.vertex.attribute(baseline, "time")
  ages <- if (identical(mark_decay, "activity")) {
    time - pmax(times[tails], times[heads])
  } else time - times[heads]
  if (any(!is.finite(ages)) || any(ages < 0)) stop("Candidate endpoint ages must be finite and nonnegative")
  tab$ages <- ages
  tab
}

.cs_edge_keys <- function(net) {
  e <- network::as.edgelist(net)
  if (!nrow(e)) return(character(0))
  paste(pmin(e[, 1L], e[, 2L]), pmax(e[, 1L], e[, 2L]), sep = "-")
}

.cs_attribute_log <- function(params, observed, levels) {
  ans <- 0
  for (attr in names(observed)) {
    p <- expand_vertex_categorical_probs(params$vertex_categorical[[attr]], levels[[attr]])
    idx <- match(observed[[attr]], names(p))
    if (anyNA(idx)) return(-Inf)
    ans <- ans + sum(log(p[idx]))
  }
  ans
}

.pmf_cs_joint <- function(time, params, mark_filtration, mark,
                          generate_mark, formula_RHS, truncation, mark_decay,
                          growth_only, max_node_time, max_candidates,
                          condition_nonempty) {
  if (is.null(truncation)) truncation <- 4L
  if (length(truncation) != 1L || !is.finite(truncation) || truncation < 1L ||
      truncation != as.integer(truncation)) stop("truncation must be a positive integer node-window size")
  mark_decay <- match.arg(mark_decay, c("node_entrance", "activity"))
  if (is.null(formula_RHS)) stop("formula_RHS is required for whole-update CS")
  if (is.null(max_node_time)) max_node_time <- Inf
  empty <- network::network.initialize(0, directed = FALSE)
  network::set.vertex.attribute(empty, "time", numeric(0))
  if (is.null(mark_filtration)) mark_filtration <- empty
  old <- if (generate_mark && !is.null(mark)) network::network.copy(mark) else
    filtration_to_net(mark_filtration, time, equals = FALSE)
  if (is.null(old)) old <- empty
  if (network::is.directed(old)) stop("Whole-update CS currently supports undirected edge additions")
  old_n <- network::network.size(old)
  zero_tab <- .cs_joint_context(old, old, time, params, formula_RHS,
                                truncation, mark_decay, growth_only, max_candidates)
  zero_dist <- .cs_table_distribution(zero_tab, params)
  lambda <- if (time > max_node_time) 0 else params$node_lambda
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) stop("node_lambda must be finite and nonnegative")
  log_nonempty <- .cs_nonempty_log(lambda, zero_dist$log_probs)

  if (generate_mark) {
    if (condition_nonempty) {
      if (!is.finite(log_nonempty)) stop("No nonempty CS update is possible from this state")
      log_mass_zero_nodes <- -lambda + .cs_logsumexp(zero_dist$log_probs[-1L])
      choose_zero <- stats::runif(1) < exp(log_mass_zero_nodes - log_nonempty)
      b <- if (choose_zero) 0L else as.integer(stats::qpois(
        max(stats::runif(1), .Machine$double.xmin) * (-expm1(-lambda)),
        lambda, lower.tail = FALSE))
    } else b <- stats::rpois(1, lambda)
    proposed <- network::network.copy(old)
    if (b) {
      proposed <- network::add.vertices(proposed, b)
      network::set.vertex.attribute(proposed, "time", c(
        network::get.vertex.attribute(old, "time"), rep(time, b)))
      proposed <- sample_vertex_attrs(params, old, proposed, old_n, b)
    }
    tab <- if (!b) zero_tab else .cs_joint_context(
      old, proposed, time, params, formula_RHS, truncation, mark_decay,
      growth_only, max_candidates)
    dist <- if (!b) zero_dist else .cs_table_distribution(tab, params)
    draw_log_probs <- dist$log_probs
    if (condition_nonempty && b == 0L) draw_log_probs[1L] <- -Inf
    draw_probs <- exp(draw_log_probs - max(draw_log_probs))
    selected_state <- sample.int(length(draw_probs), 1L, prob = draw_probs)
    selected <- which(tab$states[selected_state, ])
    if (length(selected)) {
      proposed <- network::add.edges(proposed, tab$tails[selected], tab$heads[selected])
      new_ids <- unlist(lapply(selected, function(j) network::get.edgeIDs(
        proposed, tab$tails[j], tab$heads[j])), use.names = FALSE)
      network::set.edge.attribute(proposed, "time", rep(time, length(new_ids)), e = new_ids)
    }
    mark <- proposed
  } else {
    if (is.null(mark)) mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
    b <- network::network.size(mark) - old_n
    if (b < 0L || network::is.directed(mark)) stop("CS marks must preserve old vertices and be undirected")
    tab <- if (!b) zero_tab else .cs_joint_context(
      old, mark, time, params, formula_RHS, truncation, mark_decay,
      growth_only, max_candidates)
    dist <- if (!b) zero_dist else .cs_table_distribution(tab, params)
    old_keys <- .cs_edge_keys(old)
    mark_keys <- .cs_edge_keys(mark)
    added <- setdiff(mark_keys, old_keys)
    candidate_keys <- paste(pmin(tab$tails, tab$heads), pmax(tab$tails, tab$heads), sep = "-")
    valid <- all(old_keys %in% mark_keys) && all(added %in% candidate_keys)
    selected <- which(candidate_keys %in% added)
    selected_state <- if (valid) 1L + sum(2^(selected - 1L)) else NA_integer_
  }
  valid <- !is.na(selected_state) && (!condition_nonempty || b > 0L || length(selected) > 0L)
  if (b > 0L) for (attr in names(params$vertex_categorical)) {
    if (!(attr %in% network::list.vertex.attributes(mark)) ||
        anyNA(network::get.vertex.attribute(mark, attr)[seq.int(old_n + 1L, old_n + b)])) {
      stop(sprintf("New vertices must include nonmissing categorical attribute '%s'", attr))
    }
  }
  cat_res <- log_categorical_density(params, mark, old_n, old_n + b)
  observed <- cat_res$observed
  levels <- lapply(names(observed), function(attr) vertex_categorical_level_names(params, attr, mark))
  names(levels) <- names(observed)
  # The cache contains numerical state tables and observed attributes, not an
  # ERNM external pointer or an arbitrary within-mark edge order.
  # Use a minimal environment so likelihood caches do not retain the entire
  # working network and all transient objects for every historical event.
  cache <- list2env(list(table_local = tab, zero_local = zero_tab,
    state_local = selected_state, b_local = b, valid_local = valid,
    no_births = time > max_node_time, observed_local = observed,
    levels_local = levels, nonempty_local = condition_nonempty,
    .cs_table_distribution = .cs_table_distribution,
    .cs_attribute_log = .cs_attribute_log,
    .cs_nonempty_log = .cs_nonempty_log), parent = baseenv())
  log_density_func <- eval(quote(function(params) {
      if (!valid_local) return(-Inf)
      lam <- if (no_births) 0 else params$node_lambda
      if (length(lam) != 1L || !is.finite(lam) || lam < 0) return(-Inf)
      prob <- .cs_table_distribution(table_local, params)
      value <- stats::dpois(b_local, lam, log = TRUE) + prob$log_probs[state_local] +
        .cs_attribute_log(params, observed_local, levels_local)
      if (nonempty_local) {
        zprob <- .cs_table_distribution(zero_local, params)
        normalizer <- .cs_nonempty_log(lam, zprob$log_probs)
        if (!is.finite(normalizer)) return(-Inf)
        value <- value - normalizer
      }
      value
    }), envir = cache)
  density_func <- eval(quote(function(params) exp(log_density_func(params))),
                       envir = list2env(list(log_density_func = log_density_func),
                                        parent = baseenv()))
  value <- log_density_func(params)
  marginal_log_probs <- dist$log_probs
  if (condition_nonempty && b == 0L) {
    marginal_log_probs[1L] <- -Inf
    norm <- .cs_logsumexp(marginal_log_probs)
    if (is.finite(norm)) marginal_log_probs <- marginal_log_probs - norm
  }
  out <- list(log_mark_density = value, mark_density = exp(value),
              log_density_func = log_density_func,
              density_func = density_func,
              edge_probs = as.vector(crossprod(exp(marginal_log_probs), tab$states)),
              reference_edge_probs = dist$reference_probs,
              mark_change_stats = if (valid) tab$changes[selected_state, ] + tab$birth_change else NULL,
              n_mark_states = nrow(tab$states), mark_model = "whole_update_cs",
              combined_inputs = NULL)
  if (generate_mark) {
    out$mark_sample <- mark
    out$mark_sample_density <- exp(value)
    out$log_mark_sample_density <- value
  }
  out
}
