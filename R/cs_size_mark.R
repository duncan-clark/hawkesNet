# CS-2: draw a distinct-edge count, then an unordered whole-update mark.
# The Poisson count is conditioned on K <= D, not capped at D. Given B and K,
# structural parameters affect the selected set but cannot affect its count.

.cs_size_log_count <- function(m, d) {
  if (length(m) != 1L || !is.finite(m) || m < 0) {
    stop("m must be a finite nonnegative scalar")
  }
  if (length(d) != 1L || !is.finite(d) || d < 0 || d != as.integer(d)) {
    stop("The number of candidate edges must be a nonnegative integer")
  }
  if (m == 0 || d == 0L) return(c(0, rep(-Inf, d)))
  k <- seq.int(0L, d)
  # Removing the common exp(-m) factor avoids underflow for large m.
  log_weight <- k * log(m) - lgamma(k + 1)
  log_weight - .cs_logsumexp(log_weight)
}

.cs_size_nonempty_log <- function(lambda, m, d_zero) {
  zero_count <- .cs_size_log_count(m, d_zero)
  .cs_logsumexp(c(log(-expm1(-lambda)),
                 -lambda + .cs_logsumexp(zero_count[-1L])))
}

.cs_size_context <- function(old, proposed, time, params, formula_RHS,
                             truncation, mark_decay, growth_only, max_candidates) {
  if (length(max_candidates) != 1L || !is.finite(max_candidates) ||
      max_candidates < 0L || max_candidates > 20L ||
      max_candidates != as.integer(max_candidates)) {
    stop("max_candidates must be an integer between 0 and 20")
  }
  old_n <- network::network.size(old)
  n <- network::network.size(proposed)
  baseline <- network::network.copy(old)
  if (n > old_n) {
    baseline <- network::add.vertices(baseline, n - old_n)
    network::set.vertex.attribute(baseline, "time", c(
      network::get.vertex.attribute(old, "time"), rep(time, n - old_n)))
    for (attr in names(params$vertex_categorical)) {
      network::set.vertex.attribute(baseline, attr,
                                    network::get.vertex.attribute(proposed, attr))
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
  d <- length(tails)
  if (d > max_candidates) {
    stop(sprintf(paste0("Exact size-conditional CS enumerates choose(D,K) states; ",
                        "D=%d exceeds max_candidates=%d. Reduce the candidate ",
                        "window. No independent-edge approximation is substituted."),
                 d, max_candidates))
  }
  times <- if (identical(mark_decay, "activity")) get_latest_times(baseline) else
    network::get.vertex.attribute(baseline, "time")
  ages <- if (identical(mark_decay, "activity")) {
    time - pmax(times[tails], times[heads])
  } else time - times[heads]
  if (any(!is.finite(ages)) || any(ages < 0)) {
    stop("Candidate endpoint ages must be finite and nonnegative")
  }
  list(baseline = baseline, tails = tails, heads = heads, ages = ages, d = d)
}

.cs_size_state_table <- function(context, k, formula_RHS) {
  d <- context$d
  if (length(k) != 1L || !is.finite(k) || k < 0 || k > d ||
      k != as.integer(k)) stop("The edge count must be an integer between 0 and D")
  choices <- if (k == 0L) matrix(integer(0), nrow = 0L, ncol = 1L) else
    utils::combn(d, k)
  states <- matrix(FALSE, ncol(choices), d)
  if (k) {
    states[cbind(rep(seq_len(nrow(states)), each = k), as.vector(choices))] <- TRUE
  }
  model <- .cs_change_model(context$baseline, formula_RHS)
  binary <- model$getNetwork()
  before <- model$statistics()
  changes <- matrix(0, nrow(states), length(before),
                    dimnames = list(NULL, names(before)))
  previous <- rep(FALSE, d)
  for (i in seq_len(nrow(states))) {
    # Telescoping ERNM toggles are only a device to compute each endpoint's
    # exact statistics. No toggle order appears in the mark distribution.
    for (j in which(states[i, ] != previous)) {
      model$dyadUpdate(context$tails[j], context$heads[j])
      binary$setDyads(context$tails[j], context$heads[j], states[i, j])
    }
    changes[i, ] <- model$statistics() - before
    previous <- states[i, ]
  }
  list(states = states, changes = changes, tails = context$tails,
       heads = context$heads, ages = context$ages,
       summed_ages = as.vector(states %*% context$ages),
       baseline_stats = before, k = k, d = d)
}

.cs_size_distribution <- function(tab, params) {
  theta <- as.numeric(params$CS_params)
  if (length(theta) != ncol(tab$changes) || any(!is.finite(theta))) {
    stop("CS_params must match the finite whole-update statistic vector")
  }
  tau <- params$beta_edges
  if (length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("beta_edges must be a finite nonnegative scalar")
  }
  # Edge counts (and any other statistics constant at fixed B,K) cancel.
  # Centering first also makes this cancellation numerically explicit.
  centered <- sweep(tab$changes, 2L, tab$changes[1L, ], "-")
  log_weights <- as.vector(centered %*% theta) -
    tau * (tab$summed_ages - tab$summed_ages[1L])
  peak <- max(log_weights)
  if (!is.finite(peak)) stop("CS whole-update scores must be finite")
  # Normalize after subtracting the peak. Adding log(sum(exp(.))) to an
  # enormous peak first can round that correction away and break total mass.
  centered_scores <- log_weights - peak
  log_sum <- log(sum(exp(centered_scores)))
  log_probs <- centered_scores - log_sum
  list(log_probs = log_probs, probs = exp(log_probs),
       log_normalizer_centered = peak + log_sum)
}

.pmf_cs_size_conditional <- function(time, params, mark_filtration, mark,
                                     generate_mark, formula_RHS, truncation,
                                     mark_decay, growth_only, max_node_time,
                                     max_candidates, condition_nonempty) {
  if (is.null(truncation)) truncation <- 4L
  if (length(truncation) != 1L || !is.finite(truncation) || truncation < 1L ||
      truncation != as.integer(truncation)) {
    stop("truncation must be a positive integer node-window size")
  }
  mark_decay <- match.arg(mark_decay, c("node_entrance", "activity"))
  if (is.null(formula_RHS)) stop("formula_RHS is required for size-conditional CS")
  if (is.null(max_node_time)) max_node_time <- Inf
  empty <- network::network.initialize(0, directed = FALSE)
  network::set.vertex.attribute(empty, "time", numeric(0))
  if (is.null(mark_filtration)) mark_filtration <- empty
  old <- if (generate_mark && !is.null(mark)) network::network.copy(mark) else
    filtration_to_net(mark_filtration, time, equals = FALSE)
  if (is.null(old)) old <- empty
  if (network::is.directed(old)) {
    stop("Size-conditional CS currently supports undirected edge additions")
  }
  old_n <- network::network.size(old)
  zero_context <- .cs_size_context(old, old, time, params, formula_RHS,
                                   truncation, mark_decay, growth_only, max_candidates)
  lambda <- if (time > max_node_time) 0 else params$node_lambda
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("node_lambda must be finite and nonnegative")
  }
  zero_count <- .cs_size_log_count(params$m, zero_context$d)
  log_nonempty <- .cs_size_nonempty_log(lambda, params$m, zero_context$d)

  if (generate_mark) {
    if (condition_nonempty) {
      if (!is.finite(log_nonempty)) stop("No nonempty CS update is possible from this state")
      log_mass_zero_nodes <- -lambda + .cs_logsumexp(zero_count[-1L])
      choose_zero <- stats::runif(1) < exp(log_mass_zero_nodes - log_nonempty)
      b <- if (choose_zero) 0L else as.integer(stats::qpois(
        log(max(stats::runif(1), .Machine$double.xmin)) + log(-expm1(-lambda)),
        lambda, lower.tail = FALSE, log.p = TRUE))
    } else b <- stats::rpois(1, lambda)
    proposed <- network::network.copy(old)
    if (b) {
      proposed <- network::add.vertices(proposed, b)
      network::set.vertex.attribute(proposed, "time", c(
        network::get.vertex.attribute(old, "time"), rep(time, b)))
      proposed <- sample_vertex_attrs(params, old, proposed, old_n, b)
    }
    context <- if (!b) zero_context else .cs_size_context(
      old, proposed, time, params, formula_RHS, truncation, mark_decay,
      growth_only, max_candidates)
    count_log_probs <- .cs_size_log_count(params$m, context$d)
    if (condition_nonempty && b == 0L) count_log_probs[1L] <- -Inf
    k <- sample.int(length(count_log_probs), 1L,
                    prob = exp(count_log_probs - max(count_log_probs))) - 1L
    tab <- .cs_size_state_table(context, k, formula_RHS)
    dist <- .cs_size_distribution(tab, params)
    selected_state <- sample.int(nrow(tab$states), 1L, prob = dist$probs)
    selected <- which(tab$states[selected_state, ])
    if (length(selected)) {
      proposed <- network::add.edges(proposed, tab$tails[selected], tab$heads[selected])
      new_ids <- unlist(lapply(selected, function(j) network::get.edgeIDs(
        proposed, tab$tails[j], tab$heads[j])), use.names = FALSE)
      network::set.edge.attribute(proposed, "time", rep(time, length(new_ids)), e = new_ids)
    }
    mark <- proposed
    valid <- TRUE
  } else {
    if (is.null(mark)) mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
    b <- network::network.size(mark) - old_n
    if (b < 0L || network::is.directed(mark)) {
      stop("CS marks must preserve old vertices and be undirected")
    }
    if (b > 0L) for (attr in names(params$vertex_categorical)) {
      if (!(attr %in% network::list.vertex.attributes(mark)) ||
          anyNA(network::get.vertex.attribute(mark, attr)[seq.int(old_n + 1L, old_n + b)])) {
        stop(sprintf("New vertices must include nonmissing categorical attribute '%s'", attr))
      }
    }
    context <- if (!b) zero_context else .cs_size_context(
      old, mark, time, params, formula_RHS, truncation, mark_decay,
      growth_only, max_candidates)
    old_keys <- .cs_edge_keys(old)
    mark_keys <- .cs_edge_keys(mark)
    added <- setdiff(mark_keys, old_keys)
    candidate_keys <- paste(pmin(context$tails, context$heads),
                            pmax(context$tails, context$heads), sep = "-")
    valid <- all(old_keys %in% mark_keys) && all(added %in% candidate_keys) &&
      !anyDuplicated(mark_keys) && network::network.edgecount(mark) == length(mark_keys)
    selected <- which(candidate_keys %in% added)
    k <- length(selected)
    tab <- .cs_size_state_table(context, k, formula_RHS)
    dist <- .cs_size_distribution(tab, params)
    selected_state <- if (valid) which(rowSums(tab$states[, selected, drop = FALSE]) == k)[1L] else
      NA_integer_
  }
  valid <- valid && !is.na(selected_state) && .cs_preserves_old_vertex_data(old, mark, params) &&
    (!condition_nonempty || b > 0L || k > 0L)
  # Birth-only changes and edge-count changes are constant within this B,K
  # support. Retain them for reporting, although they cancel from its PMF.
  tab$birth_change <- if (b > 0L) {
    before_births <- strip_vertex_attrs_for_ernm(old, formula_RHS, params)
    tab$baseline_stats - .cs_change_model(before_births, formula_RHS)$statistics()
  } else rep(0, length(tab$baseline_stats))
  cat_res <- log_categorical_density(params, mark, old_n, old_n + b)
  observed <- cat_res$observed
  levels <- lapply(names(observed), function(attr) vertex_categorical_level_names(params, attr, mark))
  names(levels) <- names(observed)
  # Density re-evaluation needs scores and counts, not the candidate networks,
  # ERNM objects, dyad lists, or even the enumerated Boolean subset matrix.
  cache <- list2env(list(table_local = tab[c("changes", "summed_ages", "k", "d")],
    d_zero = zero_context$d,
    state_local = selected_state, b_local = b, valid_local = valid,
    no_births = time > max_node_time, observed_local = observed,
    levels_local = levels, nonempty_local = condition_nonempty,
    .cs_size_distribution = .cs_size_distribution,
    .cs_size_log_count = .cs_size_log_count,
    .cs_size_nonempty_log = .cs_size_nonempty_log,
    .cs_attribute_log = .cs_attribute_log), parent = baseenv())
  log_density_func <- eval(quote(function(params) {
    if (!valid_local) return(-Inf)
    lam <- if (no_births) 0 else params$node_lambda
    if (length(lam) != 1L || !is.finite(lam) || lam < 0) return(-Inf)
    prob <- .cs_size_distribution(table_local, params)
    log_count <- .cs_size_log_count(params$m, table_local$d)
    value <- stats::dpois(b_local, lam, log = TRUE) +
      log_count[table_local$k + 1L] + prob$log_probs[state_local] +
      .cs_attribute_log(params, observed_local, levels_local)
    if (nonempty_local) {
      normalizer <- .cs_size_nonempty_log(lam, params$m, d_zero)
      if (!is.finite(normalizer)) return(-Inf)
      value <- value - normalizer
    }
    value
  }), envir = cache)
  density_func <- eval(quote(function(params) exp(log_density_func(params))),
                       envir = list2env(list(log_density_func = log_density_func),
                                        parent = baseenv()))
  value <- log_density_func(params)
  count_log_probs <- .cs_size_log_count(params$m, context$d)
  if (condition_nonempty && b == 0L) {
    count_log_probs[1L] <- -Inf
    norm <- .cs_logsumexp(count_log_probs)
    if (is.finite(norm)) count_log_probs <- count_log_probs - norm
  }
  out <- list(log_mark_density = value, mark_density = exp(value),
              log_density_func = log_density_func, density_func = density_func,
              edge_probs = as.vector(crossprod(dist$probs, tab$states)),
              edge_probs_conditioning = "given_node_and_edge_counts",
              reference_edge_probs = NULL,
              edge_count = k, edge_count_probs = exp(count_log_probs),
              mark_change_stats = if (valid) tab$changes[selected_state, ] + tab$birth_change else NULL,
              n_mark_states = nrow(tab$states), candidate_count = context$d,
              mark_model = "size_conditional_cs", combined_inputs = NULL)
  if (generate_mark) {
    out$mark_sample <- mark
    out$mark_sample_density <- exp(value)
    out$log_mark_sample_density <- value
  }
  out
}
