# CS-1: structural preferences evaluated one edge at a time on the same
# pre-update graph, followed by Poisson attempts and duplicate collapsing.
# The helpers below never enumerate unordered subsets or edge permutations.

.cs_independent_context <- function(old, proposed, time, params, formula_RHS,
                                    truncation, mark_decay, growth_only) {
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
  model <- .cs_change_model(baseline, formula_RHS)
  binary <- model$getNetwork()
  before <- model$statistics()
  changes <- matrix(0, length(tails), length(before),
                    dimnames = list(NULL, names(before)))
  for (j in seq_along(tails)) {
    model$dyadUpdate(tails[j], heads[j])
    binary$setDyads(tails[j], heads[j], TRUE)
    changes[j, ] <- model$statistics() - before
    # Restore the baseline before calculating the next edge's statistics.
    model$dyadUpdate(tails[j], heads[j])
    binary$setDyads(tails[j], heads[j], FALSE)
  }
  birth_change <- if (n > old_n) {
    before_births <- strip_vertex_attrs_for_ernm(
      network::network.copy(old), formula_RHS, params)
    before - .cs_change_model(before_births, formula_RHS)$statistics()
  } else rep(0, length(before))
  times <- if (identical(mark_decay, "activity")) get_latest_times(baseline) else
    network::get.vertex.attribute(baseline, "time")
  ages <- if (identical(mark_decay, "activity")) {
    time - pmax(times[tails], times[heads])
  } else time - times[heads]
  if (any(!is.finite(ages)) || any(ages < 0)) {
    stop("Candidate endpoint ages must be finite and nonnegative")
  }
  list(tails = tails, heads = heads, changes = changes, ages = ages,
       baseline_stats = before, birth_change = birth_change)
}

.cs_independent_distribution <- function(tab, params) {
  theta <- as.numeric(params$CS_params)
  if (length(theta) != ncol(tab$changes) || any(!is.finite(theta))) {
    stop("CS_params must match the finite single-edge statistic vector")
  }
  m <- params$m
  tau <- params$beta_edges
  if (length(m) != 1L || !is.finite(m) || m < 0 ||
      length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("m and beta_edges must be finite nonnegative scalars")
  }
  if (!length(tab$tails)) {
    return(list(probs = numeric(0), log_probs = numeric(0), rates = numeric(0),
                selection_probs = numeric(0), log_nonempty = -Inf))
  }
  score <- as.vector(tab$changes %*% theta)
  if (any(!is.finite(score))) stop("Single-edge structural scores must be finite")
  # Subtracting the common minimum age leaves normalized weights unchanged.
  # log.p avoids underflow even when all logistic scores are extremely small.
  log_weights <- -tau * (tab$ages - min(tab$ages)) +
    stats::plogis(score, log.p = TRUE)
  centered <- log_weights - max(log_weights)
  log_p <- centered - log(sum(exp(centered)))
  log_rates <- log(m) + log_p
  rates <- exp(log_rates)
  log_probs <- log(-expm1(-rates))
  # exp(log_rates) can underflow while the selected edge's log PMF remains
  # representable. In this regime log(1-exp(-x)) equals log(x) numerically.
  tiny <- log_rates < log(.Machine$double.xmin)
  log_probs[tiny] <- log_rates[tiny]
  list(probs = -expm1(-rates), log_probs = log_probs, rates = rates,
       selection_probs = exp(log_p), log_nonempty = log(-expm1(-m)))
}

.cs_independent_log_prob <- function(dist, selected) {
  take <- seq_along(dist$rates) %in% selected
  sum(dist$log_probs[take]) - sum(dist$rates[!take])
}

.cs_independent_nonempty_log <- function(lambda, has_candidates, m) {
  log(-expm1(-(lambda + if (has_candidates) m else 0)))
}

.cs_independent_draw <- function(dist, nonempty = FALSE) {
  d <- length(dist$rates)
  if (!nonempty) return(which(stats::runif(d) < dist$probs))
  if (!d || !is.finite(dist$log_nonempty)) {
    stop("No nonempty CS edge set is possible from this state")
  }
  # Sample the first included candidate given at least one inclusion. Earlier
  # candidates are absent and later ones retain their independent Bernoulli
  # laws. This is a sampling device, not an ordering of edges in the mark PMF.
  previous_rate <- c(0, head(cumsum(dist$rates), -1L))
  log_first <- dist$log_probs - previous_rate
  first <- sample.int(d, 1L, prob = exp(log_first - max(log_first)))
  remaining <- if (first < d) seq.int(first + 1L, d) else integer(0)
  c(first, remaining[stats::runif(length(remaining)) < dist$probs[remaining]])
}

.pmf_cs_independent <- function(time, params, mark_filtration, mark,
                                generate_mark, formula_RHS, truncation,
                                mark_decay, growth_only, max_node_time,
                                condition_nonempty) {
  if (is.null(truncation)) truncation <- 4L
  if (length(truncation) != 1L || !is.finite(truncation) || truncation < 1L ||
      truncation != as.integer(truncation)) {
    stop("truncation must be a positive integer node-window size")
  }
  mark_decay <- match.arg(mark_decay, c("node_entrance", "activity"))
  if (is.null(formula_RHS)) stop("formula_RHS is required for independent CS")
  if (is.null(max_node_time)) max_node_time <- Inf
  empty <- network::network.initialize(0, directed = FALSE)
  network::set.vertex.attribute(empty, "time", numeric(0))
  if (is.null(mark_filtration)) mark_filtration <- empty
  old <- if (generate_mark && !is.null(mark)) network::network.copy(mark) else
    filtration_to_net(mark_filtration, time, equals = FALSE)
  if (is.null(old)) old <- empty
  if (network::is.directed(old)) stop("Independent CS currently supports undirected edge additions")
  old_n <- network::network.size(old)
  zero_tab <- .cs_independent_context(old, old, time, params, formula_RHS,
                                     truncation, mark_decay, growth_only)
  zero_dist <- .cs_independent_distribution(zero_tab, params)
  lambda <- if (time > max_node_time) 0 else params$node_lambda
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("node_lambda must be finite and nonnegative")
  }
  has_zero_candidates <- length(zero_tab$tails) > 0L
  log_nonempty <- .cs_independent_nonempty_log(lambda, has_zero_candidates, params$m)

  if (generate_mark) {
    if (condition_nonempty) {
      if (!is.finite(log_nonempty)) stop("No nonempty CS update is possible from this state")
      choose_zero <- stats::runif(1) <
        exp(-lambda + zero_dist$log_nonempty - log_nonempty)
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
    tab <- if (!b) zero_tab else .cs_independent_context(
      old, proposed, time, params, formula_RHS, truncation, mark_decay, growth_only)
    dist <- if (!b) zero_dist else .cs_independent_distribution(tab, params)
    selected <- .cs_independent_draw(dist, condition_nonempty && b == 0L)
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
    tab <- if (!b) zero_tab else .cs_independent_context(
      old, mark, time, params, formula_RHS, truncation, mark_decay, growth_only)
    dist <- if (!b) zero_dist else .cs_independent_distribution(tab, params)
    old_keys <- .cs_edge_keys(old)
    mark_keys <- .cs_edge_keys(mark)
    added <- setdiff(mark_keys, old_keys)
    candidate_keys <- paste(pmin(tab$tails, tab$heads), pmax(tab$tails, tab$heads), sep = "-")
    valid <- all(old_keys %in% mark_keys) && all(added %in% candidate_keys) &&
      !anyDuplicated(mark_keys) && network::network.edgecount(mark) == length(mark_keys)
    selected <- which(candidate_keys %in% added)
  }
  valid <- valid && .cs_preserves_old_vertex_data(old, mark, params) &&
    (!condition_nonempty || b > 0L || length(selected) > 0L)
  cat_res <- log_categorical_density(params, mark, old_n, old_n + b)
  observed <- cat_res$observed
  levels <- lapply(names(observed), function(attr) vertex_categorical_level_names(params, attr, mark))
  names(levels) <- names(observed)
  # Retain numeric features only. No network, ERNM pointer, state table, or
  # sampled auxiliary ordering is needed for subsequent likelihood calls.
  cache <- list2env(list(table_local = tab, selected_local = selected,
    b_local = b, valid_local = valid, no_births = time > max_node_time,
    has_zero_candidates = has_zero_candidates, observed_local = observed,
    levels_local = levels, nonempty_local = condition_nonempty,
    .cs_independent_distribution = .cs_independent_distribution,
    .cs_independent_log_prob = .cs_independent_log_prob,
    .cs_independent_nonempty_log = .cs_independent_nonempty_log,
    .cs_attribute_log = .cs_attribute_log), parent = baseenv())
  log_density_func <- eval(quote(function(params) {
    if (!valid_local) return(-Inf)
    lam <- if (no_births) 0 else params$node_lambda
    if (length(lam) != 1L || !is.finite(lam) || lam < 0) return(-Inf)
    prob <- .cs_independent_distribution(table_local, params)
    value <- stats::dpois(b_local, lam, log = TRUE) +
      .cs_independent_log_prob(prob, selected_local) +
      .cs_attribute_log(params, observed_local, levels_local)
    if (nonempty_local) {
      normalizer <- .cs_independent_nonempty_log(lam, has_zero_candidates, params$m)
      if (!is.finite(normalizer)) return(-Inf)
      value <- value - normalizer
    }
    value
  }), envir = cache)
  density_func <- eval(quote(function(params) exp(log_density_func(params))),
                       envir = list2env(list(log_density_func = log_density_func),
                                        parent = baseenv()))
  value <- log_density_func(params)
  marginal_probs <- dist$probs
  if (condition_nonempty && b == 0L) {
    marginal_probs <- if (is.finite(dist$log_nonempty))
      exp(dist$log_probs - dist$log_nonempty) else rep(0, length(dist$probs))
  }
  whole_change <- if (valid) {
    final <- strip_vertex_attrs_for_ernm(network::network.copy(mark), formula_RHS, params)
    .cs_change_model(final, formula_RHS)$statistics() - tab$baseline_stats + tab$birth_change
  } else NULL
  d <- length(tab$tails)
  out <- list(log_mark_density = value, mark_density = exp(value),
              log_density_func = log_density_func, density_func = density_func,
              edge_probs = marginal_probs, reference_edge_probs = dist$probs,
              edge_selection_probs = dist$selection_probs,
              edge_change_stats = tab$changes, mark_change_stats = whole_change,
              n_candidate_edges = d, n_mark_states = 2^d,
              mark_model = "independent_cs", combined_inputs = NULL)
  if (generate_mark) {
    out$mark_sample <- mark
    out$mark_sample_density <- exp(value)
    out$log_mark_sample_density <- value
  }
  out
}
