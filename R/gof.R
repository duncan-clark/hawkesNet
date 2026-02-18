# =============================================================================
# Goodness-of-Fit (GOF) functions for hawkesNet
# =============================================================================

#' Degree distribution as vector of counts (degree min_deg, ..., max_deg).
#'
#' @param net Network object.
#' @param max_deg Maximum degree to count (default 20).
#' @param min_deg Minimum degree to count (default 0).
#' @return Numeric vector of counts for degrees min_deg through max_deg.
#' @noRd
degree_dist <- function(net, max_deg = 20, min_deg = 0) {
  degs <- degree(net, gmode = "graph")
  tab <- table(factor(degs, levels = min_deg:max_deg))
  as.vector(tab)
}

#' ESP distribution via ernm (edge-wise shared partners min_esp, ..., k_max).
#'
#' @param net Network object.
#' @param k_max Maximum ESP count (default 15).
#' @param min_esp Minimum ESP to count (default 0).
#' @return Numeric vector of ESP counts for min_esp through k_max.
#' @noRd
esp_dist <- function(net, k_max = 15, min_esp = 0) {
  tryCatch({
    as.vector(calculateStatistics(net ~ esp(min_esp:k_max)))
  }, error = function(e) rep(NA_real_, k_max - min_esp + 1))
}

#' Geodesic distance distribution (upper triangle of distance matrix, excluding Inf).
#'
#' @param net Network object.
#' @return Numeric vector of geodesic distances (excluding Inf and self-distances).
#' @noRd
geodist_dist <- function(net) {
  d <- geodist(net, inf.replace = NA)
  if (is.list(d)) d <- d$gdist
  d <- as.vector(d)
  d <- d[!is.na(d) & d > 0]
  if (length(d) == 0) return(numeric(0))
  d
}

#' Ensure vertex attribute is properly set for all nodes before ERNM operations.
#'
#' Ensures all nodes have the attribute, replaces NA/missing values with default,
#' and removes problematic 'na' attribute.
#'
#' @param net Network object.
#' @param attr_name Name of vertex attribute to ensure.
#' @param default_value Default value to use for missing/NA values (default "unknown").
#' @return Network object with attribute properly set.
#' @noRd
ensure_vertex_attribute <- function(net, attr_name, default_value = "unknown") {
  if (is.null(net)) return(net)
  nv <- network.size(net)
  if (nv == 0) return(net)
  
  # Remove problematic 'na' attribute if it exists
  if ("na" %in% list.vertex.attributes(net)) {
    delete.vertex.attribute(net, "na")
  }
  
  # Check if attribute exists
  if (!attr_name %in% list.vertex.attributes(net)) {
    set.vertex.attribute(net, attr_name, rep(default_value, nv))
  } else {
    attr_vals <- get.vertex.attribute(net, attr_name)
    if (length(attr_vals) < nv) {
      attr_vals <- c(attr_vals, rep(default_value, nv - length(attr_vals)))
      set.vertex.attribute(net, attr_name, attr_vals)
    } else if (any(is.na(attr_vals)) || any(attr_vals == "")) {
      attr_vals[is.na(attr_vals) | attr_vals == ""] <- default_value
      set.vertex.attribute(net, attr_name, attr_vals)
    }
    if (!is.numeric(attr_vals)) {
      attr_vals <- sanitize_vertex_attr_for_binarynet(attr_vals, default_val = default_value)
      set.vertex.attribute(net, attr_name, attr_vals)
    }
  }
  
  net
}

#' Vertex attributes required by an ERNM formula (nodeMatch/nodeMix).
#' @param formula_RHS Character RHS of formula.
#' @return Character vector of attribute names, or character(0) if none.
#' @noRd
formula_vertex_attrs <- function(formula_RHS) {
  if (is.null(formula_RHS) || !nzchar(trimws(formula_RHS))) return(character(0))
  m <- gregexpr("node(?:Match|Mix)\\s*\\(\\s*['\"]([^'\"]+)['\"]", formula_RHS, perl = TRUE)[[1]]
  if (m[1] == -1) return(character(0))
  s <- attr(m, "capture.start")
  l <- attr(m, "capture.length")
  unique(substring(formula_RHS, s[, 1], s[, 1] + l[, 1] - 1))
}

#' Sanitize character vertex attributes for ernm as.BinaryNet (avoids C++ segfaults).
#' Strips non-ASCII and truncates long strings; replaces NA/empty with default.
#' @param vals Character or numeric vector.
#' @param default_val Default for NA/empty (default "unknown").
#' @param max_len Max character length per element (default 200).
#' @return Sanitized vector.
#' @noRd
sanitize_vertex_attr_for_binarynet <- function(vals, default_val = "unknown", max_len = 200L) {
  if (is.numeric(vals)) {
    vals[!is.finite(vals)] <- 0
    return(vals)
  }
  vals <- as.character(vals)
  vals[is.na(vals) | nchar(vals) == 0] <- default_val
  vals <- vapply(vals, function(x) {
    x <- iconv(x, from = "UTF-8", to = "ASCII", sub = "?")
    if (is.na(x)) return(default_val)
    if (nchar(x) > max_len) x <- substring(x, 1L, max_len)
    x
  }, character(1), USE.NAMES = FALSE)
  vals
}

#' Sanitize all vertex attributes on a network before passing to ernm as.BinaryNet.
#' Prevents C++ segfaults from non-ASCII, long strings, or NA in vertex attributes.
#' Modifies net in place and returns it.
#' @param net A network object.
#' @return The same network with sanitized vertex attributes.
#' @noRd
sanitize_net_for_binarynet <- function(net) {
  if (is.null(net) || network::network.size(net) == 0L) return(net)
  nv <- network::network.size(net)
  all_attrs <- network::list.vertex.attributes(net)
  all_attrs <- setdiff(all_attrs, "na")
  for (a in all_attrs) {
    vals <- network::get.vertex.attribute(net, a)
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
    network::set.vertex.attribute(net, a, vals)
  }
  if ("na" %in% network::list.vertex.attributes(net)) {
    network::delete.vertex.attribute(net, "na")
  }
  net
}

#' Waiting times between consecutive structure formations.
#'
#' Replays the network event-by-event (grouped by event time) and records
#' when statistics (triangles, stars, etc.) increase, then returns differences
#' between formation times.
#'
#' @param net Network with vertex/edge time attributes.
#' @param time_attr Name of the time attribute (default "time").
#' @param formula_RHS Character string RHS of ERNM formula (e.g., "triangles + star(c(2,3))").
#' @return List with named vectors of waiting times for each statistic.
#' @noRd
waiting_times_between_formations <- function(net, time_attr = "time", 
                                             formula_RHS = "triangles + star(c(2,3))") {
  # Get edge list and edge times
  el <- as.edgelist(net)
  
  # Determine number of statistics by creating a temporary model
  n <- network.size(net)
  g_temp <- network.initialize(n, directed = is.directed(net))
  
  # Copy ONLY vertex attributes required by the formula (nodeMatch/nodeMix).
  # For "edges + triangles + star(c(2,3))" none are needed; copying author/title
  # etc. causes ernm C++ segfaults with non-ASCII or long strings.
  needed_attrs <- formula_vertex_attrs(formula_RHS)
  vattrs <- list.vertex.attributes(net)
  for (attr_name in vattrs) {
    if (attr_name != "na" && attr_name %in% needed_attrs) {
      attr_vals <- get.vertex.attribute(net, attr_name)
      if (length(attr_vals) == n) {
        attr_vals <- sanitize_vertex_attr_for_binarynet(attr_vals)
        set.vertex.attribute(g_temp, attr_name, attr_vals)
      }
    }
  }
  if ("na" %in% list.vertex.attributes(g_temp)) {
    delete.vertex.attribute(g_temp, "na")
  }
  
  # Create formula with g_temp in the environment
  formula_str <- paste0("g_temp ~ ", formula_RHS)
  formula_obj <- as.formula(formula_str)
  model_temp <- createCppModel(formula_obj)
  model_temp$setNetwork(as.BinaryNet(g_temp))
  model_temp$calculate()
  n_stats <- length(model_temp$statistics())
  rm(model_temp, g_temp)
  
  if (nrow(el) == 0) {
    # Return empty structure matching formula
    return(setNames(rep(list(numeric()), n_stats), paste0("stat", seq_len(n_stats))))
  }

  # Try edge times first; fall back to vertex times for edge ordering
  edge_times <- NULL
  if (time_attr %in% list.edge.attributes(net)) {
    edge_times <- get.edge.attribute(net, time_attr)
  }
  if (is.null(edge_times) || length(edge_times) != nrow(el)) {
    # Use vertex times: assign each edge the max time of its endpoints
    vtimes <- get.vertex.attribute(net, time_attr)
    if (is.null(vtimes) || all(is.na(vtimes))) {
      return(setNames(rep(list(numeric()), n_stats), paste0("stat", seq_len(n_stats))))
    }
    edge_times <- pmax(vtimes[el[, 1]], vtimes[el[, 2]], na.rm = TRUE)
  }

  # Sort edges by time
  ord <- order(edge_times)
  el <- el[ord, , drop = FALSE]
  edge_times <- edge_times[ord]

  # Group edges by unique event times
  event_times <- unique(edge_times)
  g <- network.initialize(n, directed = is.directed(net))
  
  # Copy ONLY vertex attributes required by the formula (same as g_temp)
  for (attr_name in needed_attrs) {
    if (attr_name %in% list.vertex.attributes(net)) {
      attr_vals <- get.vertex.attribute(net, attr_name)
      if (length(attr_vals) == n) {
        attr_vals <- sanitize_vertex_attr_for_binarynet(attr_vals)
        set.vertex.attribute(g, attr_name, attr_vals)
      }
    }
  }
  if ("na" %in% list.vertex.attributes(g)) {
    delete.vertex.attribute(g, "na")
  }

  # Initialize tracking for all statistics
  stat_prev <- rep(0, n_stats)
  stat_times <- rep(list(numeric()), n_stats)

  # Create ERNM model once and initialize with empty network
  # Create formula with g now that it exists
  formula_str <- paste0("g ~ ", formula_RHS)
  formula_obj <- as.formula(formula_str)
  model <- createCppModel(formula_obj)
  model$setNetwork(as.BinaryNet(g))
  model$calculate()
  
  # Pre-allocate stat_times if we have many events (optional, but good for large nets)
  # For now, we'll keep the list of vectors.

  for (t_cur in event_times) {
    idx <- which(edge_times == t_cur)
    new_edges <- el[idx, , drop = FALSE]
    
    if (nrow(new_edges) > 0) {
      # Incremental approach: add edges one at a time, get change stats, update model
      # This correctly handles interactions (e.g., triangles formed by multiple edges)
      stat_curr <- stat_prev
      for (j in seq_len(nrow(new_edges))) {
        tail_j <- new_edges[j, 1]
        head_j <- new_edges[j, 2]
        
        # Get change stats BEFORE adding this edge (from current network state)
        change_stats_j <- model$computeChangeStats(tail_j, head_j)
        # change_stats_j is a matrix: rows = 1 (single edge), columns = statistics
        delta_stats <- change_stats_j[1, ]
        
        # Add edge to network
        add.edges(g, tail = tail_j, head = head_j)
        
        # Update model to reflect this new edge (for next iteration's change stats)
        # OPTIMIZATION: Use model$toggle() if available in ERNM, otherwise calculate()
        # ERNM's computeChangeStats + manual add.edges + calculate is standard.
        model$setNetwork(as.BinaryNet(g))
        model$calculate()
        
        # Update statistics incrementally
        stat_curr <- stat_curr + as.integer(round(delta_stats))
      }
    } else {
      stat_curr <- stat_prev
    }
    
    # Track when each statistic increases
    for (i in seq_len(n_stats)) {
      if (stat_curr[i] > stat_prev[i]) {
        # Optimization: append to list of times
        stat_times[[i]] <- c(stat_times[[i]], t_cur)
      }
    }
    
    stat_prev <- stat_curr
  }

  # Compute waiting times (differences between formation times)
  out <- lapply(stat_times, function(t_vec) {
    if (length(t_vec) >= 2) diff(t_vec) else numeric()
  })
  
  # Name outputs based on formula
  # For backward compatibility, use triangle, star2, star3 if formula matches
  if (formula_RHS == "triangles + star(c(2,3))") {
    names(out) <- c("triangle", "star2", "star3")
  } else if (formula_RHS == "edges + triangles + star(c(2,3))") {
    names(out) <- c("edges", "triangle", "star2", "star3")
  } else {
    # Generic naming: try to extract statistic names from formula
    # This is a simple heuristic - may not work for all formulas
    stat_names <- trimws(strsplit(formula_RHS, "\\+")[[1]])
    if (length(stat_names) == length(out)) {
      names(out) <- stat_names
    } else {
      names(out) <- paste0("stat", seq_len(length(out)))
    }
  }
  
  out
}

#' Goodness-of-fit analysis for hawkesNet models.
#'
#' Simulates networks from a fitted model and compares observed vs simulated
#' network statistics (degree distribution, ESP, geodesic distances, waiting times).
#'
#' @param fit Fitted hawkesNet model object (from \code{fit_hawkesNet} or \code{fit_hawkesNet_inhom}).
#' @param net_obs Observed network (used for computing observed statistics).
#' @param params_init Initial parameter list used for fitting (needed to reconstruct full parameter structure).
#' @param PMF_mark Mark probability mass function (e.g., \code{PMF_mark_CS}).
#' @param cond_intensity Conditional intensity function (cached closure from fitting).
#'   If \code{inhom_bg} is provided, \code{cond_intensity_inhom} will be used automatically for simulations.
#' @param formula_RHS Character string RHS of ERNM formula (e.g., "edges + triangles + star(c(2,3))").
#' @param time_window Time window for simulations (default c(0, 0.05)).
#' @param truncation Truncation parameter for mark PMF (default 100).
#' @param mark_decay Mark decay type: "node_entrance" or "activity" (default "activity").
#' @param max_node_time Maximum node time (default 1).
#' @param inhom_bg Optional inhomogeneous background object from \code{prepare_inhomogeneous_background}.
#'   If provided, simulations use \code{cond_intensity_inhom} with time-varying background rate to match the fitted model.
#' @param n_sim Number of simulated networks to generate (default 50).
#' @param cores Number of cores for parallelization (default 7). Ignored when \code{cores_outer} is set.
#' @param cores_outer Optional. When set, use PSOCK cluster with this many workers instead of fork.
#'   Avoids BLAS/fork deadlocks on high-core machines (e.g. 100 cores). Capped at 60 for socket limits.
#' @param max_deg Maximum degree for degree distribution (default 15).
#' @param k_esp Maximum ESP count for ESP distribution (default 15).
#' @param degree Minimum degree to include in degree distribution (default 0).
#' @param esp Minimum ESP to include in ESP distribution (default 0).
#' @param mu_multiplier Multiplier for mu in simulations (default 5).
#' @param growth_only Logical; if \code{TRUE}, edges only form when a node enters the network.
#'   Default \code{FALSE}.
#' @param seed_events Optional integer. If \code{> 0}, the first \code{seed_events}
#'   observed events are used to seed the simulation. The simulation then
#'   generates subsequent events conditional on this initial history.
#'   Recommended for sparse networks where cold-start simulation is difficult.
#' @param verbose Print progress messages (default TRUE).
#' @return List with observed and simulated statistics and plots.
#' @examples
#' \donttest{
#' params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
#' set.seed(1)
#' sim <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
#'                      verbose = FALSE, mu_multiplier = 5, truncation = 30)
#' fit <- fit_hawkesNet(params, c(0, 3), sim$net, PMF_mark_BA,
#'                      maxit = 10, verbose = FALSE, truncation = 30)
#' # Run GOF with 2 simulations for speed
#' g <- gof(fit, sim$net, params, PMF_mark_BA, cond_intensity,
#'          time_window = c(0, 3), n_sim = 2, truncation = 30)
#' names(g$plots)
#' }
#' @export
gof <- function(fit, net_obs, params_init, PMF_mark, cond_intensity, formula_RHS = NULL,
                time_window = c(0, 0.05), truncation = 100L, mark_decay = "activity",
                growth_only = FALSE,
                max_node_time = 1, inhom_bg = NULL, n_sim = 50L, cores = 7L,
                cores_outer = NULL,
                max_deg = 15L, k_esp = 15L, degree = 0L, esp = 0L, mu_multiplier = 5,
                seed_events = 0L, verbose = TRUE) {
  
  # -------------------------------------------------------------------------
  # Force evaluation of user-provided arguments that may be promises.
  #
  # In scripts it's common to call e.g. growth_only = GROWTH_ONLY, mark_decay = MARK_DECAY.
  # Under PSOCK parallelism, those promises can get serialized and evaluated later inside
  # workers that do not have the originating symbols, producing errors like:
  #   "object 'GROWTH_ONLY' not found"
  # Coerce here so the values are concrete scalars captured in the closure.
  # -------------------------------------------------------------------------
  growth_only <- isTRUE(growth_only)
  mark_decay  <- as.character(mark_decay)
  truncation  <- as.integer(truncation)
  max_node_time <- as.numeric(max_node_time)
  n_sim <- as.integer(n_sim)
  cores <- as.integer(cores)
  cores_outer <- if (is.null(cores_outer)) NULL else as.integer(cores_outer)
  use_psock <- !is.null(cores_outer) && cores_outer > 0L
  n_workers <- if (use_psock) min(as.integer(cores_outer), 60L) else cores
  max_deg <- as.integer(max_deg)
  k_esp <- as.integer(k_esp)
  degree <- as.integer(degree)
  esp <- as.integer(esp)
  mu_multiplier <- as.numeric(mu_multiplier)
  seed_events <- as.integer(seed_events)
  if (!is.null(formula_RHS)) formula_RHS <- as.character(formula_RHS)
  
  # Initialize results early (will be populated even if some computations fail)
  GOF_results <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                      geodist_obs = NULL, geodist_sim = NULL, wait_obs = NULL, wait_sim = NULL,
                      nodemix_obs = NULL, nodemix_sim = NULL, nets_sim = NULL, plots = list())
  
  if (is.null(fit)) {
    if (verbose) cat("  No fit available; skipping GOF\n")
    return(GOF_results)
  }
  
  if (verbose) {
    cond_str <- if (seed_events > 0) sprintf(" (conditional on first %d events)", seed_events) else ""
    cat("  Simulating", n_sim, "networks from fitted model", cond_str, "...\n", sep = "")
  }
  
  # -------------------------------------------------------------------------
  # Reconstruct fitted params from fit$par (name-based, robust).
  # fit_hawkesNet strips fixed_params + vertex_categorical_levels from the
  # skeleton before optim, so fit$par only has the FREE parameters.
  # Instead of fragile relist() (which breaks if skeleton doesn't match
  # exactly), we map values back by name.
  # -------------------------------------------------------------------------
  
  fixed <- fit$fixed_params
  if (verbose) {
    if (!is.null(fixed)) {
      cat("  Fixed parameters:", paste(fixed, collapse = ", "), "\n")
    } else {
      cat("  WARNING: fit$fixed_params is NULL -- was this fit run with the latest code?\n")
    }
  }
  
  par_vec <- fit$fit$par
  par_names <- names(par_vec)
  
  if (verbose) {
    cat("  fit$par (", length(par_vec), " values):", paste(par_names, "=",
        round(par_vec, 6), collapse = ", "), "\n")
  }
  
  # Start from params_init and overwrite with fitted values by name
  pfit <- params_init
  pfit$vertex_categorical_levels <- params_init$vertex_categorical_levels
  
  # Map scalar parameters directly
  scalar_names <- c("mu", "beta_overall", "K", "beta_edges", "node_lambda", "m")
  for (nm in scalar_names) {
    if (nm %in% par_names) {
      pfit[[nm]] <- par_vec[nm]
    }
    # If not in par_vec, keep params_init value (it's fixed or absent)
  }
  
  # Map CS_params: look for CS_params1, CS_params2, ... in par_vec
  cs_idx <- grep("^CS_params[0-9]+$", par_names)
  if (length(cs_idx) > 0) {
    cs_nums <- as.integer(sub("^CS_params", "", par_names[cs_idx]))
    n_cs <- length(pfit$CS_params)
    for (j in seq_along(cs_idx)) {
      k <- cs_nums[j]
      if (k >= 1L && k <= n_cs) {
        pfit$CS_params[k] <- par_vec[cs_idx[j]]
      }
    }
    if (verbose) {
      cat("  Mapped", length(cs_idx), "CS_params from fit$par into",
          n_cs, "slots\n")
    }
  }
  
  # Map vertex_categorical: look for vertex_categorical.ATTR.LEVEL names
  vc_idx <- grep("^vertex_categorical\\.", par_names)
  if (length(vc_idx) > 0 && !is.null(pfit$vertex_categorical)) {
    for (j in vc_idx) {
      parts <- strsplit(par_names[j], "\\.")[[1]]
      if (length(parts) >= 3) {
        attr_name <- parts[2]
        level_name <- paste(parts[3:length(parts)], collapse = ".")
        if (!is.null(pfit$vertex_categorical[[attr_name]])) {
          pfit$vertex_categorical[[attr_name]][level_name] <- par_vec[j]
        }
      }
    }
    if (verbose) {
      cat("  Mapped", length(vc_idx), "vertex_categorical values from fit$par\n")
    }
  }
  
  # Handle inhomogeneous background
  use_inhom <- !is.null(inhom_bg) && !is.null(inhom_bg$mu_fit) && !is.null(inhom_bg$mu_fit$mu_fun)
  
  if (use_inhom) {
    # For inhomogeneous: mu is not used directly during simulation (mu_fun is),
    # but we set it for thinning bound and compatibility
    Tval <- time_window[2] - time_window[1]
    pfit$mu <- inhom_bg$integral_bg / Tval  # Average mu for compatibility
    if (verbose) {
      cat("  Using inhomogeneous background for simulations (matches fitted model)\n")
    }
  }
  
  # Ensure parameters are within valid ranges
  pfit$K <- min(max(pfit$K, 0.001), 0.999)
  pfit$mu <- max(pfit$mu, 0.001)
  pfit$node_lambda <- max(pfit$node_lambda, 0.1)
  pfit$beta_overall <- max(pfit$beta_overall, 0.001)
  pfit$beta_edges <- max(pfit$beta_edges, 0.001)
  
  # Log the actual parameters being used for GOF simulations
  if (verbose) {
    cat("  GOF simulation parameters:\n")
    scalar_params <- c("mu", "beta_overall", "K", "beta_edges", "node_lambda")
    for (p in scalar_params) {
      if (!is.null(pfit[[p]])) {
        src <- if (!is.null(fixed) && p %in% fixed) "(FIXED)" else "(fitted)"
        cat(sprintf("    %s = %.6f %s\n", p, pfit[[p]], src))
      }
    }
    if (!is.null(pfit$CS_params)) {
      cat("    CS_params =", paste(round(pfit$CS_params, 4), collapse = ", "), "(fitted)\n")
    }
  }
  
  # -------------------------------------------------------------------------
  # Handle Seeding (Conditional Simulation)
  # -------------------------------------------------------------------------
  seed_net <- NULL
  seed_times <- NULL
  if (seed_events > 0) {
    all_times <- get_times(net_obs)$times
    if (length(all_times) >= seed_events) {
      t_seed <- all_times[seed_events]
      seed_net <- filtration_to_net(net_obs, t_seed, equals = TRUE)
      seed_times <- all_times[1:seed_events]
      if (verbose) {
        cat(sprintf("  Seeding simulation with first %d events (up to t=%.4f)\n", 
                    seed_events, t_seed))
      }
    } else {
      if (verbose) cat(sprintf("  WARNING: net_obs only has %d events; cannot seed with %d. Starting from scratch.\n",
                               length(all_times), seed_events))
    }
  }
  
  # Set vertex_categorical if it exists (use defaults if needed)
  if (!is.null(params_init$vertex_categorical)) {
    if (is.null(pfit$vertex_categorical)) {
      pfit$vertex_categorical <- params_init$vertex_categorical
    }
  }
  
  # Validate and repair parameters before simulation
  if (!point_process_params_valid(pfit)) {
    if (verbose) cat("  WARNING: Parameters invalid after reconstruction; attempting repair...\n")
    pfit <- tryCatch({
      repair_vertex_categorical_params(pfit, eps = 1e-6)
    }, error = function(e) {
      if (verbose) cat("  WARNING: Failed to repair parameters:", e$message, "\n")
      pfit
    })
    if (!point_process_params_valid(pfit)) {
      if (verbose) cat("  ERROR: Parameters still invalid after repair; GOF may fail\n")
    }
  }
  
  # Parallelize GOF simulations
  if (verbose) {
    cat("  Using", n_workers, if (use_psock) "PSOCK workers" else "cores", "for parallel GOF simulations...\n")
    if (use_inhom) {
      cat("  Simulations will use inhomogeneous background (cond_intensity_inhom)\n")
    } else {
      cat("  Simulations will use homogeneous background (cond_intensity)\n")
    }
  }
  sim_fun <- function(i) {
    s <- tryCatch({
      sim_hawkesNet(
          params = pfit,
          time_window = time_window,
          PMF_mark = PMF_mark,
          cond_intensity = cond_intensity,
          formula_RHS = formula_RHS,
          truncation = truncation,
          mark_decay = mark_decay,
          growth_only = growth_only,
          max_node_time = max_node_time,
          hashed_edges = TRUE,
          verbose = FALSE,
          mu_multiplier = mu_multiplier,
          stop_on_full_network = FALSE,
          inhom_bg = inhom_bg,
          seed_net = seed_net,
          seed_times = seed_times
      )
    }, error = function(e) {
      return(list(net = NULL, error = paste0("Sim ", i, ": ", e$message)))
    })
    if (is.null(s$net)) {
      return(list(net = NULL, error = ifelse(is.null(s$error), paste0("Sim ", i, ": unknown error"), s$error)))
    }
    list(net = s$net, error = NULL)
  }
  cl_gof <- NULL
  if (use_psock) {
    message(sprintf("  [GOF] Creating PSOCK cluster (%d workers) at %s", n_workers, format(Sys.time(), "%H:%M:%S")))
    # Find package root so workers load dev version (avoids "unused argument" when installed pkg is stale)
    pkg_path <- tryCatch({
      p <- getwd()
      if (file.exists(file.path(p, "DESCRIPTION"))) p
      else if (file.exists(file.path(p, "..", "DESCRIPTION"))) normalizePath(file.path(p, ".."), mustWork = TRUE)
      else NULL
    }, error = function(e) NULL)
    cl_gof <- tryCatch({
      cl <- parallel::makeCluster(n_workers, outfile = "")
      if (!is.null(pkg_path)) {
        env_export <- new.env()
        env_export$pkg_path <- pkg_path
        parallel::clusterExport(cl, "pkg_path", envir = env_export)
      }
      # library() is needed here to attach packages on PSOCK workers (separate R processes)
      parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages({
          if (!is.null(pkg_path) && nzchar(pkg_path) && requireNamespace("devtools", quietly = TRUE)) {
            devtools::load_all(pkg_path, quiet = TRUE)
          } else {
            library(hawkesNet)
          }
          library(network)
          library(ernm)
          library(sna)
        })
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
          RhpcBLASctl::blas_set_num_threads(1L)
          RhpcBLASctl::omp_set_num_threads(1L)
        }
        Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
      })
      parallel::clusterExport(cl, c("pfit", "time_window", "PMF_mark", "cond_intensity",
                                   "formula_RHS", "truncation", "mark_decay", "growth_only",
                                   "max_node_time", "inhom_bg", "seed_net", "seed_times",
                                   "mu_multiplier"), envir = environment())
      cl
    }, error = function(e) {
      if (verbose) cat("  WARNING: PSOCK cluster failed, falling back to fork:", e$message, "\n")
      NULL
    })
  }
  message(sprintf("  [GOF] Memory before simulation: %.1f Mb", gc()[2, 2]))
  message(sprintf("  [GOF] Starting %d parallel simulations on %d %s at %s",
                  n_sim, n_workers, if (use_psock) "workers" else "cores", format(Sys.time(), "%H:%M:%S")))
  sim_results <- tryCatch({
    if (!is.null(cl_gof)) {
      parallel::parLapply(cl_gof, seq_len(n_sim), sim_fun)
    } else {
      safe_parallel_lapply(seq_len(n_sim), sim_fun, mc.cores = n_workers, parallel_type = "auto")
    }
  }, error = function(e) {
    if (verbose) cat("  ERROR: Failed to run simulations:", e$message, "\n")
    message(sprintf("  [GOF] Simulation FAILED: %s", e$message))
    list()
  })
  message(sprintf("  [GOF] Simulations complete at %s", format(Sys.time(), "%H:%M:%S")))
  
  if (is.null(sim_results) || length(sim_results) == 0) {
    if (verbose) cat("  No simulation results; returning empty GOF results\n")
    return(GOF_results)
  }
  
  sim_nets <- lapply(sim_results, function(x) x$net)
  sim_nets <- sim_nets[!sapply(sim_nets, is.null)]
  n_success <- length(sim_nets)
  n_fail <- n_sim - n_success
  GOF_results$nets_sim <- sim_nets
  
  # Report errors if any
  if (n_fail > 0 && verbose) {
    errors <- sapply(sim_results, function(x) if (!is.null(x$error)) x$error else NULL)
    errors <- errors[!sapply(errors, is.null)]
    if (length(errors) > 0) {
      cat("  GOF simulation errors (showing first 3):\n")
      for (i in seq_len(min(3, length(errors)))) {
        cat("    ", errors[[i]], "\n")
      }
    }
  }
  
  if (verbose) cat("  GOF simulations:", n_success, "succeeded,", n_fail, "failed\n")
  
  if (length(sim_nets) > 0) {
    if (verbose) cat("  Computing GOF statistics (parallelized)...\n")
    t_stats <- proc.time()
    
    # Observed statistics (single network)
    if (verbose) cat("    Computing observed statistics...\n")
    # Ensure vertex attributes are set if formula includes nodeMatch/nodeMix
    if (any(grepl("nodeMix|nodeMatch", formula_RHS))) {
      net_obs <- tryCatch({
        ensure_vertex_attribute(net_obs, "gender", default_value = "unknown")
      }, error = function(e) {
        if (verbose) cat("      Warning: Could not ensure gender attribute:", e$message, "\n")
        net_obs
      })
    }
    
    # Compute each observed statistic independently (failures don't stop others)
    GOF_results$degree_obs <- tryCatch({
      degree_dist(net_obs, max_deg, min_deg = degree)
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute observed degree distribution:", e$message, "\n")
      NULL
    })
    
    GOF_results$esp_obs <- tryCatch({
      esp_dist(net_obs, k_esp, min_esp = esp)
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute observed ESP distribution:", e$message, "\n")
      NULL
    })
    
    GOF_results$geodist_obs <- tryCatch({
      geodist_dist(net_obs)
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute observed geodesic distances:", e$message, "\n")
      NULL
    })
    
    GOF_results$wait_obs <- tryCatch({
      # Use a standard set of statistics for waiting times in ALL cases
      wait_formula_standard <- "edges + triangles + star(c(2,3))"
      wait_obs_raw <- waiting_times_between_formations(net_obs, formula_RHS = wait_formula_standard)
      
  # Standard statistic names for these waiting times
  stat_names_standard <- c("edges", "triangles", "star2", "star3")
  if (length(stat_names_standard) == length(wait_obs_raw)) {
    names(wait_obs_raw) <- stat_names_standard
  }
  wait_obs_raw
}, error = function(e) {
  if (verbose) cat("      Warning: Could not compute observed waiting times:", e$message, "\n")
  NULL
})

# Observed nodeMix statistics (if gender attribute exists)
# Fail gracefully: if nodal covariates cause errors (encoding, ernm C++, etc.), skip nodeMix entirely
if (verbose) cat("    Computing observed nodeMix statistics...\n")
tryCatch({
  if ("gender" %in% list.vertex.attributes(net_obs) || 
      any(grepl("nodeMix|nodeMatch", formula_RHS))) {
    net_obs_clean <- tryCatch(
      ensure_vertex_attribute(net_obs, "gender", default_value = "unknown"),
      error = function(e) {
        if (verbose) cat("      Warning: Could not ensure gender attribute:", e$message, "\n")
        NULL
      }
    )
    if (!is.null(net_obs_clean)) {
      GOF_results$nodemix_obs <- tryCatch(
        as.vector(calculateStatistics(net_obs_clean ~ nodeMix('gender'))),
        error = function(e) {
          if (verbose) cat("      Warning: Could not compute nodeMix statistics:", e$message, "\n")
          NULL
        }
      )
    }
  }
}, error = function(e) {
  if (verbose) cat("      Warning: Nodal covariate (nodeMix) skipped:", e$message, "\n")
  GOF_results$nodemix_obs <- NULL
})

# --- Optimization: Clean up net_obs before parallel stats if possible ---
# (Actually we need it for names, but we can clear some memory)
gc()

# Simulated statistics (parallelized) - each wrapped in tryCatch
n_deg_bins <- max_deg - degree + 1L
n_esp_bins <- k_esp - esp + 1L

# Pre-calculate common objects for parallel workers to reduce overhead
# We can't easily export the ERNM model object itself, but we can ensure
# the workers are as lean as possible.

if (verbose) cat("    Computing distributional statistics (Degree, ESP, Geodist, nodeMix)...\n")
message(sprintf("  [GOF] Starting distributional stats (%d nets, %d %s) at %s",
                length(sim_nets), n_workers, if (use_psock) "workers" else "cores", format(Sys.time(), "%H:%M:%S")))

# Combined distributional statistics: Degree, ESP, Geodist, nodeMix
# Pre-calculate observed nodeMix presence to avoid repeated grepl/list.vertex.attributes
has_nodemix_obs <- !is.null(GOF_results$nodemix_obs)
needs_gender <- any(grepl("nodeMix|nodeMatch", formula_RHS))

dist_stats_fun <- function(n) {
    tryCatch({
      if (is.null(n)) return(NULL)
      n_clean <- n
      if (needs_gender) {
        n_clean <- tryCatch({
          ensure_vertex_attribute(n, "gender", default_value = "unknown")
        }, error = function(e) n)
      }
      nodemix_val <- NULL
      if (has_nodemix_obs) {
        nodemix_val <- tryCatch({
          as.vector(calculateStatistics(n_clean ~ nodeMix('gender')))
        }, error = function(e) NULL)
      }
      list(
        degree = degree_dist(n, max_deg, min_deg = degree),
        esp = esp_dist(n, k_esp, min_esp = esp),
        geodist = geodist_dist(n),
        nodemix = nodemix_val
      )
    }, error = function(e) {
      list(degree = rep(NA_real_, n_deg_bins), 
           esp = rep(NA_real_, n_esp_bins), 
           geodist = numeric(0), 
           nodemix = if (has_nodemix_obs) rep(NA_real_, length(GOF_results$nodemix_obs)) else NULL)
    })
}
dist_stats_sim <- tryCatch({
  if (!is.null(cl_gof)) {
    parallel::clusterExport(cl_gof, c("max_deg", "degree", "k_esp", "esp", "has_nodemix_obs", "needs_gender",
                                     "n_deg_bins", "n_esp_bins", "GOF_results"), envir = environment())
    parallel::clusterExport(cl_gof, c("degree_dist", "esp_dist", "geodist_dist", "ensure_vertex_attribute"),
                            envir = asNamespace("hawkesNet"))
    parallel::parLapply(cl_gof, sim_nets, dist_stats_fun)
  } else {
    safe_parallel_lapply(sim_nets, dist_stats_fun, mc.cores = n_workers, parallel_type = "auto")
  }
}, error = function(e) {
  if (verbose) cat("      Warning: Parallel distributional stats failed:", e$message, "\n")
  NULL
})
    
    if (!is.null(dist_stats_sim)) {
      # Filter out NULLs and atomic vectors (error messages) if any task failed
      dist_stats_sim <- dist_stats_sim[sapply(dist_stats_sim, is.list)]
      
      if (length(dist_stats_sim) > 0) {
        GOF_results$degree_sim <- do.call(rbind, lapply(dist_stats_sim, function(x) {
          if (is.list(x) && !is.null(x$degree)) x$degree else rep(NA_real_, n_deg_bins)
        }))
        GOF_results$esp_sim <- do.call(rbind, lapply(dist_stats_sim, function(x) {
          if (is.list(x) && !is.null(x$esp)) x$esp else rep(NA_real_, n_esp_bins)
        }))
        GOF_results$geodist_sim <- lapply(dist_stats_sim, function(x) {
          if (is.list(x) && !is.null(x$geodist)) x$geodist else numeric(0)
        })
        nodemix_list <- lapply(dist_stats_sim, function(x) {
          if (is.list(x)) x$nodemix else NULL
        })
        if (!all(sapply(nodemix_list, is.null))) {
          # Filter out NULLs from nodemix_list before rbind
          nodemix_list_clean <- nodemix_list[!sapply(nodemix_list, is.null)]
          if (length(nodemix_list_clean) > 0) {
            # Check if all elements are numeric vectors of correct length
            if (has_nodemix_obs) {
              expected_len <- length(GOF_results$nodemix_obs)
              nodemix_list_clean <- lapply(nodemix_list_clean, function(x) {
                if (is.numeric(x) && length(x) == expected_len) x else rep(NA_real_, expected_len)
              })
            }
            GOF_results$nodemix_sim <- do.call(rbind, nodemix_list_clean)
          }
        }
      }
      rm(dist_stats_sim); gc()
    }

    # --- Fallback: ensure core stats exist for plotting ---
    # In rare cases (e.g. parallel worker failure), dist_stats_sim can be NULL/empty even when
    # simulations succeeded. Ensure degree/ESP/geodesic summaries are populated so plots aren't empty.
    if (length(sim_nets) > 0) {
      if (is.null(GOF_results$degree_sim)) {
        GOF_results$degree_sim <- tryCatch({
          do.call(rbind, lapply(sim_nets, function(n) degree_dist(n, max_deg, min_deg = degree)))
        }, error = function(e) NULL)
      }
      if (is.null(GOF_results$esp_sim)) {
        GOF_results$esp_sim <- tryCatch({
          do.call(rbind, lapply(sim_nets, function(n) esp_dist(n, k_esp, min_esp = esp)))
        }, error = function(e) NULL)
      }
      if (is.null(GOF_results$geodist_sim)) {
        GOF_results$geodist_sim <- tryCatch({
          lapply(sim_nets, function(n) geodist_dist(n))
        }, error = function(e) NULL)
      }
    }
    
    if (verbose) cat("    Computing waiting times (expensive)...\n")
    message(sprintf("  [GOF] Starting waiting times (%d nets, %d %s) at %s",
                    length(sim_nets), n_workers, if (use_psock) "workers" else "cores", format(Sys.time(), "%H:%M:%S")))
    stat_names_standard <- c("edges", "triangles", "star2", "star3")
    wait_fun <- function(n) {
      tryCatch({
        wait_sim_raw <- waiting_times_between_formations(n, formula_RHS = "edges + triangles + star(c(2,3))")
        if (length(stat_names_standard) == length(wait_sim_raw)) {
          names(wait_sim_raw) <- stat_names_standard
        }
        wait_sim_raw
      }, error = function(e) list())
    }
    GOF_results$wait_sim <- tryCatch({
      if (!is.null(cl_gof)) {
        parallel::clusterExport(cl_gof, c("stat_names_standard"), envir = environment())
        parallel::clusterExport(cl_gof, c("waiting_times_between_formations"), envir = asNamespace("hawkesNet"))
        parallel::parLapply(cl_gof, sim_nets, wait_fun)
      } else {
        safe_parallel_lapply(sim_nets, wait_fun, mc.cores = n_workers, parallel_type = "auto")
      }
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute simulated waiting times:", e$message, "\n")
      NULL
    })
    
    message(sprintf("  [GOF] Waiting times complete at %s", format(Sys.time(), "%H:%M:%S")))
    if (verbose) {
      cat("    Waiting times: done\n")
      cat("  GOF statistics:", round((proc.time() - t_stats)[3], 1), "s\n")
    }
    
    # Generate plots (if ggplot2 is available)
    if (verbose) cat("  Generating GOF plots...\n")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      GOF_results$plots <- tryCatch(
        create_gof_plots(GOF_results),
        error = function(e) {
          if (verbose) cat("    Warning: Could not generate plots:", e$message, "\n")
          list()
        }
      )
    } else {
      if (verbose) cat("    ggplot2 not available; skipping plots\n")
      GOF_results$plots <- list()
    }
  } else {
    if (verbose) cat("  No successful simulations; returning empty GOF results\n")
    GOF_results$plots <- list()  # ensure plots is always set (never NULL)
  }
  
  if (!is.null(cl_gof)) {
    tryCatch(parallel::stopCluster(cl_gof), error = function(e) NULL)
    cl_gof <- NULL
  }
  
  # Always return results, even if some computations failed
  GOF_results
}

#' Create GOF plots for distributional fits and waiting times.
#'
#' @param GOF_results List with observed and simulated statistics from \code{gof()}.
#' @return List of ggplot objects: \code{degree_plot}, \code{esp_plot}, \code{geodist_plot}, 
#'   \code{nodemix_plot}, \code{waiting_times_plot}.
#' @noRd
create_gof_plots <- function(GOF_results) {
  plots <- list()
  
  # Helper function to create boxplot for distributional statistics
  # Observed values shown as dots, simulated as boxplots.
  # x_values: optional (e.g. 0:(n-1)) so x-axis starts at 0 for degree/ESP.
  create_dist_plot <- function(obs, sim, stat_name, x_label = NULL, x_values = NULL) {
    if (is.null(obs) || is.null(sim) || nrow(sim) == 0) return(NULL)
    
    n_obs <- length(obs)
    n_sim <- nrow(sim)
    x_vals <- if (!is.null(x_values) && length(x_values) == n_obs) x_values else seq_len(n_obs)
    
    # Create data frame for simulated (boxplot)
    df_sim <- data.frame(
      value = as.vector(sim),
      x = rep(x_vals, each = nrow(sim)),
      type = "Simulated"
    )
    
    # Create data frame for observed (dots)
    df_obs <- data.frame(
      value = obs,
      x = x_vals,
      type = "Observed"
    )
    
    x_lab <- if (is.null(x_label)) "Index" else x_label
    x_levels <- sort(unique(x_vals))
    p <- ggplot2::ggplot(df_sim, ggplot2::aes(x = factor(x, levels = x_levels), y = value)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5, fill = "#56B4E9") +
      ggplot2::geom_point(data = df_obs, ggplot2::aes(x = factor(x, levels = x_levels), y = value),
                          color = "#E69F00", size = 2, shape = 19) +
      ggplot2::labs(
        title = paste(stat_name, "Distribution"),
        x = x_lab,
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
    
    p
  }
  
  # Degree distribution plot (x-axis starts at 0)
  if (!is.null(GOF_results$degree_obs) && !is.null(GOF_results$degree_sim)) {
    n_deg <- length(GOF_results$degree_obs)
    plots$degree_plot <- create_dist_plot(
      GOF_results$degree_obs,
      GOF_results$degree_sim,
      "Degree",
      "Degree",
      x_values = seq(0L, length.out = n_deg)
    )
  }
  
  # ESP distribution plot (x-axis starts at 0)
  if (!is.null(GOF_results$esp_obs) && !is.null(GOF_results$esp_sim)) {
    n_esp <- length(GOF_results$esp_obs)
    plots$esp_plot <- create_dist_plot(
      GOF_results$esp_obs,
      GOF_results$esp_sim,
      "ESP (Edge-wise Shared Partners)",
      "ESP Count",
      x_values = seq(0L, length.out = n_esp)
    )
  }
  
  # Geodesic distance plot (boxplot with relative proportions)
  if (!is.null(GOF_results$geodist_obs) && !is.null(GOF_results$geodist_sim) && 
      length(GOF_results$geodist_obs) > 0 && length(GOF_results$geodist_sim) > 0) {
    # Calculate relative proportions for each distance
    # Observed: proportion of pairs at each distance
    obs_dist <- GOF_results$geodist_obs
    obs_tab <- table(factor(obs_dist, levels = sort(unique(c(obs_dist, unlist(GOF_results$geodist_sim))))))
    obs_prop <- as.numeric(obs_tab) / sum(obs_tab)
    obs_dist_levels <- as.numeric(names(obs_tab))
    
    # Simulated: proportion for each simulation, then average
    sim_dist_list <- GOF_results$geodist_sim
    all_dist_levels <- sort(unique(c(obs_dist, unlist(sim_dist_list))))
    
    sim_prop_mat <- do.call(rbind, lapply(sim_dist_list, function(sim_dist) {
      if (length(sim_dist) == 0) return(rep(0, length(all_dist_levels)))
      sim_tab <- table(factor(sim_dist, levels = all_dist_levels))
      as.numeric(sim_tab) / sum(sim_tab)
    }))
    
    # Create data frame
    df_geod <- data.frame(
      distance = rep(all_dist_levels, 2),
      proportion = c(
        obs_prop[match(all_dist_levels, obs_dist_levels)],
        colMeans(sim_prop_mat, na.rm = TRUE)
      ),
      type = rep(c("Observed", "Simulated"), each = length(all_dist_levels))
    )
    df_geod$proportion[is.na(df_geod$proportion)] <- 0
    
    # Observed as dot, simulated as boxplot
    df_sim_box <- data.frame(
      distance = rep(all_dist_levels, each = nrow(sim_prop_mat)),
      proportion = as.vector(sim_prop_mat)
    )
    df_obs_dot <- df_geod[df_geod$type == "Observed", ]
    
    plots$geodist_plot <- ggplot2::ggplot(df_sim_box, ggplot2::aes(x = factor(distance), y = proportion)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5, fill = "#56B4E9") +
      ggplot2::geom_point(data = df_obs_dot, ggplot2::aes(x = factor(distance), y = proportion),
                         color = "#E69F00", size = 2, shape = 19) +
      ggplot2::labs(
        title = "Geodesic Distance Distribution (Relative Proportions)",
        x = "Geodesic Distance",
        y = "Proportion of Pairs"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }
  
  # nodeMix plot
  if (!is.null(GOF_results$nodemix_obs) && !is.null(GOF_results$nodemix_sim) && 
      nrow(GOF_results$nodemix_sim) > 0) {
    plots$nodemix_plot <- create_dist_plot(
      GOF_results$nodemix_obs,
      GOF_results$nodemix_sim,
      "nodeMix (Gender Mixing)",
      "Mixing Pattern"
    )
  }
  
  # Waiting times plot (Boxplots of waiting times by statistic and type)
  if (!is.null(GOF_results$wait_obs) && !is.null(GOF_results$wait_sim) && 
      length(GOF_results$wait_sim) > 0 && length(GOF_results$wait_obs) > 0) {
    # Extract waiting times for each statistic
    wait_names <- names(GOF_results$wait_obs)
    if (is.null(wait_names) || length(wait_names) == 0) {
      wait_names <- paste0("stat", seq_along(GOF_results$wait_obs))
    }
    
    df_wait_list <- list()
    for (i in seq_along(wait_names)) {
      stat_name <- wait_names[i]
      obs_wait <- GOF_results$wait_obs[[i]]
      if (is.null(obs_wait)) obs_wait <- numeric(0)
      
      # For simulated, we want to keep the individual simulation IDs to allow boxplots
      # if we had multiple simulations. However, wait_sim is a list of lists.
      # Let's flatten it but keep the simulation ID.
      sim_wait_df <- do.call(rbind, lapply(seq_along(GOF_results$wait_sim), function(sim_id) {
        x <- GOF_results$wait_sim[[sim_id]]
        vals <- if (is.list(x) && stat_name %in% names(x)) {
          x[[stat_name]]
        } else if (is.list(x) && i <= length(x)) {
          x[[i]]
        } else {
          NULL
        }
        if (length(vals) > 0) {
          data.frame(waiting_time = vals, sim_id = sim_id)
        } else {
          NULL
        }
      }))
      
      if (length(obs_wait) > 0 || (!is.null(sim_wait_df) && nrow(sim_wait_df) > 0)) {
        if (length(obs_wait) > 0) {
          df_obs <- data.frame(
            waiting_time = obs_wait,
            type = "Observed",
            statistic = stat_name
          )
          df_wait_list[[paste0(stat_name, "_obs")]] <- df_obs
        }
        if (!is.null(sim_wait_df) && nrow(sim_wait_df) > 0) {
          df_sim <- data.frame(
            waiting_time = sim_wait_df$waiting_time,
            type = "Simulated",
            statistic = stat_name
          )
          df_wait_list[[paste0(stat_name, "_sim")]] <- df_sim
        }
      }
    }
    
    if (length(df_wait_list) > 0) {
      df_wait <- do.call(rbind, df_wait_list)
      
      plots$waiting_times_plot <- ggplot2::ggplot(df_wait, ggplot2::aes(x = type, y = waiting_time, fill = type)) +
        ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
        ggplot2::scale_fill_manual(values = c("Observed" = "#E69F00", "Simulated" = "#56B4E9")) +
        ggplot2::scale_y_log10() +
        ggplot2::facet_wrap(~statistic, scales = "free_y") +
        ggplot2::labs(
          title = "Waiting Times Between Structure Formations",
          x = "",
          y = "Waiting Time (log scale)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          strip.text = ggplot2::element_text(face = "bold")
        )
    }
  }
  
  plots
}
