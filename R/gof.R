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
  degs <- sna::degree(net, gmode = "graph")
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
    as.vector(ernm::calculateStatistics(net ~ esp(min_esp:k_max)))
  }, error = function(e) rep(NA_real_, k_max - min_esp + 1))
}

#' Geodesic distance distribution (upper triangle of distance matrix, excluding Inf).
#'
#' @param net Network object.
#' @return Numeric vector of geodesic distances (excluding Inf and self-distances).
#' @noRd
geodist_dist <- function(net) {
  d <- sna::geodist(net, inf.replace = NA)
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
  nv <- network::network.size(net)
  if (nv == 0) return(net)
  
  # Remove problematic 'na' attribute if it exists
  if ("na" %in% network::list.vertex.attributes(net)) {
    network::delete.vertex.attribute(net, "na")
  }
  
  # Check if attribute exists
  if (!attr_name %in% network::list.vertex.attributes(net)) {
    # Attribute doesn't exist: set all nodes to default
    network::set.vertex.attribute(net, attr_name, rep(default_value, nv))
  } else {
    # Attribute exists: check and fix missing/NA values
    attr_vals <- network::get.vertex.attribute(net, attr_name)
    if (length(attr_vals) < nv) {
      # Not enough values: pad with default
      attr_vals <- c(attr_vals, rep(default_value, nv - length(attr_vals)))
      network::set.vertex.attribute(net, attr_name, attr_vals)
    } else if (any(is.na(attr_vals)) || any(attr_vals == "")) {
      # Has NA or empty values: replace with default
      attr_vals[is.na(attr_vals) | attr_vals == ""] <- default_value
      network::set.vertex.attribute(net, attr_name, attr_vals)
    }
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
  el <- network::as.edgelist(net)
  
  # Determine number of statistics by creating a temporary model
  n <- network::network.size(net)
  g_temp <- network::network.initialize(n, directed = network::is.directed(net))
  
  # Copy vertex attributes from net to g_temp (required for nodeMatch/nodeMix)
  vattrs <- network::list.vertex.attributes(net)
  for (attr_name in vattrs) {
    if (attr_name != "na") {
      attr_vals <- network::get.vertex.attribute(net, attr_name)
      if (length(attr_vals) == n) {
        network::set.vertex.attribute(g_temp, attr_name, attr_vals)
      }
    }
  }
  if ("na" %in% network::list.vertex.attributes(g_temp)) {
    network::delete.vertex.attribute(g_temp, "na")
  }
  
  # Create formula with g_temp in the environment
  formula_str <- paste0("g_temp ~ ", formula_RHS)
  formula_obj <- as.formula(formula_str)
  model_temp <- ernm::createCppModel(formula_obj)
  model_temp$setNetwork(ernm::as.BinaryNet(g_temp))
  model_temp$calculate()
  n_stats <- length(model_temp$statistics())
  rm(model_temp, g_temp)
  
  if (nrow(el) == 0) {
    # Return empty structure matching formula
    return(setNames(rep(list(numeric()), n_stats), paste0("stat", seq_len(n_stats))))
  }

  # Try edge times first; fall back to vertex times for edge ordering
  edge_times <- NULL
  if (time_attr %in% network::list.edge.attributes(net)) {
    edge_times <- network::get.edge.attribute(net, time_attr)
  }
  if (is.null(edge_times) || length(edge_times) != nrow(el)) {
    # Use vertex times: assign each edge the max time of its endpoints
    vtimes <- network::get.vertex.attribute(net, time_attr)
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
  g <- network::network.initialize(n, directed = network::is.directed(net))
  
  # Copy vertex attributes from net to g (required for nodeMatch/nodeMix)
  vattrs <- network::list.vertex.attributes(net)
  for (attr_name in vattrs) {
    if (attr_name != "na") {
      attr_vals <- network::get.vertex.attribute(net, attr_name)
      if (length(attr_vals) == n) {
        network::set.vertex.attribute(g, attr_name, attr_vals)
      }
    }
  }
  if ("na" %in% network::list.vertex.attributes(g)) {
    network::delete.vertex.attribute(g, "na")
  }

  # Initialize tracking for all statistics
  stat_prev <- rep(0, n_stats)
  stat_times <- rep(list(numeric()), n_stats)

  # Create ERNM model once and initialize with empty network
  # Create formula with g now that it exists
  formula_str <- paste0("g ~ ", formula_RHS)
  formula_obj <- as.formula(formula_str)
  model <- ernm::createCppModel(formula_obj)
  model$setNetwork(ernm::as.BinaryNet(g))
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
        network::add.edges(g, tail = tail_j, head = head_j)
        
        # Update model to reflect this new edge (for next iteration's change stats)
        # OPTIMIZATION: Use model$toggle() if available in ERNM, otherwise calculate()
        # ERNM's computeChangeStats + manual add.edges + calculate is standard.
        model$setNetwork(ernm::as.BinaryNet(g))
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
#' @param cores Number of cores for parallelization (default 7).
#' @param max_deg Maximum degree for degree distribution (default 15).
#' @param k_esp Maximum ESP count for ESP distribution (default 15).
#' @param degree Minimum degree to include in degree distribution (default 0).
#' @param esp Minimum ESP to include in ESP distribution (default 0).
#' @param mu_multiplier Multiplier for mu in simulations (default 5).
#' @param verbose Print progress messages (default TRUE).
#' @return List with observed and simulated statistics and plots:
#'   \itemize{
#'     \item \code{degree_obs}: Observed degree distribution
#'     \item \code{degree_sim}: Matrix of simulated degree distributions (rows = simulations)
#'     \item \code{esp_obs}: Observed ESP distribution
#'     \item \code{esp_sim}: Matrix of simulated ESP distributions (rows = simulations)
#'     \item \code{geodist_obs}: Observed geodesic distances
#'     \item \code{geodist_sim}: List of simulated geodesic distances
#'     \item \code{wait_obs}: Observed waiting times between structure formations
#'     \item \code{wait_sim}: List of simulated waiting times
#'     \item \code{nodemix_obs}: Observed nodeMix statistics
#'     \item \code{nodemix_sim}: Matrix of simulated nodeMix statistics (rows = simulations)
#'     \item \code{nets_sim}: List of simulated network objects (for further comparison)
#'     \item \code{plots}: List of ggplot objects (if ggplot2 available):
#'       \itemize{
#'         \item \code{degree_plot}: Boxplot of degree distributions
#'         \item \code{esp_plot}: Boxplot of ESP distributions
#'         \item \code{geodist_plot}: Histogram of geodesic distances
#'         \item \code{nodemix_plot}: Boxplot of nodeMix statistics
#'         \item \code{waiting_times_plot}: Faceted histograms of waiting times
#'       }
#'   }
#' @export
gof <- function(fit, net_obs, params_init, PMF_mark, cond_intensity, formula_RHS,
                time_window = c(0, 0.05), truncation = 100L, mark_decay = "activity",
                max_node_time = 1, inhom_bg = NULL, n_sim = 50L, cores = 7L,
                max_deg = 15L, k_esp = 15L, degree = 0L, esp = 0L, mu_multiplier = 5, verbose = TRUE) {
  
  # Initialize results early (will be populated even if some computations fail)
  GOF_results <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                      geodist_obs = NULL, geodist_sim = NULL, wait_obs = NULL, wait_sim = NULL,
                      nodemix_obs = NULL, nodemix_sim = NULL, nets_sim = NULL, plots = list())
  
  if (is.null(fit)) {
    if (verbose) cat("  No fit available; skipping GOF\n")
    return(GOF_results)
  }
  
  if (verbose) cat("  Simulating", n_sim, "networks from fitted model...\n")
  
  # Reconstruct fitted params from fit$par.
  # fit_hawkesNet strips fixed_params and vertex_categorical_levels from the

  # skeleton before optim, so fit$par only contains the FREE parameters.
  # We must strip the same fields from our skeleton before relisting, then
  # merge the fitted values back into the full params_init structure.
  
  fixed <- fit$fixed_params
  if (verbose && !is.null(fixed)) cat("  Fixed parameters:", paste(fixed, collapse = ", "), "\n")
  
  # Build skeleton matching exactly what optim saw (no fixed params, no levels)
  skel <- params_init
  if (!is.null(fixed)) {
    for (p in fixed) skel[[p]] <- NULL
  }
  skel$vertex_categorical_levels <- NULL
  
  pfit_vals <- tryCatch({
    relist(fit$fit$par, skeleton = skel)
  }, error = function(e) {
    if (verbose) cat("  ERROR: relist failed:", e$message, "\n")
    if (verbose) cat("    fit$par length:", length(fit$fit$par), "| skeleton length:", length(unlist(skel)), "\n")
    if (verbose) cat("    fit$par names:", paste(names(fit$fit$par), collapse = ", "), "\n")
    if (verbose) cat("    skeleton names:", paste(names(unlist(skel)), collapse = ", "), "\n")
    if (verbose) cat("  FALLING BACK to params_init — GOF simulations may be wrong!\n")
    NULL
  })

  # Merge fitted values back into full parameter structure
  pfit <- params_init
  if (!is.null(pfit_vals)) {
    for (n in names(pfit_vals)) pfit[[n]] <- pfit_vals[[n]]
  } else {
    if (verbose) cat("  WARNING: Using params_init as fallback (relist failed)\n")
  }

  # Restore metadata
  pfit$vertex_categorical_levels <- params_init$vertex_categorical_levels
  
  # Handle inhomogeneous background
  use_inhom <- !is.null(inhom_bg) && !is.null(inhom_bg$mu_fit) && !is.null(inhom_bg$mu_fit$mu_fun)
  
  if (use_inhom) {
    # For inhomogeneous: mu is not used directly, but we set it for compatibility
    # The actual mu_at_t will be computed from mu_fun during simulation
    Tval <- time_window[2] - time_window[1]
    pfit$mu <- inhom_bg$integral_bg / Tval  # Average mu for compatibility
    if (verbose) {
      cat("  Using inhomogeneous background for simulations (matches fitted model)\n")
    }
  }
  # If mu is fixed, it stays at params_init$mu (already in pfit from the merge above)
  # If mu is free, the fitted value is already in pfit from pfit_vals
  
  # Ensure parameters are within valid ranges
  pfit$K <- min(max(pfit$K, 0.001), 0.999)  # K must be in (0,1) for stability
  pfit$mu <- max(pfit$mu, 0.001)  # mu must be positive
  pfit$node_lambda <- max(pfit$node_lambda, 0.1)  # node_lambda must be positive
  pfit$beta_overall <- max(pfit$beta_overall, 0.001)
  pfit$beta_edges <- max(pfit$beta_edges, 0.001)
  
  # Log the actual parameters being used for GOF simulations
  if (verbose) {
    cat("  GOF simulation parameters:\n")
    scalar_params <- c("mu", "beta_overall", "K", "beta_edges", "node_lambda")
    for (p in scalar_params) {
      if (!is.null(pfit[[p]])) {
        src <- if (!is.null(fixed) && p %in% fixed) "(fixed)" else "(fitted)"
        cat(sprintf("    %s = %.6f %s\n", p, pfit[[p]], src))
      }
    }
    if (!is.null(pfit$CS_params)) {
      cat("    CS_params =", paste(round(pfit$CS_params, 4), collapse = ", "), "(fitted)\n")
    }
  }
  
  # Set vertex_categorical if it exists (use defaults if needed)
  if (!is.null(params_init$vertex_categorical)) {
    if (is.null(pfit$vertex_categorical)) {
      pfit$vertex_categorical <- params_init$vertex_categorical
    } else {
      # Repair vertex_categorical parameters (may be invalid from optimization)
      pfit <- tryCatch({
        repair_vertex_categorical_params(pfit, eps = 1e-6)
      }, error = function(e) {
        if (verbose) cat("  WARNING: Failed to repair vertex_categorical:", e$message, "\n")
        pfit
      })
    }
  }
  
  # Validate parameters before simulation
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
    cat("  Using", cores, "cores for parallel GOF simulations...\n")
    if (use_inhom) {
      cat("  Simulations will use inhomogeneous background (cond_intensity_inhom)\n")
    } else {
      cat("  Simulations will use homogeneous background (cond_intensity)\n")
    }
  }
  # CRITICAL: Use safe_parallel_lapply instead of raw parallel::mclapply.
  # Raw mclapply doesn't guard BLAS/OpenMP threads before forking and doesn't
  # clean up zombie children from prior mclapply calls (e.g. the fitting step).
  # After sequential fits, BLAS threads may be restored to multi-threaded state;
  # forking with active BLAS threads causes deadlock on the second fork.
  cat(sprintf("  [GOF] Memory before simulation fork: %.1f Mb\n", gc()[2, 2]), file = stderr())
  cat(sprintf("  [GOF] Starting %d parallel simulations on %d cores at %s\n",
              n_sim, cores, format(Sys.time(), "%H:%M:%S")), file = stderr())
  sim_results <- tryCatch({
    safe_parallel_lapply(seq_len(n_sim), function(i) {
      s <- tryCatch(
        sim_hawkesNet(
          params = pfit,
          time_window = time_window,
          PMF_mark = PMF_mark,
          cond_intensity = cond_intensity,
          formula_RHS = formula_RHS,
          truncation = truncation,
          mark_decay = mark_decay,
          max_node_time = max_node_time,
          hashed_edges = TRUE,
          verbose = FALSE,
          mu_multiplier = mu_multiplier,
          stop_on_full_network = FALSE,
          inhom_bg = inhom_bg
        ),
        error = function(e) { 
          return(list(net = NULL, error = paste0("Sim ", i, ": ", e$message))) 
        }
      )
      if (is.null(s$net)) {
        return(list(net = NULL, error = ifelse(is.null(s$error), paste0("Sim ", i, ": unknown error"), s$error)))
      }
      return(list(net = s$net, error = NULL))
    }, mc.cores = cores, parallel_type = "auto")
  }, error = function(e) {
    if (verbose) cat("  ERROR: Failed to run simulations:", e$message, "\n")
    cat(sprintf("  [GOF] Simulation FAILED: %s\n", e$message), file = stderr())
    list()  # Return empty list if simulations fail completely
  })
  cat(sprintf("  [GOF] Simulations complete at %s\n", format(Sys.time(), "%H:%M:%S")), file = stderr())
  
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
      wait_obs_raw <- waiting_times_between_formations(net_obs, formula_RHS = formula_RHS)
      # Get ERNM statistic names from formula
      exp_cs <- tryCatch({
        expected_params_PMF_mark_CS(net_obs, formula_RHS)
      }, error = function(e) NULL)
      if (!is.null(exp_cs) && !is.null(exp_cs$CS_params_names) && 
          length(exp_cs$CS_params_names) == length(wait_obs_raw)) {
        names(wait_obs_raw) <- exp_cs$CS_params_names
      }
      wait_obs_raw
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute observed waiting times:", e$message, "\n")
      NULL
    })
    
    # Observed nodeMix statistics (if gender attribute exists)
    if (verbose) cat("    Computing observed nodeMix statistics...\n")
    tryCatch({
      if ("gender" %in% network::list.vertex.attributes(net_obs) || 
          any(grepl("nodeMix|nodeMatch", formula_RHS))) {
        # Ensure gender attribute is properly set for all nodes before ERNM operations
        net_obs_clean <- ensure_vertex_attribute(net_obs, "gender", default_value = "unknown")
        GOF_results$nodemix_obs <- as.vector(ernm::calculateStatistics(net_obs_clean ~ nodeMix('gender')))
      }
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute observed nodeMix:", e$message, "\n")
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
    cat(sprintf("  [GOF] Starting distributional stats (%d nets, %d cores) at %s\n",
                length(sim_nets), cores, format(Sys.time(), "%H:%M:%S")), file = stderr())
    
    # Combined distributional statistics: Degree, ESP, Geodist, nodeMix
    if (verbose) cat("    Computing distributional statistics (Degree, ESP, Geodist, nodeMix)...\n")
    cat(sprintf("  [GOF] Starting distributional stats (%d nets, %d cores) at %s\n",
                length(sim_nets), cores, format(Sys.time(), "%H:%M:%S")), file = stderr())
    
    # Pre-calculate observed nodeMix presence to avoid repeated grepl/list.vertex.attributes
    has_nodemix_obs <- !is.null(GOF_results$nodemix_obs)
    needs_gender <- any(grepl("nodeMix|nodeMatch", formula_RHS))

    dist_stats_sim <- tryCatch({
      safe_parallel_lapply(sim_nets, function(n) {
        tryCatch({
          if (is.null(n)) return(NULL)
          # Ensure vertex attributes for nodeMix if needed
          n_clean <- n
          if (needs_gender) {
            n_clean <- tryCatch(ensure_vertex_attribute(n, "gender", default_value = "unknown"), 
                               error = function(e) n)
          }
          
          list(
            degree = degree_dist(n, max_deg, min_deg = degree),
            esp = esp_dist(n, k_esp, min_esp = esp),
            geodist = geodist_dist(n),
            nodemix = if (has_nodemix_obs) {
              as.vector(ernm::calculateStatistics(n_clean ~ nodeMix('gender')))
            } else NULL
          )
        }, error = function(e) {
          # Return a structure with NAs so rbind doesn't fail on atomic vectors
          list(degree = rep(NA_real_, n_deg_bins), 
               esp = rep(NA_real_, n_esp_bins), 
               geodist = numeric(0), 
               nodemix = if (has_nodemix_obs) rep(NA_real_, length(GOF_results$nodemix_obs)) else NULL)
        })
      }, mc.cores = cores, parallel_type = "auto")
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
    
    if (verbose) cat("    Computing waiting times (expensive)...\n")
    cat(sprintf("  [GOF] Starting waiting times (%d nets, %d cores) at %s\n",
                length(sim_nets), cores, format(Sys.time(), "%H:%M:%S")), file = stderr())
    GOF_results$wait_sim <- tryCatch({
      # Get ERNM statistic names from formula (use first simulated network or observed)
      exp_cs <- tryCatch({
        net_for_names <- if (!is.null(net_obs)) net_obs else if (length(sim_nets) > 0 && !is.null(sim_nets[[1]])) sim_nets[[1]] else NULL
        if (!is.null(net_for_names)) {
          expected_params_PMF_mark_CS(net_for_names, formula_RHS)
        } else NULL
      }, error = function(e) NULL)
      stat_names <- if (!is.null(exp_cs) && !is.null(exp_cs$CS_params_names)) exp_cs$CS_params_names else NULL
      
      # Use safe_parallel_lapply for stability
      safe_parallel_lapply(sim_nets, function(n) {
        tryCatch({
          wait_sim_raw <- waiting_times_between_formations(n, formula_RHS = formula_RHS)
          if (!is.null(stat_names) && length(stat_names) == length(wait_sim_raw)) {
            names(wait_sim_raw) <- stat_names
          }
          wait_sim_raw
        }, error = function(e) list())
      }, mc.cores = cores, parallel_type = "auto")
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute simulated waiting times:", e$message, "\n")
      NULL
    })
    
    cat(sprintf("  [GOF] Waiting times complete at %s\n", format(Sys.time(), "%H:%M:%S")), file = stderr())
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
        y = "Count",
        fill = "Type"
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
        y = "Proportion of Pairs",
        fill = "Type"
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
  
  # Waiting times plot (2-column format: observed in one column, simulated in another, row by row)
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
      
      sim_wait <- unlist(lapply(GOF_results$wait_sim, function(x) {
        if (is.list(x) && stat_name %in% names(x) && length(x[[stat_name]]) > 0) {
          x[[stat_name]]
        } else if (is.list(x) && i <= length(x) && length(x[[i]]) > 0) {
          x[[i]]
        } else {
          NULL
        }
      }))
      if (is.null(sim_wait)) sim_wait <- numeric(0)
      
      if (length(obs_wait) > 0 || length(sim_wait) > 0) {
        # Create separate data frames for observed and simulated
        df_obs <- data.frame(
          waiting_time = obs_wait,
          type = "Observed",
          statistic = stat_name
        )
        df_sim <- data.frame(
          waiting_time = sim_wait,
          type = "Simulated",
          statistic = stat_name
        )
        df_wait_list[[stat_name]] <- list(obs = df_obs, sim = df_sim)
      }
    }
    
    if (length(df_wait_list) > 0) {
      # Combine all observed and simulated separately
      df_obs_all <- do.call(rbind, lapply(df_wait_list, function(x) x$obs))
      df_sim_all <- do.call(rbind, lapply(df_wait_list, function(x) x$sim))
      
      # Create combined data frame with type and statistic
      df_wait <- rbind(df_obs_all, df_sim_all)
      
      # Common x-axis range so all panels use the same waiting-time scale
      x_range <- range(df_wait$waiting_time, na.rm = TRUE, finite = TRUE)
      if (diff(x_range) < .Machine$double.eps) x_range <- x_range + c(-0.5, 0.5)
      # facet_wrap so each panel (statistic x type) has its own y-scale; otherwise Observed
      # density spike compresses Simulated in the same row when using facet_grid.
      plots$waiting_times_plot <- ggplot2::ggplot(df_wait, ggplot2::aes(x = waiting_time, fill = type)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                alpha = 0.7, bins = 30, position = "identity") +
        ggplot2::scale_fill_manual(values = c("Observed" = "#E69F00", "Simulated" = "#56B4E9")) +
        ggplot2::facet_wrap(ggplot2::vars(statistic, type), scales = "free_y", ncol = 2L) +
        ggplot2::coord_cartesian(xlim = x_range) +
        ggplot2::labs(
          title = "Waiting Times Between Structure Formations",
          x = "Waiting Time",
          y = "Density"
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
