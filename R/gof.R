# =============================================================================
# Goodness-of-Fit (GOF) functions for hawkesGrowthNet
# =============================================================================

#' Degree distribution as vector of counts (degree 0, 1, 2, ... up to max_deg).
#'
#' @param net Network object.
#' @param max_deg Maximum degree to count (default 20).
#' @return Numeric vector of counts for degrees 0 through max_deg.
#' @noRd
degree_dist <- function(net, max_deg = 20) {
  degs <- sna::degree(net, gmode = "graph")
  tab <- table(factor(degs, levels = 0:max_deg))
  as.vector(tab)
}

#' ESP distribution via ernm (edge-wise shared partners 0, 1, ... k).
#'
#' @param net Network object.
#' @param k_max Maximum ESP count (default 15).
#' @return Numeric vector of ESP counts for 0 through k_max.
#' @noRd
esp_dist <- function(net, k_max = 15) {
  tryCatch({
    as.vector(ernm::calculateStatistics(net ~ esp(0:k_max)))
  }, error = function(e) rep(NA_real_, k_max + 1))
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

#' Goodness-of-fit analysis for hawkesGrowthNet models.
#'
#' Simulates networks from a fitted model and compares observed vs simulated
#' network statistics (degree distribution, ESP, geodesic distances, waiting times).
#'
#' @param fit Fitted hawkesGrowthNet model object (from \code{fit_hawkesGrowthNet} or \code{fit_hawkesGrowthNet_inhom}).
#' @param net_obs Observed network (used for computing observed statistics).
#' @param params_init Initial parameter list used for fitting (needed to reconstruct full parameter structure).
#' @param PMF_mark Mark probability mass function (e.g., \code{PMF_mark_CS}).
#' @param cond_intensity Conditional intensity function (cached closure from fitting).
#' @param formula_RHS Character string RHS of ERNM formula (e.g., "edges + triangles + star(c(2,3))").
#' @param time_window Time window for simulations (default c(0, 0.05)).
#' @param truncation Truncation parameter for mark PMF (default 100).
#' @param mark_decay Mark decay type: "node_entrance" or "activity" (default "activity").
#' @param max_node_time Maximum node time (default 1).
#' @param inhom_bg Optional inhomogeneous background object (for computing mu from integral_bg).
#' @param n_sim Number of simulated networks to generate (default 50).
#' @param cores Number of cores for parallelization (default 7).
#' @param max_deg Maximum degree for degree distribution (default 15).
#' @param k_esp Maximum ESP count for ESP distribution (default 15).
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
                max_deg = 15L, k_esp = 15L, mu_multiplier = 5, verbose = TRUE) {
  
  # Initialize results early (will be populated even if some computations fail)
  GOF_results <- list(degree_obs = NULL, degree_sim = NULL, esp_obs = NULL, esp_sim = NULL,
                      geodist_obs = NULL, geodist_sim = NULL, wait_obs = NULL, wait_sim = NULL,
                      nodemix_obs = NULL, nodemix_sim = NULL, plots = list())
  
  if (is.null(fit)) {
    if (verbose) cat("  No fit available; skipping GOF\n")
    return(GOF_results)
  }
  
  # Wrap entire function body in tryCatch to ensure we always return results
  tryCatch({
    if (verbose) cat("  Simulating", n_sim, "networks from fitted model...\n")
    
    # Reconstruct fitted params: strip vertex_categorical_levels from skeleton
    # (the fitter stripped it before unlist, so fit$par doesn't include it)
    skel <- params_init
    skel$vertex_categorical_levels <- NULL
    pfit <- tryCatch({
      relist(fit$fit$par, skeleton = skel)
    }, error = function(e) {
      if (verbose) cat("  ERROR: Failed to reconstruct parameters:", e$message, "\n")
      if (verbose) cat("  Using params_init as fallback\n")
      params_init
    })
  
    # Restore metadata and fixed params
    pfit$vertex_categorical_levels <- params_init$vertex_categorical_levels
    pfit$K <- params_init$K
    
    # Compute mu from inhomogeneous background if provided
    if (!is.null(inhom_bg)) {
      Tval <- time_window[2] - time_window[1]
      pfit$mu <- inhom_bg$integral_bg / Tval
    } else {
      pfit$mu <- params_init$mu
    }
    
    # Ensure parameters are within valid ranges
    pfit$K <- min(max(pfit$K, 0.001), 0.999)  # K must be in (0,1) for stability
    pfit$mu <- max(pfit$mu, 0.001)  # mu must be positive
    pfit$node_lambda <- max(pfit$node_lambda, 1)  # node_lambda must be >= 1 for stability
    pfit$beta_overall <- max(pfit$beta_overall, 0.001)
    pfit$beta_edges <- max(pfit$beta_edges, 0.001)
    
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
    if (verbose) cat("  Using", cores, "cores for parallel GOF simulations...\n")
    sim_results <- tryCatch({
      parallel::mclapply(seq_len(n_sim), function(i) {
    s <- tryCatch(
      sim_hawkesGrowthNet(
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
        stop_on_full_network = FALSE
      ),
      error = function(e) { 
        return(list(net = NULL, error = paste0("Sim ", i, ": ", e$message))) 
      }
    )
    if (is.null(s$net)) {
      return(list(net = NULL, error = ifelse(is.null(s$error), paste0("Sim ", i, ": unknown error"), s$error)))
    }
        return(list(net = s$net, error = NULL))
      }, mc.cores = cores)
    }, error = function(e) {
      if (verbose) cat("  ERROR: Failed to run simulations:", e$message, "\n")
      list()  # Return empty list if simulations fail completely
    })
    
    if (is.null(sim_results) || length(sim_results) == 0) {
      if (verbose) cat("  No simulation results; returning empty GOF results\n")
      return(GOF_results)
    }
    
    sim_nets <- lapply(sim_results, function(x) x$net)
    sim_nets <- sim_nets[!sapply(sim_nets, is.null)]
    n_success <- length(sim_nets)
    n_fail <- n_sim - n_success
    
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
      # Wrap entire statistics computation in tryCatch to ensure we always return partial results
      tryCatch({
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
      degree_dist(net_obs, max_deg)
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute observed degree distribution:", e$message, "\n")
      NULL
    })
    
    GOF_results$esp_obs <- tryCatch({
      esp_dist(net_obs, k_esp)
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
      waiting_times_between_formations(net_obs, formula_RHS = formula_RHS)
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
    
    # Simulated statistics (parallelized) - each wrapped in tryCatch
    if (verbose) cat("    Computing degree distributions...\n")
    GOF_results$degree_sim <- tryCatch({
      do.call(rbind, parallel::mclapply(sim_nets, function(n) {
        tryCatch({
          degree_dist(n, max_deg)
        }, error = function(e) rep(NA_real_, max_deg + 1))
      }, mc.cores = cores))
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute simulated degree distributions:", e$message, "\n")
      NULL
    })
    
    if (verbose) cat("    Computing ESP distributions...\n")
    GOF_results$esp_sim <- tryCatch({
      do.call(rbind, parallel::mclapply(sim_nets, function(n) {
        tryCatch({
          esp_dist(n, k_esp)
        }, error = function(e) rep(NA_real_, k_esp + 1))
      }, mc.cores = cores))
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute simulated ESP distributions:", e$message, "\n")
      NULL
    })
    
    if (verbose) cat("    Computing geodesic distances...\n")
    GOF_results$geodist_sim <- tryCatch({
      parallel::mclapply(sim_nets, function(n) {
        tryCatch({
          geodist_dist(n)
        }, error = function(e) numeric(0))
      }, mc.cores = cores)
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute simulated geodesic distances:", e$message, "\n")
      NULL
    })
    
    if (verbose) cat("    Computing waiting times...\n")
    GOF_results$wait_sim <- tryCatch({
      parallel::mclapply(sim_nets, function(n) {
        tryCatch({
          waiting_times_between_formations(n, formula_RHS = formula_RHS)
        }, error = function(e) {
          if (verbose && length(sim_nets) <= 5) cat("        Warning: Could not compute waiting times for one sim:", e$message, "\n")
          list()
        })
      }, mc.cores = cores)
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute simulated waiting times:", e$message, "\n")
      NULL
    })
    
    # Simulated nodeMix statistics (parallelized)
    if (verbose) cat("    Computing simulated nodeMix statistics...\n")
    GOF_results$nodemix_sim <- tryCatch({
      # Determine expected length from observed nodeMix
      nodemix_obs_len <- if (!is.null(GOF_results$nodemix_obs)) length(GOF_results$nodemix_obs) else {
        # Try to compute length from observed network
        tryCatch({
          if ("gender" %in% network::list.vertex.attributes(net_obs) || 
              any(grepl("nodeMix|nodeMatch", formula_RHS))) {
            net_obs_clean <- ensure_vertex_attribute(net_obs, "gender", default_value = "unknown")
            length(ernm::calculateStatistics(net_obs_clean ~ nodeMix('gender')))
          } else {
            0
          }
        }, error = function(e) 0)
      }
      
      if (nodemix_obs_len > 0) {
        do.call(rbind, parallel::mclapply(sim_nets, function(n) {
          tryCatch({
            if ("gender" %in% network::list.vertex.attributes(n) || 
                any(grepl("nodeMix|nodeMatch", formula_RHS))) {
              # Ensure gender attribute is properly set before ERNM operations
              n_clean <- ensure_vertex_attribute(n, "gender", default_value = "unknown")
              as.vector(ernm::calculateStatistics(n_clean ~ nodeMix('gender')))
            } else {
              rep(NA_real_, nodemix_obs_len)
            }
          }, error = function(e) {
            rep(NA_real_, nodemix_obs_len)
          })
        }, mc.cores = cores))
      } else {
        NULL
      }
    }, error = function(e) {
      if (verbose) cat("      Warning: Could not compute simulated nodeMix statistics:", e$message, "\n")
      NULL
    })
    
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
      }, error = function(e) {
        # If statistics computation fails completely, return whatever we have
        if (verbose) cat("  ERROR in GOF statistics computation:", e$message, "\n")
        if (verbose) cat("  Returning partial results (some statistics may be NULL)\n")
      })
    } else {
      if (verbose) cat("  No successful simulations; returning empty GOF results\n")
    }
  }, error = function(e) {
    # If anything fails at the top level (parameter reconstruction, simulations, etc.)
    if (verbose) cat("  ERROR in GOF function:", e$message, "\n")
    if (verbose) cat("  Returning empty GOF results\n")
  })
  
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
  create_dist_plot <- function(obs, sim, stat_name, x_label = NULL) {
    if (is.null(obs) || is.null(sim) || nrow(sim) == 0) return(NULL)
    
    # Prepare data
    n_obs <- length(obs)
    n_sim <- nrow(sim)
    
    # Create data frame
    df_obs <- data.frame(
      value = obs,
      x = seq_along(obs),
      type = "Observed"
    )
    
    df_sim <- data.frame(
      value = as.vector(sim),
      x = rep(seq_len(ncol(sim)), each = nrow(sim)),
      type = "Simulated"
    )
    
    df <- rbind(df_obs, df_sim)
    
    # Use x_label if provided, otherwise "Index"
    x_lab <- if (is.null(x_label)) "Index" else x_label
    
    # Create plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(x), y = value, fill = type)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      ggplot2::scale_fill_manual(values = c("Observed" = "#E69F00", "Simulated" = "#56B4E9")) +
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
  
  # Degree distribution plot
  if (!is.null(GOF_results$degree_obs) && !is.null(GOF_results$degree_sim)) {
    plots$degree_plot <- create_dist_plot(
      GOF_results$degree_obs,
      GOF_results$degree_sim,
      "Degree",
      "Degree"
    )
  }
  
  # ESP distribution plot
  if (!is.null(GOF_results$esp_obs) && !is.null(GOF_results$esp_sim)) {
    plots$esp_plot <- create_dist_plot(
      GOF_results$esp_obs,
      GOF_results$esp_sim,
      "ESP (Edge-wise Shared Partners)",
      "ESP Count"
    )
  }
  
  # Geodesic distance plot (histogram style)
  if (!is.null(GOF_results$geodist_obs) && !is.null(GOF_results$geodist_sim) && 
      length(GOF_results$geodist_obs) > 0 && length(GOF_results$geodist_sim) > 0) {
    # Combine observed and simulated geodesic distances
    df_geod <- data.frame(
      distance = c(
        GOF_results$geodist_obs,
        unlist(GOF_results$geodist_sim)
      ),
      type = c(
        rep("Observed", length(GOF_results$geodist_obs)),
        rep("Simulated", sum(sapply(GOF_results$geodist_sim, length)))
      )
    )
    
    plots$geodist_plot <- ggplot2::ggplot(df_geod, ggplot2::aes(x = distance, fill = type)) +
      ggplot2::geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
      ggplot2::scale_fill_manual(values = c("Observed" = "#E69F00", "Simulated" = "#56B4E9")) +
      ggplot2::labs(
        title = "Geodesic Distance Distribution",
        x = "Geodesic Distance",
        y = "Frequency",
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
  
  # Waiting times plot (faceted histograms)
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
        df_wait_list[[stat_name]] <- data.frame(
          waiting_time = c(obs_wait, sim_wait),
          type = c(rep("Observed", length(obs_wait)), rep("Simulated", length(sim_wait))),
          statistic = stat_name
        )
      }
    }
    
    if (length(df_wait_list) > 0) {
      df_wait <- do.call(rbind, df_wait_list)
      
      # Determine number of columns for faceting (max 2)
      n_stats <- length(unique(df_wait$statistic))
      ncol_facet <- min(2, n_stats)
      
      plots$waiting_times_plot <- ggplot2::ggplot(df_wait, ggplot2::aes(x = waiting_time, fill = type)) +
        ggplot2::geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
        ggplot2::facet_wrap(~ statistic, scales = "free", ncol = ncol_facet) +
        ggplot2::scale_fill_manual(values = c("Observed" = "#E69F00", "Simulated" = "#56B4E9")) +
        ggplot2::labs(
          title = "Waiting Times Between Structure Formations",
          x = "Waiting Time",
          y = "Frequency",
          fill = "Type"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "bottom",
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          strip.text = ggplot2::element_text(face = "bold")
        )
    }
  }
  
  plots
}
