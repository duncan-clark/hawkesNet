

#' Create an animation of network growth over time
#'
#' Builds a dynamic network from a list of network snapshots and optional video file.
#' This function is central for visualizing how the network grows with Hawkes process events.
#'
#' @param net_list List of \code{network} objects (one per event time).
#' @param times Numeric vector of event times, same length as \code{net_list}.
#' @param adjust Scaling factor for animation timing (default 10).
#' @param file Optional character path to save the animation as a video file.
#' @return A \code{networkDynamic} object (from package \pkg{networkDynamic}).
#' @details
#' Requires the suggested packages \pkg{networkDynamic} and \pkg{animation} to be installed.
#' Install with \code{install.packages(c("networkDynamic", "animation"))}.
#' The animation displays the network at each event time; if \code{file} is provided,
#' \code{animation::saveVideo} is used to write a video file.
#' @examples
#' \dontrun{
#' nets <- list(network::network(2, directed = FALSE))
#' make_network_growth_animation(nets, times = 0, file = "growth.mp4")
#' }
#' @rdname make_network_growth_animation
#' @export
make_network_growth_animation <- function(net_list,
                                          times,
                                          adjust = 10,
                                          file = NULL){
  if (!requireNamespace("networkDynamic", quietly = TRUE)) {
    stop("Package 'networkDynamic' is required for make_network_growth_animation(). ",
         "Install it with install.packages('networkDynamic').")
  }
  if (!requireNamespace("animation", quietly = TRUE) && !is.null(file)) {
    stop("Package 'animation' is required to save animation to file. ",
         "Install it with install.packages('animation').")
  }

  for(i in seq_along(net_list)){
    set.vertex.attribute(net_list[[i]],'vertex.names',1:(net_list[[i]] %n% 'n'))
  }

  animate <- networkDynamic::networkDynamic(network.list=net_list,
                            onsets = adjust*times,
                            termini = adjust * c(times[seq_len(length(times))[-1]], max(times)),
                            vertex.pid = "vertex.names"
  )
  ndtv::render.animation(animate,
                   render.par = list(tween.frames = 1,
                                     show.time = TRUE,
                                     show.stats = NULL,
                                     extraPlotCmds=NULL,
                                     initial.coords=0),
                   displaylabels = FALSE,
                   displayisolates = TRUE
  )
  if(!is.null(file)){
    animation::saveVideo(animation::ani.replay(),
              video.name=file,
              other.opts="-b 5000k",
              clean=TRUE)

  }
  return(animate)
}

# =============================================================================
# Safe parallel lapply — avoids fork deadlocks on macOS interactive sessions
# =============================================================================
#
# On macOS (Darwin), fork() can deadlock on the *second* call (not just interactive):
# after the first fork+join, ObjC runtime state can make the next mclapply() hang.
# So we avoid fork on Darwin entirely.
#
# This helper routes to:
#   - PSOCK cluster (parallel::parLapply) on macOS and Windows → no fork.
#   - parallel::mclapply on Linux (fast fork).
#
# Users can override with parallel_type = "fork" / "psock" / "auto".
# @noRd
safe_parallel_lapply <- function(X, FUN, mc.cores,
                                 mc.preschedule = FALSE,
                                 parallel_type = "auto") {

  os <- Sys.info()[["sysname"]]

  if (parallel_type == "auto") {
    # PSOCK (no fork) when: Windows (no fork), macOS (fork often deadlocks on 2nd call),
    # or any interactive session (RStudio / GUI are multi-threaded; fork can deadlock).
    # Fork (mclapply) only for non-interactive Linux (e.g. SLURM Rscript).
    use_psock <- (os == "Darwin") || (os == "Windows") || interactive() ||
                 isTRUE(getOption("hawkesNet.force_psock", FALSE))
  } else {
    use_psock <- (parallel_type == "psock")
  }

  if (use_psock) {
    n_workers <- mc.cores
    # Cap PSOCK workers to avoid OOM. Each worker loads packages + a copy of the network.
    # Default: 16 for interactive sessions (laptops), higher for SLURM batch jobs.
    default_cap <- if (!interactive() && nzchar(Sys.getenv("SLURM_JOB_ID"))) 64L else 16L
    max_psock <- getOption("hawkesNet.max_psock_workers", default_cap)
    if (n_workers > max_psock) {
      message("  [parallel] Capping PSOCK workers to ", max_psock, " (avoid OOM; set options(hawkesNet.max_psock_workers = N) to override)")
      n_workers <- as.integer(max_psock)
    } else {
      message("  [parallel] PSOCK cluster (", n_workers, " workers)")
    }
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Load packages that intens_func's captured closures depend on
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(hawkesNet)
        library(ernm)
        library(network)
        library(sna)
      })
    })
    result <- parallel::parLapply(cl, X, FUN)
    return(result)
  }

  # Fork-based: plain mclapply (no pbmcapply progress pipe)
  # Clean up any zombie child processes from previous mclapply calls.
  # Without this, the second mclapply in the same session can hang because
  # stale pipes/signal handlers from the first call interfere.
  tryCatch({
    children_fn  <- get("children",  envir = asNamespace("parallel"))
    mccollect_fn <- get("mccollect", envir = asNamespace("parallel"))
    active_children <- children_fn()
    if (length(active_children) > 0L) {
      message("  [parallel] Cleaning up ", length(active_children), " zombie child processes...")
      mccollect_fn(active_children, wait = FALSE, timeout = 2)
      # Force kill if they still persist
      still_active <- children_fn()
      if (length(still_active) > 0L) {
        message("  [parallel] WARNING: ", length(still_active), " children still active, sending SIGKILL")
        tools::pskill(still_active, tools::SIGKILL)
        mccollect_fn(still_active, wait = FALSE)
      }
    }
  }, error = function(e) {
    message("  [parallel] Error during zombie cleanup: ", e$message)
  })
  
  # Log memory state before forking
  if (!interactive()) {
    mem <- gc()
    message(sprintf("  [parallel] Parent memory: %.1f Mb (Vcells used)", mem[2, 2]))
  }
  
  # CRITICAL: Disable multi-threading in BLAS/OpenMP before forking.
  # Many BLAS libraries (OpenBLAS, MKL) are not fork-safe and will deadlock
  # on the second fork if a thread pool was initialized in the first.
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas <- RhpcBLASctl::blas_get_num_procs()
    old_omp  <- RhpcBLASctl::omp_get_max_threads()
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    on.exit({
      RhpcBLASctl::blas_set_num_threads(old_blas)
      RhpcBLASctl::omp_set_num_threads(old_omp)
    }, add = TRUE)
  } else {
    # Fallback: set environment variables (only works if set before library load,
    # but some libraries check them dynamically or we can at least try).
    Sys.setenv(OMP_NUM_THREADS = "1")
    Sys.setenv(MKL_NUM_THREADS = "1")
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  }
  
  gc()  # reclaim memory before forking so children inherit a lean process
  message("  [parallel] mclapply (", mc.cores, " cores, fork",
          if (mc.preschedule) ", preschedule" else "", ")")
  
  # Wrap FUN to log task start/end in children
  # We use indices to allow logging task numbers regardless of X's content
  FUN_wrapped_idx <- function(idx) {
    i_val <- X[[idx]]
    # Only log for a subset of tasks to avoid flooding
    # For smaller lists (like GOF simulations), log more frequently
    should_log <- (length(X) <= 100) || (idx == 1L || idx == length(X) || (idx %% 50 == 0))
    if (should_log) {
      cat(sprintf("  [parallel-child] Task %d/%d starting (pid %d) at %s\n", 
                  idx, length(X), Sys.getpid(), format(Sys.time(), "%H:%M:%S")), file = stderr())
    }
    
    out <- tryCatch({
      res_fun <- FUN(i_val)
      if (should_log) {
        cat(sprintf("  [parallel-child] Task %d/%d complete (pid %d) at %s\n", 
                    idx, length(X), Sys.getpid(), format(Sys.time(), "%H:%M:%S")), file = stderr())
      }
      res_fun
    }, error = function(e) {
      cat(sprintf("  [parallel-child] Task %d FAILED (pid %d): %s\n", idx, Sys.getpid(), e$message), file = stderr())
      # Return a try-error so the parent can detect it
      structure(list(message = e$message, call = e$call), class = "try-error")
    })
    
    out
  }

  res <- parallel::mclapply(seq_along(X), FUN_wrapped_idx, mc.cores = mc.cores,
                            mc.preschedule = mc.preschedule)
  
  # Check for errors in results (mclapply returns try-error or NULL on some failures)
  if (!is.list(res)) {
    # If res is not a list (e.g. character vector of errors), wrap it
    message("  [parallel] CRITICAL: mclapply did not return a list. This usually indicates a major fork failure.")
    res <- lapply(res, function(x) structure(list(message = as.character(x)), class = "try-error"))
  }
  
  errors <- vapply(res, function(x) inherits(x, "try-error"), logical(1))
  if (any(errors)) {
    message("  [parallel] WARNING: ", sum(errors), " tasks failed in mclapply")
  }
  nulls <- vapply(res, is.null, logical(1))
  if (all(nulls) && length(res) > 0) {
    message("  [parallel] CRITICAL: All tasks returned NULL. This often indicates a fork crash (OOM or deadlock).")
  }
  
  return(res)
}

# Helper function for edge hash maps
#' Check whether an edge exists in a hash-backed edge set
#'
#' @param i Vertex index (head).
#' @param j Vertex index (tail).
#' @param edge_hash A \code{hash} object (from package \pkg{hash}) with keys of the form \code{"i-j"}.
#' @return Logical: \code{TRUE} if the edge (i, j) is in the hash.
#' @rdname has_edge
#' @export
has_edge <- function(i, j, edge_hash) {
  existing_keys <- keys(edge_hash)
  keys_to_check <- paste(i, j, sep = "-")
  keys_to_check %in% existing_keys
}


#' Convert an event list (edge additions over time) to a network
#'
#' Builds a \code{network} object by applying edge events in order; optionally starts from an existing network.
#'
#' @param events_list List with components \code{i} (head vertex indices), \code{j} (tail), and optionally \code{t} (times).
#' @param net Optional starting \code{network} object; if \code{NULL}, one is created.
#' @param directed Logical; whether the network is directed (default \code{FALSE}).
#' @return A \code{network} object with vertices and edges from \code{events_list}.
#' @seealso \code{\link[network]{network}}, \code{\link[network]{add.vertices}}, \code{\link[network]{set.edge.attribute}}
#' @rdname events_to_net
#' @export
events_to_net <- function(events_list,
                          net = NULL,
                          directed = FALSE){

  for(k in seq_along(events_list$i)){
    if(k==1 & is.null(net)){
      net <- network::network(matrix(c(events_list$i[k],events_list$j[k]),nrow = 1),directed = directed)
      network::set.vertex.attribute(net,"time",events_list$t[k],v=events_list$i[k])
      network::set.vertex.attribute(net,"time",events_list$t[k],v=events_list$j[k])
    }else{
      N <- net %n% 'n'
      over_i <- events_list$i[k] - N
      over_j <- events_list$j[k] - N
      if(over_i > 0){
        net <- network::add.vertices(net,over_i)
        network::set.vertex.attribute(net,"time",events_list$t[k],v=events_list$i[k])
        # since network is now bigger amend the over j
        over_j <- over_j - 1
      }
      if(over_j > 0){
        net <- network::add.vertices(net,over_j)
        network::set.vertex.attribute(net,"time",events_list$t[k],v=events_list$j[k])
      }
      add.edge(net, events_list$i[k], events_list$j[k])
    }
    e <- get.dyads.eids(net,
                        events_list$i[k],
                        events_list$j[k])
    network::set.edge.attribute(net,
                                "time",
                                events_list$t[k],
                                e = e[[1]]
    )
  }
  set.network.attribute(net,'n',max(c(events_list$i,events_list$j)))
  return(net)
}

#' Extract the filtration (subnetwork) at or before a given time
#'
#' Removes edges and vertices with time greater than \code{t}; optionally excludes those exactly at \code{t}.
#'
#' @param net A \code{network} object with vertex and edge attribute \code{time}.
#' @param t Numeric; cutoff time.
#' @param equals If \code{FALSE}, also remove edges/vertices with time exactly equal to \code{t} (default \code{FALSE}).
#' @return The filtered \code{network} object.
#' @rdname filtration_to_net
#' @export
filtration_to_net <- function(net,
                              t,
                              equals = FALSE){
  # CRITICAL: network objects can be modified in-place by delete.edges/vertices.
  # We must work on a copy to avoid corrupting the original network, especially
  # when this is called inside parallel workers or loops.
  net <- network::network.copy(net)
  
  # make sure to leave one less edge or vertex that if equals
  delete.edges(net, which(get.edge.attribute(net,"time")>t))
  delete.vertices(net,which(get.vertex.attribute(net,"time")>t))
  if(!equals){
    e_times <- get.edge.attribute(net,"time")
    n_times <- get.vertex.attribute(net,"time")
    delete.edges(net,which(e_times == t))
    delete.vertices(net,which(n_times == t))
  }
  # no need for vertex names
  delete.vertex.attribute(net,'vertex.names')
  return(net)
}

filtration_to_net_both <- function(net,t){
  e_times <- get.edge.attribute(net,"time")
  n_times <- get.vertex.attribute(net,"time")
  
  times <- c(e_times,n_times)
  t_to_delete <- max(times[times <= t])
  net_less <- net
  
  del_e_less <- which(e_times>t)
  del_n_less <- which(n_times>t)
  
  del_e_eq <- setdiff(which(e_times>=t_to_delete),del_e_less) 
  del_n_eq <- setdiff(which(n_times>=t_to_delete),del_n_less)
  
  
  # make sure to leave one less edge or vertex that if equals
  delete.edges(net_less, which(e_times>t))
  delete.vertices(net_less,which(n_times>t))
  
  net_eq <- net_less
  delete.edges(net_eq, del_e_eq)
  delete.vertices(net_eq, del_n_eq)
  # no need for vertex names
  delete.vertex.attribute(net,'vertex.names')
  return(list(eq = net_eq,
              less = net_less))
}

#' Get event times from a network
#'
#' Extracts vertex and edge \code{time} attributes and returns a sorted unique vector of all event times.
#'
#' @param net A \code{network} object with \code{time} attribute on vertices and edges.
#' @param time_name Character; name of the time attribute (default \code{"time"}).
#' @return List with \code{node_times}, \code{edge_times}, and \code{times} (sorted unique).
#' @rdname get_times
#' @export
get_times <- function(net, time_name = 'time'){
  node_times <- network::get.vertex.attribute(net,time_name)
  edge_times <- network::get.edge.attribute(net,time_name)
  return(list(node_times = node_times,
              edge_times = edge_times,
              times = sort(unique(c(node_times,edge_times)))
  ))
}

# function to plot pp on line:
pp_line_plot <- function(t,title=NULL){
  plot(c(min(t),max(t)), c(-1, 1), type = "n", yaxt = "n",
     xlab = "Value", ylab = "", main = paste0("Vector on a Number Line: ",title))
abline(h = 0, col = "gray", lwd = 2)
points(t, rep(0, length(t)), pch = 19, col = "blue", cex = 1.5)
}

plot_kde_intensity <- function(event_times,
                               bw_adjust = 0.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_kde_intensity(). ",
         "Install it with install.packages('ggplot2').")
  }
  # event_times: numeric vector (0 to 1)
  df <- data.frame(time = event_times)
  
  ggplot2::ggplot(df, ggplot2::aes(x = time)) +
    # 1. The KDE Line (Raw Intensity)
    ggplot2::stat_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      geom = "line",
      color = "#2c3e50",
      linewidth = 1.2,
      adjust = bw_adjust
    ) +
    # 2. Add rug marks to see the actual event locations
    ggplot2::geom_rug(alpha = 0.4, color = "firebrick") +
    # 3. Formatting
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title = "Estimated Event Intensity over Time",
      subtitle = paste0("KDE Line Plot (Bandwidth Adjust = ", bw_adjust, ")"),
      x = "Normalized Time (0 = Oldest, 1 = Newest)",
      y = "Intensity (Event Density)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold")
    )
}

# --- Execution ---
# times <- network_data %v% "time_scaled"
# plot_kde_intensity(times, bw_adjust = 0.3)

get_latest_times <- function(nw) {
  el <- as.matrix.network.edgelist(nw, names = FALSE)
  times <- get.edge.attribute(nw, "time")
  
  n <- network.size(nw)
  best <- rep(-Inf, n)
  
  # Treat NA times as -Inf so they don't win
  t2 <- times
  t2[is.na(t2)] <- -Inf
  
  # One vectorized pass
  best[el[, 1]] <- pmax(best[el[, 1]], t2)
  best[el[, 2]] <- pmax(best[el[, 2]], t2)
  
  # Convert untouched nodes to NA
  best[is.infinite(best)] <- NA_real_
  
  # Fallback to node times
  node_times <- get_times(nw)$node_times
  idx <- is.na(best)
  best[idx] <- node_times[idx]
  
  best
}

#' Normalize node + edge "time" attributes to `[0, 1]`
#'
#' Uses get_times() to compute the global min/max across BOTH node and edge times,
#' then rescales vertex and edge "time" attributes into `[0, 1]`.
#'
#' @param net A `network` object.
#' @param attr Name of the time attribute (default "time").
#' @param keep_na If TRUE, leave NA times as NA (default TRUE). If FALSE, error on NA.
#' @param constant_value Value to assign when all non-NA times are identical (default 0).
#' @return The network with normalized times.
#' @export
normalize_times_01 <- function(net, attr = "time", keep_na = TRUE, constant_value = 0) {
  # Pull times using your existing helper
  times_obj <- get_times(net)
  
  # Global range across node+edge times (unique & sorted already)
  all_times <- times_obj$times
  
  if (length(all_times) == 0L) return(net)
  
  if (!keep_na && anyNA(all_times)) {
    stop("NA times found; set keep_na = TRUE to keep them as NA.")
  }
  
  rng <- range(all_times, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(net)
  
  # Helper
  norm01 <- function(x) (x - rng[1]) / (rng[2] - rng[1])
  
  # Vertex times
  v_times <- network::get.vertex.attribute(net, attr)
  if (!is.null(v_times)) {
    v_times <- as.numeric(v_times)
    if (rng[2] == rng[1]) {
      v_times[!is.na(v_times)] <- constant_value
    } else {
      idx <- !is.na(v_times)
      v_times[idx] <- norm01(v_times[idx])
    }
    network::set.vertex.attribute(net, attr, v_times)
  }
  
  # Edge times
  e_times <- network::get.edge.attribute(net, attr)
  if (!is.null(e_times)) {
    e_times <- as.numeric(e_times)
    if (rng[2] == rng[1]) {
      e_times[!is.na(e_times)] <- constant_value
    } else {
      idx <- !is.na(e_times)
      e_times[idx] <- norm01(e_times[idx])
    }
    network::set.edge.attribute(net, attr, e_times)
  }
  
  net
}


#' Check whether point process parameters are in valid regions
#'
#' Returns \code{FALSE} if any of \code{mu}, \code{beta_overall}, \code{K},
#' \code{beta_edges}, \code{node_lambda} are present but not strictly positive
#' and finite, or if \code{vertex_categorical} probabilities are not valid
#' (non-negative, finite, and sum < 1 so the reference level gets a positive
#' probability). Uses the n-1 parametrization: user supplies n-1 probabilities
#' and the last level's probability is \code{1 - sum(p)}.
#'
#' @param params List of parameters (e.g. from \code{relist} or passed to \code{sim_hawkesNet}).
#' @param eps Scalar rate/scale params must be \code{> eps} (default \code{1e-10}).
#' @return \code{TRUE} if all present parameters are valid, \code{FALSE} otherwise.
#' @noRd
point_process_params_valid <- function(params, eps = 1e-10) {
  if (is.null(params) || length(params) == 0) return(TRUE)
  eps <- max(eps, .Machine$double.eps)
  scalar_ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > eps
  if (!is.null(params$mu) && !scalar_ok(params$mu)) return(FALSE)
  if (!is.null(params$beta_overall) && !scalar_ok(params$beta_overall)) return(FALSE)
  if (!is.null(params$K) && !scalar_ok(params$K)) return(FALSE)
  if (!is.null(params$beta_edges) && !scalar_ok(params$beta_edges)) return(FALSE)
  if (!is.null(params$node_lambda) && !scalar_ok(params$node_lambda)) return(FALSE)
  if (!is.null(params$m) && !scalar_ok(params$m)) return(FALSE)
  if (!is.null(params$vertex_categorical) && is.list(params$vertex_categorical)) {
    for (attr_name in names(params$vertex_categorical)) {
      p <- params$vertex_categorical[[attr_name]]
      if (!is.numeric(p) || length(p) == 0L) return(FALSE)
      if (any(!is.finite(p)) || any(p < 0)) return(FALSE)
      s <- sum(p)
      # n-1 parametrization: sum must be < 1 (reference level gets 1 - sum)
      if (!is.finite(s) || s <= 0 || s >= 1) return(FALSE)
    }
  }
  TRUE
}

#' Reconstruct vertex_categorical parameters with correct names after relist
#'
#' When parameters are flattened and then relisted, vertex_categorical names may be lost.
#' This function ensures names are correctly restored from vertex_categorical_levels.
#'
#' @param params Parameter list (after relist).
#' @param vertex_categorical_levels Level names for each attribute.
#' @return Parameter list with correctly named vertex_categorical.
#' @noRd
reconstruct_vertex_categorical_names <- function(params, vertex_categorical_levels) {
  if (is.null(params$vertex_categorical) || !is.list(params$vertex_categorical)) {
    return(params)
  }
  if (is.null(vertex_categorical_levels) || !is.list(vertex_categorical_levels)) {
    return(params)
  }
  
  for (attr_name in names(params$vertex_categorical)) {
    p <- params$vertex_categorical[[attr_name]]
    levs <- vertex_categorical_levels[[attr_name]]
    
    if (!is.null(levs) && length(levs) > 1L && is.numeric(p)) {
      n_levs <- length(levs) - 1L
      if (length(p) == n_levs) {
        # Restore names from levels (n-1 parametrization: exclude reference level)
        names(p) <- levs[1:n_levs]
        params$vertex_categorical[[attr_name]] <- p
      }
    }
  }
  
  params
}

repair_vertex_categorical_params <- function(params, eps = 1e-6) {
  if (is.null(params$vertex_categorical) || !is.list(params$vertex_categorical)) {
    return(params)
  }
  eps <- max(eps, .Machine$double.eps)
  
  for (attr_name in names(params$vertex_categorical)) {
    p <- params$vertex_categorical[[attr_name]]
    
    # If p has no names but we have levels, restore names first
    if (is.numeric(p) && length(p) > 0L && is.null(names(p))) {
      levs <- params$vertex_categorical_levels[[attr_name]]
      if (!is.null(levs) && length(levs) > 1L) {
        n_levs <- length(levs) - 1L
        if (length(p) == n_levs) {
          names(p) <- levs[1:n_levs]
        }
      }
    }
    
    if (!is.numeric(p) || length(p) == 0L) {
      # Invalid: set to default (equal probabilities)
      levs <- params$vertex_categorical_levels[[attr_name]]
      if (!is.null(levs) && length(levs) > 1L) {
        n_levs <- length(levs) - 1L
        default_val <- (1 - eps * n_levs) / n_levs
        p <- setNames(rep(default_val, n_levs), levs[1:n_levs])
      } else {
        next
      }
    }
    
    # Only repair if parameters are INVALID (negative, non-finite, or sum >= 1)
    # Do NOT clamp valid small positive values - they might be legitimate estimates
    s <- sum(p)
    needs_repair <- FALSE
    
    # Check for invalid values
    if (any(!is.finite(p)) || any(p < 0) || !is.finite(s) || s <= 0 || s >= 1) {
      needs_repair <- TRUE
    }
    
    if (needs_repair) {
      # Repair invalid parameters
      # First, handle non-finite and negative values
      p[!is.finite(p)] <- eps
      p[p < 0] <- eps
      
      # Recompute sum after fixing non-finite/negative
      s <- sum(p)
      
      if (!is.finite(s) || s <= 0) {
        # Invalid sum: set to default (equal probabilities)
        levs <- params$vertex_categorical_levels[[attr_name]]
        if (!is.null(levs) && length(levs) > 1L) {
          n_levs <- length(levs) - 1L
          default_val <- (1 - eps * n_levs) / n_levs
          p <- setNames(rep(default_val, n_levs), levs[1:n_levs])
        }
      } else if (s >= 1) {
        # Scale down proportionally to ensure sum < 1 (n-1 parametrization)
        # Reference level gets probability 1 - sum(p), so we need sum(p) <= 1 - eps
        max_sum <- 1 - eps
        if (max_sum > eps && s > 0) {
          p <- p * (max_sum / s)
          # After scaling, ensure individual values are still >= eps
          p <- pmax(p, eps)
          # Re-scale if needed to maintain sum <= max_sum
          s_new <- sum(p)
          if (s_new > max_sum) {
            p <- p * (max_sum / s_new)
          }
        } else {
          # If even max_sum is too small, use equal probabilities
          levs <- params$vertex_categorical_levels[[attr_name]]
          if (!is.null(levs) && length(levs) > 1L) {
            n_levs <- length(levs) - 1L
            default_val <- max_sum / n_levs
            p <- setNames(rep(default_val, n_levs), names(p)[1:n_levs])
          }
        }
      }
      # After repair, ensure values are in valid range
      p <- pmax(p, eps)
      p <- pmin(p, 1 - eps)
    }
    # If parameters are valid (finite, non-negative, sum < 1), leave them as-is
    # This preserves legitimate small values from optimization
    
    params$vertex_categorical[[attr_name]] <- p
  }
  
  params
}

#' Validate point process parameters and stop if invalid
#'
#' Checks that \code{mu}, \code{beta_overall}, \code{K}, \code{beta_edges},
#' \code{node_lambda} are strictly positive and finite when present, and that
#' \code{vertex_categorical} probabilities (n-1 parametrization) are non-negative,
#' finite, and sum to strictly less than 1 (so the reference level gets positive
#' probability).  If any check fails, \code{stop()} is called with a message.
#'
#' @param params List of parameters (e.g. passed to \code{sim_hawkesNet}).
#' @param eps Scalar params must be \code{> eps} (default \code{1e-10}).
#' @return \code{invisible(params)} if valid.
#' @export
validate_point_process_params <- function(params, eps = 1e-10) {
  if (is.null(params) || length(params) == 0) return(invisible(params))
  eps <- max(eps, .Machine$double.eps)
  msg <- character(0L)
  scalar_check <- function(name, x) {
    if (is.null(x)) return(invisible(NULL))
    if (!is.numeric(x) || length(x) != 1L) {
      msg <<- c(msg, paste0(name, " must be a numeric scalar"))
      return(invisible(NULL))
    }
    if (!is.finite(x)) {
      msg <<- c(msg, paste0(name, " must be finite (got ", x, ")"))
      return(invisible(NULL))
    }
    if (x <= eps) {
      msg <<- c(msg, paste0(name, " must be > ", eps, " (got ", x, ")"))
      return(invisible(NULL))
    }
    invisible(NULL)
  }
  scalar_check("mu", params$mu)
  scalar_check("beta_overall", params$beta_overall)
  scalar_check("K", params$K)
  scalar_check("beta_edges", params$beta_edges)
  scalar_check("node_lambda", params$node_lambda)
  scalar_check("m", params$m)
  if (!is.null(params$vertex_categorical) && is.list(params$vertex_categorical)) {
    for (attr_name in names(params$vertex_categorical)) {
      p <- params$vertex_categorical[[attr_name]]
      if (!is.numeric(p) || length(p) == 0L) {
        msg <- c(msg, paste0("vertex_categorical$", attr_name, " must be a non-empty numeric vector"))
        next
      }
      if (any(!is.finite(p))) {
        msg <- c(msg, paste0("vertex_categorical$", attr_name, " must have finite values"))
        next
      }
      if (any(p < 0)) {
        msg <- c(msg, paste0("vertex_categorical$", attr_name, " must have non-negative probabilities"))
        next
      }
      s <- sum(p)
      # n-1 parametrization: sum of n-1 probs must be < 1 (reference level gets 1 - sum)
      if (!is.finite(s) || s <= 0 || s >= 1) {
        msg <- c(msg, paste0("vertex_categorical$", attr_name, " must have sum in (0, 1) for n-1 parametrization (got ", s, ")"))
      }
    }
  }
  if (length(msg) > 0L) {
    stop("Invalid point process parameters: ", paste(msg, collapse = "; "))
  }
  invisible(params)
}


# intensity_ggplot:
pp_intensity_ggplot <- function(times,
                                line_multiplier = 1,
                                base_multiplier = 1,
                                tlim = NULL,
                                dt = NULL,
                                mu = 0,
                                K = 1,
                                beta = 1,
                                spikes = TRUE,
                                smooth = TRUE,
                                title = "Point process intensity-style plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for pp_intensity_ggplot(). ",
         "Install it with install.packages('ggplot2').")
  }
  times <- sort(as.numeric(times))
  times <- times[is.finite(times)]
  if (!length(times)) stop("times is empty after removing non-finite values.")
  
  if (is.null(tlim)) tlim <- range(times)
  t0 <- tlim[1]; t1 <- tlim[2]
  
  if (is.null(dt)) dt <- (t1 - t0) / 2000
  dt <- max(dt, .Machine$double.eps)
  
  df_events <- data.frame(time = times)
  
  plot_layers <- list()
  
  # ---- Smooth Hawkes-style intensity ----
  if (smooth) {
    grid <- seq(t0, t1, by = dt)
    
    # Efficient recursion
    decay <- exp(-beta * dt)
    s <- numeric(length(grid))
    
    idx <- findInterval(times, grid, rightmost.closed = TRUE)
    counts <- tabulate(pmax(1L, idx), nbins = length(grid))
    
    for (k in 2:length(grid)) {
      s[k] <- decay * s[k - 1] + counts[k - 1]
    }
    
    lambda <- mu + K * s
    df_lambda <- data.frame(time = grid, intensity = lambda)
    
    plot_layers <- c(
      plot_layers,
      list(
        ggplot2::geom_line(
          data = df_lambda,
          ggplot2::aes(time, intensity),
          linewidth = 1*line_multiplier
        )
      )
    )
    
    ymax <- max(lambda)
  } else {
    ymax <- 1
  }
  
  # ---- Event spikes ----
  if (spikes) {
    plot_layers <- c(
      plot_layers,
      list(
        ggplot2::geom_linerange(
          data = df_events,
          ggplot2::aes(x = time, ymin = 0, ymax = ymax),
          alpha = 0.6
        )
      )
    )
  }
  
  ggplot2::ggplot() +
    plot_layers +
    ggplot2::coord_cartesian(xlim = c(t0, t1), ylim = c(0, ymax)) +
    ggplot2::labs(
      x = "time",
      y = "intensity",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 1*base_multiplier)
}