# Helper: build lower/upper bounds for L-BFGS-B from flattened parameter names.
# Parameters whose base name is in `positive_params` get lower = eps; all others get -Inf.
# Upper bounds: vertex_categorical params get upper = 1 - eps; all others get Inf.
build_optim_bounds <- function(par_names, eps = 1e-6) {
  positive_params <- c("mu", "beta_overall", "beta_edges", "K", "node_lambda", "m")
  lower <- rep(-Inf, length(par_names))
  upper <- rep(Inf, length(par_names))
  for (i in seq_along(par_names)) {
    # Strip trailing digits for vector params (e.g. "CS_params1" -> "CS_params")
    base_name <- sub("[0-9]+$", "", par_names[i])
    # Also handle nested names like "vertex_categorical.gender1"
    top_name <- strsplit(par_names[i], "\\.")[[1L]][1L]
    if (base_name %in% positive_params || top_name %in% positive_params) {
      lower[i] <- eps
    }
    if (top_name == "vertex_categorical") {
      lower[i] <- eps
      upper[i] <- 1 - eps
    }
  }
  list(lower = lower, upper = upper)
}

# Helper: rename CS_params1, CS_params2, ... in a fit table to actual ERNM stat names
# e.g. "CS_params1" -> "edges", "CS_params2" -> "triangles", etc.
rename_CS_params_in_table <- function(fit_table, mark_filtration, dot_args) {
  formula_RHS <- dot_args$formula_RHS
  if (is.null(formula_RHS) || is.null(mark_filtration)) return(fit_table)
  exp_cs <- tryCatch(
    expected_params_PMF_mark_CS(mark_filtration, formula_RHS),
    error = function(e) NULL
  )
  if (is.null(exp_cs) || is.null(exp_cs$CS_params_names)) return(fit_table)
  stat_names <- exp_cs$CS_params_names
  n_cs <- length(stat_names)
  # Find all CS_params in the fit table (may not be sequential if other params are interspersed)
  cs_indices <- grep("^CS_params[0-9]+$", fit_table$parameter)
  if (length(cs_indices) != n_cs) {
    # Mismatch: try to match by extracting numbers
    warning("rename_CS_params_in_table: Found ", length(cs_indices), 
            " CS_params in fit table but expected ", n_cs, 
            " statistics from formula. Some parameters may not be renamed.")
    
  }
  # Extract parameter numbers and rename
  for (idx in cs_indices) {
    old_name <- fit_table$parameter[idx]
    # Extract number: CS_params1 -> 1, CS_params10 -> 10
    param_num <- as.integer(sub("^CS_params", "", old_name))
    if (!is.na(param_num) && param_num >= 1L && param_num <= n_cs) {
      fit_table$parameter[idx] <- stat_names[param_num]
    }
  }
  fit_table
}

#' Conditional intensity for Hawkes network growth model
#'
#' Evaluates the conditional intensity at time \code{t} given the mark (network) and past events.
#'
#' @param new_net Current network (mark) state.
#' @param t Current time.
#' @param mark_filtration Observed network filtration (history).
#' @param PMF_mark Mark PMF function (e.g. \code{PMF_mark_BA} or \code{PMF_mark_CS}).
#' @param params List of parameters (\code{mu}, \code{beta_overall}, \code{K}, etc.).
#' @param new_edge_hash Optional hash of existing edges for fast lookup.
#' @param times Optional precomputed event times; if \code{NULL}, taken from \code{mark_filtration}.
#' @param ... Arguments passed to \code{PMF_mark}.
#' @return List with \code{result} (intensity value), \code{func} (function to evaluate intensity at new params), and optional debug fields.
#' @rdname cond_intensity
#' @export
cond_intensity <- function(new_net,
                           t,
                           mark_filtration = NULL,
                           PMF_mark,
                           params,
                           new_edge_hash = NULL,
                           times = NULL,
                           ...) {
  tmp <- PMF_mark(time = t,
                  params = params,
                  mark_filtration = mark_filtration,
                  mark = new_net,
                  generate_mark = FALSE,
                  new_edge_hash = new_edge_hash,
                  ...)
  if (is.null(times)) times <- get_times(mark_filtration)
  tt <- times$times
  tt <- tt[tt < t]
  diffs <- t - tt
  
  # Pull ONLY what we need out of tmp (do NOT keep tmp in the closure)
  log_mark_density0 <- tmp$log_mark_density
  log_density_func  <- tmp$log_density_func
  
  # Initial value at current params (optional but you do it)
  decays0 <- exp(-params$beta_overall * diffs)
  log_result0 <- log_mark_density0 + log(params$mu + params$K * sum(decays0))
  result0 <- exp(log_result0)
  
  # ========================
  # 1. Create the 'tiny' environment first
  # Using emptyenv() as parent is safest for "tiny", but baseenv() is needed 
  # for functions like exp() and log() to work inside the closure.
  e_tiny <- new.env(parent = baseenv()) 
  
  # 2. Manually assign ONLY what you need
  e_tiny$diffs_local <- diffs
  e_tiny$ldf_local   <- log_density_func
  
  # 3. Define the function
  # Note: We define it normally, then swap the environment.
  func_template <- function(params) {
    decays <- exp(-params$beta_overall * diffs_local)
    
    # Note: logic checks out, baseenv contains exp/log/sum
    exp(ldf_local(params) + log(params$mu + params$K * sum(decays)))
  }
  
  # 4. Attach the tiny environment
  func <- func_template
  environment(func) <- e_tiny
  # ==============================
  
  # Optional: sanity check what it captured (comment out in production)
  # stopifnot(identical(ls(environment(func)), c("diffs_local","ldf_local")))
  
  list(
    result = result0,
    func   = func,
    lambda = params$mu,
    kernel_sum = params$K * sum(decays0),
    decays = decays0,
    diffs  = diffs
    # (keep other outputs if you truly need them; they increase memory)
  )
}

#' Simulate a Hawkes-driven network growth process
#'
#' Uses thinning to simulate event times and marks (network edges) from the Hawkes growth model.
#' Supports both homogeneous (constant mu) and inhomogeneous (time-varying mu) background rates.
#'
#' @param params List of parameters (\code{mu}, \code{beta_overall}, \code{K}, \code{beta_edges}, and any mark-specific).
#'   For inhomogeneous background, \code{mu} is ignored and \code{mu_at_t} is used instead.
#' @param time_window Numeric vector \code{c(t0, t1)}.
#' @param PMF_mark Mark PMF function (e.g. \code{PMF_mark_BA} or \code{PMF_mark_CS}).
#' @param cond_intensity Conditional intensity function (e.g. \code{cond_intensity} or \code{cond_intensity_inhom}).
#'   If \code{inhom_bg} is provided, \code{cond_intensity_inhom} will be used automatically.
#' @param hashed_edges If \code{TRUE}, use a hash for edge lookup (default \code{FALSE}).
#' @param verbose Print progress (default \code{FALSE}).
#' @param mu_multiplier Multiplier for thinning upper bound (default 10).
#' @param joint_accept Logical (default \code{FALSE}); joint acceptance for mark and time.
#' @param n_mark_sample Optional number of mark samples per proposal.
#' @param stop_on_full_network If \code{TRUE} (default), stop when there are no candidate edges (full network); if \code{FALSE}, issue a warning and continue with no new edges added for that event.
#' @param inhom_bg Optional inhomogeneous background object from \code{prepare_inhomogeneous_background}.
#'   If provided, uses \code{cond_intensity_inhom} with time-varying background rate.
#' @param ... Passed to \code{PMF_mark} or \code{cond_intensity} (e.g. \code{truncation}, \code{formula_RHS}).
#' @return List with \code{events}, \code{net}, \code{accept_probs}.
#' @seealso \code{\link[network]{as.edgelist}}, \code{\link[hash]{hash}}, \code{\link{cond_intensity_inhom}}, \code{\link{prepare_inhomogeneous_background}}
#' @rdname sim_hawkesNet
#' @export
sim_hawkesNet <- function(params,
                                time_window,
                                PMF_mark, # function that both generates new mark and calculates the density of existing mark
                                cond_intensity, # function to calcualte condiational_intensity, takes in a kernel_func
                                hashed_edges = FALSE,
                                verbose = FALSE,
                                mu_multiplier = 10,
                                joint_accept = FALSE,
                                n_mark_sample = NULL,
                                stop_on_full_network = TRUE, # if TRUE, stop when no candidate edges (full network); if FALSE, warn and continue with no new edges
                                inhom_bg = NULL, # optional inhomogeneous background object
                                ... # to be past to PMF_mark


){
  t1 <- proc.time()
  validate_point_process_params(params)
  
  # Determine if using inhomogeneous background
  use_inhom <- !is.null(inhom_bg) && !is.null(inhom_bg$mu_fit) && !is.null(inhom_bg$mu_fit$mu_fun)
  
  if (use_inhom) {
    # Inhomogeneous: get mu_fun and compute max mu for thinning bound
    mu_fun <- inhom_bg$mu_fit$mu_fun
    # Evaluate mu_fun on a grid to find maximum for thinning bound
    t_grid <- seq(time_window[1], time_window[2], length.out = 1000)
    mu_grid <- pmax(mu_fun(t_grid), 1e-12)
    mu_max <- max(mu_grid, na.rm = TRUE)
    if (!is.finite(mu_max) || mu_max <= 0) {
      mu_max <- params$mu  # Fallback to params$mu if mu_fun fails
    }
    lambda <- mu_multiplier * mu_max
    if (verbose) {
      cat("Using inhomogeneous background: mu_max =", mu_max, "for thinning bound\n")
    }
  } else {
    # Homogeneous: use constant mu
    mu <- params$mu
    lambda <- mu_multiplier * mu
  }
  
  theta <- params$theta
  beta <- params$beta
  K <- params$K

  # propose points to be thinned:
  n_bg <- rpois(1, lambda * (time_window[2] - time_window[1]))

  # Pre-allocate buffers to avoid O(n^2) vector append in the while loop
  event_times_buf <- numeric(n_bg)
  mark_density_buf <- numeric(n_bg)
  accept_probs_buf <- numeric(n_bg)
  n_accepted <- 0L
  n_proposed <- 0L
  n_mark_dens <- 0L
  event_queue <- data.table(time = sort(runif(n_bg, min=0, max=time_window[2])))
  # maintain order so no need to sort
  setkey(event_queue, time)
  # Initialize the list to store new events
  new_events_list <- list()
  list_index <- 1
  current_net <- network::network(matrix(1),directed = FALSE)
  delete.vertices(current_net,1)

  # --- O(1) kernel recurrence state ---
  # Instead of recomputing sum(exp(-beta*(t - t_i))) = O(N) each event,

  # maintain: kernel_R such that kernel_sum = K * kernel_R.
  # Update: kernel_R = (kernel_R + 1) * exp(-beta * dt) when event accepted.
  kernel_R <- 0       # running sum of exp(-beta * (t_last - t_i)) for accepted events
  t_last_accepted <- time_window[1]  # time of last accepted event (or start)
  beta_overall <- params$beta_overall

  while (nrow(event_queue) > 0) {
    t <- proc.time()
    current_event <- event_queue[1, ,drop = FALSE]
    event_queue <- event_queue[-1, ,drop = FALSE]  # Remove the processed event

    if(is.null(current_net %v% 'n')){
      mark_sample <- PMF_mark(time = current_event$time,
                              params = params,
                              mark_filtration = current_net,
                              mark = NULL,
                              generate_mark = TRUE,
                              new_edge_hash = NULL,
                              stop_on_full_network = stop_on_full_network,
                              ...)
      net <- mark_sample$mark_sample
      accept <- 1
    }
    else{
      # get the mark samples
      mark_sample <- PMF_mark(time = current_event$time,
                              params = params,
                              mark_filtration = current_net,
                              mark = NULL,
                              generate_mark = TRUE,
                              generate_density = FALSE,
                              new_edge_hash = NULL,
                              stop_on_full_network = stop_on_full_network,
                              ...)
      net <- mark_sample$mark_sample
      if(hashed_edges && length(net$mel)!=0){
        # hash the network edge list for fast lookup:
        edges <- network::as.edgelist(net)
        keys_vec <- paste(edges[,1], edges[,2], sep = "-")
        edge_hash <- hash::hash(keys = keys_vec, values = rep(TRUE, length(keys_vec)))
      }else{
        edge_hash <- NULL
      }
      
      if(joint_accept){
        # Joint acceptance: need full conditional intensity (mark density * ground intensity)
        # Get mu_at_t for inhomogeneous case
        mu_at_t <- if (use_inhom) {
          pmax(mu_fun(current_event$time), 1e-12)
        } else { NULL }
        if (use_inhom) {
          intensity <- cond_intensity_inhom(new_net = net,
                                           t = current_event$time,
                                           mark_filtration = current_net,
                                           PMF_mark = PMF_mark,
                                           params = params,
                                           mu_at_t = mu_at_t,
                                           new_edge_hash = edge_hash,
                                           ...
          )$result
        } else {
          intensity <- cond_intensity(new_net = net,
                                     t = current_event$time,
                                     mark_filtration = current_net,
                                     PMF_mark = PMF_mark,
                                     params = params,
                                     new_edge_hash = edge_hash,
                                     ...
          )$result
        }
      }else{
        if(!is.null(n_mark_sample)){
          # Importance sampling path: need full cond_intensity
          mu_at_t <- if (use_inhom) {
            pmax(mu_fun(current_event$time), 1e-12)
          } else { NULL }
          imp_sample <- sapply(1:n_mark_sample,function(i){
            mark_sample <- PMF_mark(time = current_event$time,
                                    params = params,
                                    mark_filtration = current_net,
                                    mark = NULL,
                                    generate_mark = TRUE,
                                    new_edge_hash = TRUE,
                                    stop_on_full_network = stop_on_full_network,
                                    ...
            )
            net <- mark_sample$mark_sample
            if(hashed_edges){
              edges <- network::as.edgelist(net)
              keys_vec <- paste(edges[,1], edges[,2], sep = "-")
              edge_hash <- hash::hash(keys = keys_vec, values = rep(TRUE, length(keys_vec)))
            }else{
              edge_hash <- NULL
            }
            if (use_inhom) {
              intensity <- cond_intensity_inhom(new_net = net,
                                               t = current_event$time,
                                               mark_filtration = current_net,
                                               PMF_mark = PMF_mark,
                                               params = params,
                                               mu_at_t = mu_at_t,
                                               new_edge_hash = edge_hash,
                                               ...
              )
            } else {
              intensity <- cond_intensity(new_net = net,
                                         t = current_event$time,
                                         mark_filtration = current_net,
                                         PMF_mark = PMF_mark,
                                         params = params,
                                         new_edge_hash = edge_hash,
                                         ...
              )
            }
            return(intensity$result/mark_sample$mark_sample_density)
          })
          intensity <- mean(imp_sample)
        }else{
          # === FAST PATH: ground intensity only (no PMF_mark_CS recomputation) ===
          # For non-joint thinning, acceptance uses ground intensity = mu + K * kernel_sum.
          # The mark is already sampled above; no need to call cond_intensity (which would
          # redundantly call PMF_mark again just to compute mark density we don't use).
          # Use O(1) kernel recurrence instead of O(N) sum.
          dt <- current_event$time - t_last_accepted
          kernel_sum_at_t <- kernel_R * exp(-beta_overall * dt)
          mu_ground <- if (use_inhom) pmax(mu_fun(current_event$time), 1e-12) else params$mu
          intensity <- mu_ground + params$K * kernel_sum_at_t
        }
        # Use the same proposed mark for acceptance and for updating (do not resample)
      }

      accept <- intensity/lambda
    }
    n_proposed <- n_proposed + 1L
    accept_probs_buf[n_proposed] <- accept

    # if we accept the point add it in
    if(verbose){
      print(paste0("Number of edges proposed is ", length(net$mel)))
      print(paste0("Number of nodes proposed is ", net %n% "n"))
      print(paste0("accept prob is: ",accept))
      
      }
    if(runif(1) < accept){
      if(verbose){
        print('accepted!')
      }
      current_net <- net
      n_accepted <- n_accepted + 1L
      event_times_buf[n_accepted] <- current_event$time
      # --- Update O(1) kernel recurrence ---
      dt_acc <- current_event$time - t_last_accepted
      kernel_R <- (kernel_R + 1) * exp(-beta_overall * dt_acc)
      t_last_accepted <- current_event$time
      if(n_accepted > 2L){
        n_mark_dens <- n_mark_dens + 1L
        mark_density_buf[n_mark_dens] <- mark_sample$mark_density
      }
    }
    if(verbose){
      print(paste0("time is ",current_event$time, " size of net is ",current_net %n% 'n',' number of edges is ',length(current_net$mel)))
      print(paste0("time is ",current_event$time, " this iteration of while loop took ", round((proc.time()-t)[3],2)," seconds"))
    }
    # Concatenate new events to event_queue only if we have only one event left to go
    event_queue <- rbindlist(list(event_queue, rbindlist(new_events_list, use.names = TRUE)))
    new_events_list <- vector("list", length(new_events_list))  # Reset the list
    list_index <- 1
  }
  t1 <- proc.time() - t1
  print(paste0("simulation took ",round(t1[3],2)," seconds"))
  # Trim pre-allocated buffers to actual size
  events <- list(
    n = n_accepted,
    t = event_times_buf[seq_len(n_accepted)],
    mark_density = mark_density_buf[seq_len(n_mark_dens)]
  )
  accept_probs <- accept_probs_buf[seq_len(n_proposed)]
  return(list(events = events,
              net = current_net,
              accept_probs = accept_probs))
}

#' Log-likelihood for the Hawkes network growth model
#'
#' Computes the log-likelihood for a given parameter vector and observed network.
#' Supports both homogeneous (constant \code{mu}) and inhomogeneous (time-varying
#' \code{mu(t)}) background rates.  For an inhomogeneous background, supply
#' \code{mu_vec} and \code{integral_bg} (from
#' \code{\link{prepare_inhomogeneous_background}}).
#'
#' @param params List of parameters (\code{mu}, \code{beta_overall}, \code{K}, etc.).
#' @param time_window Numeric \code{c(t0, t1)}.
#' @param mark_filtration Observed network (filtration).
#' @param PMF_mark Mark PMF function (e.g. \code{\link{PMF_mark_BA}} or \code{\link{PMF_mark_CS}}).
#' @param mu_vec Optional numeric vector of background rates at each event time
#'   (same length as the event times extracted from \code{mark_filtration}).
#'   When \code{NULL} (default), a homogeneous rate \code{params$mu} is used.
#' @param integral_bg Optional scalar: integral of the background rate over the
#'   observation window.
#'   Required (non-\code{NULL}) whenever \code{mu_vec} is supplied.
#' @param edge_hash_list Optional list of edge hashes per event (for internal use).
#' @param verbose Print timing (default \code{FALSE}).
#' @param intens_funcs Precomputed intensity functions (for fitting).
#' @param ... Passed to \code{\link{cond_intensity}} or \code{\link{cond_intensity_inhom}}
#'   and to \code{PMF_mark} (e.g. \code{truncation}, \code{formula_RHS}, \code{cores},
#'   \code{combine_intensity}, \code{parallel_type}).
#' @return List with \code{loglik} (scalar) and \code{intens_funcs} (list of closures).
#' @seealso \code{\link{fit_hawkesNet}}, \code{\link{prepare_inhomogeneous_background}},
#'   \code{\link{cond_intensity}}, \code{\link{cond_intensity_inhom}}
#' @examples
#' \donttest{
#' params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
#' set.seed(1)
#' sim <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
#'                      verbose = FALSE, mu_multiplier = 5, truncation = 30)
#' ll <- loglik_hawkesNet(params, c(0, 3), sim$net, PMF_mark_BA, truncation = 30)
#' ll$loglik
#' }
#' @rdname loglik_hawkesNet
#' @export
loglik_hawkesNet = function(params,
                                  time_window,
                                  mark_filtration,
                                  PMF_mark,
                                  mu_vec = NULL,
                                  integral_bg = NULL,
                                  edge_hash_list = NULL,
                                  verbose = FALSE,
                                  intens_funcs = NULL,
                                  ...
){
  t<-proc.time()
  use_inhom <- !is.null(mu_vec) && !is.null(integral_bg)
  # don't allow negative parameters in first 2 (homogeneous only)
  if(!use_inhom && any(sapply(params[1:min(length(params),4)],function(x){x<0}))){
    return(list(loglik = -(10**(100)),
                intens_funcs = NULL))
  }
  times <- get_times(mark_filtration)
  times <- times$times

  tval <- time_window[2]-time_window[1]
  max_t <- max(times)
  if(tval < max(times) - min(times)){
    stop("realization has points outside time window")
  }
  if(use_inhom && length(mu_vec) != length(times)){
    stop("mu_vec must have same length as event times")
  }
  
  if(is.null(intens_funcs)){
    # Build intensity cache
    times_precalc <- list(times = times)
    dot_args <- list(...)
    formula_rhs <- dot_args$formula_RHS
    cores <- dot_args$cores
    parallel_type <- if (!is.null(dot_args$parallel_type)) dot_args$parallel_type else "auto"
    use_parallel <- !is.null(cores) && is.numeric(cores) && cores > 1
    combine_intensity <- use_inhom && isTRUE(dot_args$combine_intensity)
    # Reuse one ERNM model when running sequentially (avoids createCppModel per event)
    shared_model <- NULL
    if (!use_parallel && !is.null(formula_rhs)) {
      g0 <- network::network.initialize(0L, directed = FALSE)
      if ("na" %in% network::list.vertex.attributes(g0)) network::delete.vertex.attribute(g0, "na")
      shared_model <- createCppModel(as.formula(paste("g0 ~ ", formula_rhs)))
      shared_model$setNetwork(ernm::as.BinaryNet(g0))
    }
    if (use_parallel) message("Intensity cache: using ", cores, " cores")
    else message("Intensity cache: using 1 core")
    if (combine_intensity) message("Intensity cache: will combine CS closures into one (saves closure envs)")
    # Materialise ... into a concrete list so the closure serialises cleanly
    # for PSOCK workers (promises from ... cannot survive serialisation).
    extra_args <- dot_args
    extra_args[c("cores", "formula_RHS", "combine_intensity",
                 "parallel_type", "cache_intensity")] <- NULL
    intens_func <- function(i){
      # Log start of task in child
      # Only log for a subset of tasks to avoid flooding
      should_log <- (i == 1L || i == length(times) || (i %% 50 == 0))
      if (should_log) {
        cat(sprintf("  [intens_func] Task %d/%d starting (pid %d)\n", i, length(times), Sys.getpid()), file = stderr())
      }
      
      current_net <- filtration_to_net(mark_filtration, times[i], equals = TRUE)
      model <- if (!is.null(shared_model)) shared_model else if (!is.null(formula_rhs)) {
        createCppModel(as.formula(paste("current_net ~ ", formula_rhs)))
      } else NULL
      if (use_inhom) {
        call_args <- c(
          list(new_net = current_net, t = times[i],
               mark_filtration = current_net, PMF_mark = PMF_mark,
               params = params, mu_at_t = mu_vec[i], model = model,
               times = times_precalc,
               return_combined_inputs = combine_intensity),
          extra_args)
        intensity <- do.call(cond_intensity_inhom, call_args)
        out <- list(result = intensity$result, func = intensity$func)
        if (combine_intensity) {
          out$combined_inputs <- intensity$combined_inputs
          out$diffs_kernel <- intensity$diffs
        }
      } else {
        call_args <- c(
          list(new_net = current_net, t = times[i],
               mark_filtration = current_net, PMF_mark = PMF_mark,
               params = params, model = model, times = times_precalc),
          extra_args)
        intensity <- do.call(cond_intensity, call_args)
        out <- list(result = intensity$result, func = intensity$func)
      }
      
      if (should_log) {
        cat(sprintf("  [intens_func] Task %d/%d complete\n", i, length(times)), file = stderr())
      }
      out
    }
    if(use_parallel){
      message("starting intens list calculation")
      t <- proc.time()
      # Process last events first (biggest networks) so slow jobs start first.
      # Use mc.preschedule = TRUE for stability: avoids per-task fork management
      # that can cause mclapply to hang on the second call in the same session.
      intens_list <- safe_parallel_lapply(
        rev(seq_along(times)), intens_func,
        mc.cores = cores, mc.preschedule = TRUE,
        parallel_type = parallel_type
      )
      intens_list <- rev(intens_list)
      message("intens list ", round((proc.time()-t)[3],2)," seconds")
      
      # Check if any tasks failed or returned NULL (indicates fork issues)
      if (is.list(intens_list)) {
        failed_tasks <- vapply(intens_list, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
        if (any(failed_tasks)) {
          message(sprintf("  [loglik] WARNING: %d tasks failed during intensity calculation. Check for OOM or deadlocks.", sum(failed_tasks)))
        }
      }
      
      # Clean up any zombie child processes immediately after parallel call.
      # stale pipes/signal handlers from the first call can interfere with the second.
      tryCatch({
        children_fn  <- get("children",  envir = asNamespace("parallel"))
        mccollect_fn <- get("mccollect", envir = asNamespace("parallel"))
        while (length(children_fn()) > 0L) {
          mccollect_fn(wait = FALSE, timeout = 1)
        }
      }, error = function(e) NULL)
      gc()
    }else{
      t1 <- proc.time()
      intens_list <- lapply(seq_along(times), intens_func)
      message("intens list took ", round((proc.time()-t1)[3],2)," seconds")
    }
    intens_vec <- vapply(intens_list, function(x) x$result, numeric(1))
    intens_funcs <- lapply(intens_list, function(x) x$func)
    # Combine per-event closures into one vectorized closure (inhom + CS only)
    if (combine_intensity) {
      ci <- lapply(intens_list, function(x) x$combined_inputs)
      dk <- lapply(intens_list, function(x) x$diffs_kernel)
      # Free the massive intens_list BEFORE combine (each element holds a full
      # network filtration + ERNM model in its closure env). Without this, the
      # parent process stays bloated and any subsequent mclapply fork inherits
      # the bloated address space, causing hangs or OOM on the second fit.
      rm(intens_list); gc()
      if (all(vapply(ci, function(x) !is.null(x), NA))) {
        combined <- build_combined_intensity_funcs(ci, dk, mu_vec, times)
        if (!is.null(combined)) {
          intens_funcs <- combined
          message("Intensity cache: combined ", length(times),
                  " closures into 1 (memory: one list of change_stats)")
        }
      }
      # Free combine intermediates
      rm(ci, dk)
    } else {
      # Even without combine, free intens_list (per-event closures stay in intens_funcs)
      rm(intens_list)
    }
    # Reclaim memory so subsequent parallel calls don't fork a bloated process
    gc()
  }else{
    # Evaluate pre-cached closures
    t1 <- proc.time()
    if (length(intens_funcs) == 1L) {
      v <- intens_funcs[[1]](params)
      intens_vec <- if (length(v) > 1L) v else v[1]
    } else {
      intens_vec <- numeric(length(intens_funcs))
      for (i in seq_along(intens_funcs)) intens_vec[i] <- intens_funcs[[i]](params)
    }
    if(verbose){
      print(paste0("evaluating intens list with intens funcs took ", round((proc.time()-t1)[3],2)," seconds"))
    }
  }

  tmp <- intens_vec
  if(any(is.na(tmp))){
    tmp[is.na(tmp)] <- min(tmp[!is.na(tmp)])/2
  }
  
  if(sum(tmp<=0)!=0){
    warning("some of the intens lists have zero")
    tmp[tmp<=0] <- min(tmp[tmp>0])/2
  }
  intens_sum <- sum(log(tmp))

  # Integral (compensator):
  max_t <- max(times)
  pieces <- 1 - exp(-params$beta_overall * (tval - times))
  # Guard against division by zero or very small beta_overall
  if (!is.finite(params$beta_overall) || params$beta_overall <= 0 || params$beta_overall < 1e-10) {
    return(list(loglik = -1e10, intens_funcs = intens_funcs))
  }
  kernel_integral <- (1/params$beta_overall)*params$K*sum(pieces)
  integral <- if (use_inhom) integral_bg + kernel_integral else params$mu * tval + kernel_integral
  loglik <- intens_sum - integral
  
  # Guard: ensure loglik is finite for L-BFGS-B
  if (!is.finite(loglik) || !is.finite(intens_sum) || !is.finite(integral)) {
    if (verbose) {
      warning("loglik_hawkesNet: non-finite loglik detected. intens_sum=", intens_sum, 
              ", integral=", integral, ", returning -1e10")
    }
    return(list(loglik = -1e10, intens_funcs = intens_funcs))
  }
  
  t<-proc.time() - t
  if(verbose){
    print(paste0("this iteration of loglik took ", round(t[3],2)," seconds"))
  }
  
  return(list(loglik = loglik,
              intens_funcs = intens_funcs))
}

#' Fit the Hawkes network growth model by maximum likelihood
#'
#' Estimates parameters of a Hawkes-driven network growth model via numerical
#' optimisation of the log-likelihood.
#' Supports both homogeneous (constant \code{mu}) and inhomogeneous (time-varying
#' \code{mu(t)}) background rates; supply \code{mu_vec} and \code{integral_bg}
#' (from \code{\link{prepare_inhomogeneous_background}}) for the latter.
#'
#' @param params_init List of initial parameter values.
#' @param time_window Numeric \code{c(t0, t1)}.
#' @param mark_filtration Observed network.
#' @param PMF_mark Mark PMF function (e.g. \code{\link{PMF_mark_BA}} or \code{\link{PMF_mark_CS}}).
#' @param maxit Maximum number of iterations for the optimizer.
#' @param trace Trace level (default 0).
#' @param REPORT Reporting interval for \code{optim} (default 10).
#' @param reltol Relative convergence tolerance (default 1e-8).
#' @param parscale Scale vector for parameters (default all 1).
#' @param fixed_params Character vector of parameter names to hold fixed.
#' @param cache_intensity If \code{TRUE}, cache intensity functions (default \code{TRUE}).
#' @param combine_intensity If \code{TRUE} (default) and fitting with an
#'   inhomogeneous background (\code{mu_vec} supplied), combine per-event closures
#'   into one vectorized closure for faster evaluation.
#' @param method Optimization method: \code{"Nelder-Mead"} (default, derivative-free) or
#'   \code{"L-BFGS-B"} (gradient-based with box constraints).
#' @param verbose Print diagnostics to \code{stderr} (default \code{TRUE}).
#' @param parallel_type Parallelization strategy for intensity cache build:
#'   \code{"auto"} (default), \code{"psock"}, or \code{"fork"}.
#' @param mu_vec Optional numeric vector of background rates at each event time.
#'   When supplied together with \code{integral_bg}, the model uses an
#'   inhomogeneous background rate instead of constant \code{params$mu}.
#' @param integral_bg Optional scalar: integral of the background rate over
#'   the observation window.  Required when \code{mu_vec} is supplied.
#' @param ... Passed to \code{\link{loglik_hawkesNet}} and \code{PMF_mark}
#'   (e.g. \code{truncation}, \code{formula_RHS}, \code{cores}).
#' @return List with:
#'   \describe{
#'     \item{fit}{Output of \code{\link[stats]{optim}}.}
#'     \item{intens_funcs}{Cached intensity closures (can be large; NULL out and \code{gc()} when done).}
#'     \item{params_init_old}{Original \code{params_init} (for relist / GOF).}
#'     \item{fit_table}{Data frame of parameter estimates and standard errors.}
#'     \item{hessian}{Numerical Hessian of negative log-likelihood at MLE.}
#'   }
#' @seealso \code{\link{loglik_hawkesNet}}, \code{\link{sim_hawkesNet}},
#'   \code{\link{prepare_inhomogeneous_background}}, \code{\link{gof}}
#' @examples
#' \donttest{
#' params <- list(mu = 0.5, beta_overall = 1, K = 0.3, beta_edges = 0.5, m = 1)
#' set.seed(1)
#' sim <- sim_hawkesNet(params, c(0, 3), PMF_mark_BA, cond_intensity,
#'                      verbose = FALSE, mu_multiplier = 5, truncation = 30)
#' fit <- fit_hawkesNet(params, c(0, 3), sim$net, PMF_mark_BA,
#'                      maxit = 50, verbose = FALSE, truncation = 30)
#' fit$fit_table
#' }
#' @rdname fit_hawkesNet
#' @export
fit_hawkesNet <- function(params_init,
                                time_window,
                                mark_filtration,
                                PMF_mark,
                                maxit,
                                trace = 0,
                                REPORT = 10,
                                reltol = 1e-8,
                                parscale = NULL,
                                fixed_params = NULL,
                                cache_intensity = TRUE,
                                combine_intensity = TRUE,
                                method = "Nelder-Mead",
                                verbose = TRUE,
                                parallel_type = "auto",
                                mu_vec = NULL,
                                integral_bg = NULL,
                                ...){
  use_inhom <- !is.null(mu_vec) && !is.null(integral_bg)
  # Helper: write to stderr (unbuffered even inside optim's C code) and flush
  vcat <- function(...) if (verbose) { cat(..., file = stderr()); flush(stderr()) }

  params_init_old <- params_init
  # Shallow copy so stripping levels does not modify params_init_old (needed for loglik and relist restore)
  params_init <- as.list(params_init_old)
  if(!is.null(fixed_params)){
    for(k in fixed_params){
      params_init[[k]] <- NULL
    }
  }
  # Strip vertex_categorical_levels (character metadata, not numeric parameters)
  # so unlist() yields a purely numeric vector for optim.
  params_init$vertex_categorical_levels <- NULL
  
  if (is.null(parscale)) {
    flat_params <- unlist(params_init)
    parscale <- rep(1, length(flat_params))
  }
  
  # Validate that params match the mark PMF (required names and, for CS, CS_params length)
  validate_params_for_PMF(params_init_old, PMF_mark, mark_filtration, ...)
  
  t_fit_start <- proc.time()[3]
  
  # pre-calculate the param -> conditonal intensity mapping
  # since the observation never changes - not need to do expensive network processes every iteration
  # then param -> likelihood should be very fast
  # 1. Pre-calculate ONLY if cache_intensity is TRUE
  cached_funcs <- NULL
  if(cache_intensity){
    vcat("[fit] Pre-calculating intensity closures",
         if (use_inhom) " (inhomogeneous background)" else "", "...\n")
    t_cache_start <- proc.time()[3]
    init_lik <- loglik_hawkesNet(params = params_init_old,
                                     time_window = time_window,
                                     mark_filtration = mark_filtration,
                                     PMF_mark = PMF_mark,
                                     mu_vec = mu_vec,
                                     integral_bg = integral_bg,
                                     combine_intensity = combine_intensity,
                                     parallel_type = parallel_type,
                                     ...)
    cached_funcs <- init_lik$intens_funcs
    # Free init_lik immediately: it holds the full loglik result including
    # intens_vec and the log-likelihood. We only need cached_funcs.
    # This prevents the second fit's mclapply from forking a bloated process.
    rm(init_lik); gc()
    vcat("[fit] Intensity cache built: ", round(proc.time()[3] - t_cache_start, 1), " s\n")
  } else {
    vcat("[fit] Caching disabled (Safe Mode).\n")
    init_lik <- NULL
  }
  # --- Precompute values used every iteration ---
  times_cached <- get_times(mark_filtration)$times
  tval_cached  <- time_window[2] - time_window[1]
  is_combined  <- !is.null(cached_funcs) && length(cached_funcs) == 1L

  # --- Evaluation diagnostics (accumulate timing, report periodically) ---
  eval_env <- new.env(parent = emptyenv())
  eval_env$n_eval <- 0L
  eval_env$t_relist_total   <- 0
  eval_env$t_validate_total <- 0
  eval_env$t_closure_total  <- 0
  eval_env$t_cleanup_total  <- 0
  eval_env$t_integral_total <- 0
  eval_env$t_total_total    <- 0
  eval_env$t_last_report    <- proc.time()[3]
  eval_env$best_ll          <- -Inf
  DIAG_INTERVAL <- 10L  # report every N evaluations

  dot_args <- list(...)
  # Fast optim_func: calls cached closure directly when available, computes integral inline.
  optim_func <- function(params){
    t0 <- proc.time()[3]

    # 1. relist + restore fixed params
    t1 <- proc.time()[3]
    params_curr <- relist(params, skeleton = params_init)
    params_curr$vertex_categorical_levels <- params_init_old$vertex_categorical_levels
    if (!is.null(fixed_params)) {
      for (k in fixed_params) params_curr[[k]] <- params_init_old[[k]]
    }
    t_relist <- proc.time()[3] - t1

    # 2. validate
    t1 <- proc.time()[3]
    if (!point_process_params_valid(params_curr)) {
      eval_env$n_eval <- eval_env$n_eval + 1L
      return(-1e10)
    }
    t_validate <- proc.time()[3] - t1

    if (!is.null(cached_funcs)) {
      # Fast path: evaluate cached closures directly
      t1 <- proc.time()[3]
      intens_vec <- tryCatch({
        if (is_combined) cached_funcs[[1L]](params_curr)
        else {
          v <- numeric(length(cached_funcs))
          for (j in seq_along(cached_funcs)) v[j] <- cached_funcs[[j]](params_curr)
          v
        }
      }, error = function(e) NULL)
      t_closure <- proc.time()[3] - t1
      if (is.null(intens_vec)) {
        eval_env$n_eval <- eval_env$n_eval + 1L
        return(-1e10)
      }

      t1 <- proc.time()[3]
      bad <- !is.finite(intens_vec) | intens_vec <= 0
      if (any(bad)) intens_vec[bad] <- 1e-10
      intens_sum <- sum(log(intens_vec))
      t_cleanup <- proc.time()[3] - t1

      t1 <- proc.time()[3]
      b <- params_curr$beta_overall
      if (!is.finite(b) || b < 1e-10) {
        eval_env$n_eval <- eval_env$n_eval + 1L
        return(-1e10)
      }
      pieces <- 1 - exp(-b * (tval_cached - times_cached))
      kernel_int <- (1 / b) * params_curr$K * sum(pieces)
      integral <- if (use_inhom) integral_bg + kernel_int
                  else params_curr$mu * tval_cached + kernel_int
      t_integral <- proc.time()[3] - t1

      ll <- intens_sum - integral
      if (!is.finite(ll)) {
        eval_env$n_eval <- eval_env$n_eval + 1L
        return(-1e10)
      }

      # --- Accumulate timing ---
      t_total <- proc.time()[3] - t0
      eval_env$n_eval <- eval_env$n_eval + 1L
      eval_env$t_relist_total   <- eval_env$t_relist_total   + t_relist
      eval_env$t_validate_total <- eval_env$t_validate_total + t_validate
      eval_env$t_closure_total  <- eval_env$t_closure_total  + t_closure
      eval_env$t_cleanup_total  <- eval_env$t_cleanup_total  + t_cleanup
      eval_env$t_integral_total <- eval_env$t_integral_total + t_integral
      eval_env$t_total_total    <- eval_env$t_total_total    + t_total
      if (ll > eval_env$best_ll) eval_env$best_ll <- ll

      # --- Periodic report (first eval + every DIAG_INTERVAL) ---
      n <- eval_env$n_eval
      if (n == 1L || n %% DIAG_INTERVAL == 0L) {
        elapsed <- proc.time()[3] - eval_env$t_last_report
        ms_fn <- function(x) round(x / n * 1000, 2)
        if (n == 1L) {
          vcat(sprintf("  [eval 1] first eval: %.0f ms | ll=%.2f | closure=%.0f ms\n",
                    t_total * 1000, ll, t_closure * 1000))
        } else {
          batch <- min(n, DIAG_INTERVAL)
          vcat(sprintf("  [eval %d] %.1fs wall for last %d evals (%.0f ms/eval) | best_ll=%.2f\n",
                    n, elapsed, batch, elapsed / batch * 1000, eval_env$best_ll))
          vcat(sprintf("    avg (ms/eval): relist=%.2f validate=%.2f closure=%.2f cleanup=%.2f integral=%.2f total=%.2f\n",
                    ms_fn(eval_env$t_relist_total), ms_fn(eval_env$t_validate_total),
                    ms_fn(eval_env$t_closure_total), ms_fn(eval_env$t_cleanup_total),
                    ms_fn(eval_env$t_integral_total), ms_fn(eval_env$t_total_total)))
        }
        eval_env$t_last_report <- proc.time()[3]
      }

      return(ll)
    } else {
      # Slow path: full loglik_hawkesNet (cache disabled)
      result <- tryCatch({
        do.call(loglik_hawkesNet, c(
          list(params = params_curr,
               time_window = time_window,
               mark_filtration = mark_filtration,
               PMF_mark = PMF_mark,
               mu_vec = mu_vec,
               integral_bg = integral_bg,
               intens_funcs = NULL),
          dot_args
        ))
      }, error = function(e) {
        return(list(loglik = -1e10, intens_funcs = NULL))
      })
      ll <- result$loglik
      if (!is.finite(ll)) return(-1e10)
      eval_env$n_eval <- eval_env$n_eval + 1L
      return(ll)
    }
  }
  flat_par <- unlist(params_init)
  n_par <- length(flat_par)
  n_events_actual <- length(times_cached)

  if (!is.null(cached_funcs)) {
    if (is_combined) {
      vcat("[fit] Optimizing ", n_par, " params | 1 combined closure (",
           n_events_actual, " events, vectorized) | method=", method, " maxit=", maxit, "\n")
    } else {
      vcat("[fit] Optimizing ", n_par, " params | ", length(cached_funcs),
           " closures (", n_events_actual, " events, sequential) | method=", method,
           " maxit=", maxit, "\n")
      if (length(cached_funcs) > 20L) {
        vcat("[fit] *** WARNING: sequential eval with ", length(cached_funcs),
             " closures will be VERY slow! ***\n")
        vcat("[fit] *** Add combine_intensity = TRUE for ~1000x speedup. ***\n")
      }
    }
  }

  optim_args <- list(
    par = flat_par,
    fn = optim_func,
    method = method,
    control = list(fnscale = -1,
                   trace = trace,
                   maxit = maxit,
                   reltol = reltol,
                   parscale = parscale),
    hessian = TRUE
  )
  if (method == "L-BFGS-B") {
    bounds <- build_optim_bounds(names(flat_par))
    optim_args$lower <- bounds$lower
    optim_args$upper <- bounds$upper
    optim_args$control$REPORT <- REPORT
  }
  t_optim_start <- proc.time()[3]
  fit <- do.call(optim, optim_args)
  t_optim_elapsed <- round(proc.time()[3] - t_optim_start, 1)
  n_iter <- if (!is.null(fit$counts)) fit$counts[1L] else NA_integer_
  n_fneval <- if (!is.null(fit$counts) && length(fit$counts) >= 2L) fit$counts[2L] else NA_integer_
  vcat("[fit] Optimization done: ", t_optim_elapsed, " s | iterations: ", n_iter,
       if (is.finite(n_fneval)) paste0(" | fn evals: ", n_fneval) else "",
       " | s/iter: ", if (is.finite(n_iter) && n_iter > 0) round(t_optim_elapsed / n_iter, 2) else "n/a", "\n")
  # --- Final evaluation timing summary ---
  n_e <- eval_env$n_eval
  if (n_e > 0L && eval_env$t_total_total > 0) {
    ms <- function(x) round(x / n_e * 1000, 2)
    pct <- function(x) if (eval_env$t_total_total > 0) round(x / eval_env$t_total_total * 100, 1) else 0
    vcat(sprintf("[fit] Eval timing summary (%d evals, %.1f ms/eval avg):\n", n_e, eval_env$t_total_total / n_e * 1000))
    vcat(sprintf("  relist:   %6.2f ms/eval (%4.1f%%)\n", ms(eval_env$t_relist_total),   pct(eval_env$t_relist_total)))
    vcat(sprintf("  validate: %6.2f ms/eval (%4.1f%%)\n", ms(eval_env$t_validate_total), pct(eval_env$t_validate_total)))
    vcat(sprintf("  closure:  %6.2f ms/eval (%4.1f%%)\n", ms(eval_env$t_closure_total),  pct(eval_env$t_closure_total)))
    vcat(sprintf("  cleanup:  %6.2f ms/eval (%4.1f%%)\n", ms(eval_env$t_cleanup_total),  pct(eval_env$t_cleanup_total)))
    vcat(sprintf("  integral: %6.2f ms/eval (%4.1f%%)\n", ms(eval_env$t_integral_total), pct(eval_env$t_integral_total)))
    vcat(sprintf("  TOTAL:    %6.2f ms/eval           (wall: %.1f s)\n", ms(eval_env$t_total_total), eval_env$t_total_total))
  }
  vcat("[fit] Fitting total: ", round(proc.time()[3] - t_fit_start, 1), " s\n")

  # Results table: estimate and standard error (from numerical Hessian)
  par_names <- names(fit$par)
  fit_table <- data.frame(
    parameter = par_names,
    estimate  = fit$par,
    std.error = NA_real_,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  hessian <- fit$hessian
  if (!is.null(hessian)) {
    # optim returns hessian of fn (loglik), which is negative definite at a maximum.
    # The variance-covariance matrix is the inverse of the *negative* hessian (observed information).
    vcov <- tryCatch(solve(-hessian), error = function(e) NULL)
    if (!is.null(vcov)) {
        se <- sqrt(pmax(diag(vcov), 0))
        fit_table$std.error <- se
    }
  }
  # Replace CS_params1, CS_params2, ... with actual ERNM statistic names
  fit_table <- rename_CS_params_in_table(fit_table, mark_filtration, list(...))
  vcat("[fit] Results:\n")
  if (verbose) print(fit_table, max = NULL, row.names = TRUE)
  if (all(is.na(fit_table$std.error))) {
    vcat("(Standard errors not available; Hessian inversion failed.)\n")
  }

  # NOTE: intens_funcs can be very large (closure environments with stacked matrices
  # for all events). If running multiple fits sequentially, NULL out intens_funcs and
  # call gc() before the next fit to prevent fork()/PSOCK memory bloat.
  list(
    fit = fit,
    intens_funcs = cached_funcs,
    params_init_old = params_init_old,
    fit_table = fit_table,
    hessian = hessian
  )
}


#' Compensators (integrated intensity) for the Hawkes growth model
#'
#' @param params List of parameters.
#' @param time_window Numeric \code{c(t0, t1)}.
#' @param mark_filtration Observed network.
#' @return Numeric vector of compensator values at each event time.
#' @export
compensators_hawkesNet <- function(params,
                                         time_window,
                                         mark_filtration){
  times <- get_times(mark_filtration)
  times <- times$times

  tval <- time_window[2] - time_window[1]
  max_t <- max(times)
  if (tval < max(times) - min(times)) {
    stop("realization has points outside time window")
  }
  pieces <- 1 - exp(-params$beta_overall * (tval - times))
  incremental <- sapply(seq_along(times), function(i) {
    params$mu * times[i] + (1/params$beta_overall)*params$K*sum(pieces[1:i])
  })
  return(incremental)
}

#' Kolmogorov-Smirnov test p-value for time rescaling (Hawkes growth model)
#'
#' @param params List of parameters.
#' @param time_window Numeric \code{c(t0, t1)}.
#' @param mark_filtration Observed network.
#' @return P-value of the KS test under the null that rescaled times are uniform.
#' @export
ks_test_pval_hawkesNet <- function(params,
                                         time_window,
                                         mark_filtration){
  compensators <- compensators_hawkesNet(params = unlist(params),
                                               time_window = time_window,
                                               mark_filtration = mark_filtration
                                               )
  compensator_incs <- diff(compensators)
  test_dist <- 1 - exp(-compensator_incs)
  test <- ks.test(test_dist,"punif")
  # hist(test_dist)
  # print(test$p.value)
  return(test$p.value)
}


# =============================================================================
# Inhomogeneous background rate (KDE) support
# =============================================================================
# Replaces constant mu with time-varying mu(t) from KDE. Intensity:
# lambda(t) = mu(t) + K * sum(decays). Integral: int_0^T mu(s) ds + (1/beta)*K*sum(pieces).
# =============================================================================

#' Conditional intensity at one event time with inhomogeneous background
#'
#' Same as cond_intensity but uses mu_at_t (scalar) instead of params$mu.
#' @param new_net Network state at current time
#' @param t Current event time
#' @param mark_filtration Full mark filtration
#' @param PMF_mark Mark PMF function (e.g. PMF_mark_CS)
#' @param params Parameter list (beta_overall, K, CS_params, etc.; mu not used)
#' @param mu_at_t Scalar background rate at time t (e.g. from KDE)
#' @param new_edge_hash Optional hash of existing edges
#' @param times Precomputed times from get_times; if NULL, taken from mark_filtration
#' @param ... Passed to PMF_mark (e.g. formula_RHS, truncation)
#' @return List with result (intensity), func (function to evaluate intensity at new params), lambda, kernel_sum, decays, diffs
#' @export
cond_intensity_inhom <- function(new_net,
                                 t,
                                 mark_filtration = NULL,
                                 PMF_mark,
                                 params,
                                 mu_at_t,
                                 new_edge_hash = NULL,
                                 times = NULL,
                                 ...) {
  tmp <- PMF_mark(time = t,
                  params = params,
                  mark_filtration = mark_filtration,
                  mark = new_net,
                  generate_mark = FALSE,
                  new_edge_hash = new_edge_hash,
                  ...)
  if (is.null(times)) times <- get_times(mark_filtration)
  tt <- times$times
  tt <- tt[tt < t]
  diffs <- t - tt

  log_mark_density0 <- tmp$log_mark_density
  log_density_func  <- tmp$log_density_func

  decays0 <- exp(-params$beta_overall * diffs)
  log_result0 <- log_mark_density0 + log(mu_at_t + params$K * sum(decays0))
  result0 <- exp(log_result0)

  e_tiny <- new.env(parent = baseenv())
  e_tiny$diffs_local <- diffs
  e_tiny$ldf_local   <- log_density_func
  e_tiny$mu_at_t     <- mu_at_t

  func_template <- function(params) {
    decays <- exp(-params$beta_overall * diffs_local)
    exp(ldf_local(params) + log(mu_at_t + params$K * sum(decays)))
  }
  func <- func_template
  environment(func) <- e_tiny

  out <- list(
    result = result0,
    func   = func,
    lambda = mu_at_t,
    kernel_sum = params$K * sum(decays0),
    decays = decays0,
    diffs  = diffs
  )
  if (!is.null(tmp$combined_inputs)) out$combined_inputs <- tmp$combined_inputs
  out
}


#' Build one combined intensity closure from per-event CS ingredients
#'
#' Used when \code{combine_intensity = TRUE}: stacks change_stats into one matrix,
#' one matrix multiply per evaluation (vectorized); per-event steps (plogis, decay,
#' in_mark, node_dens, kernel) stay in a loop. Faster than N separate multiplies.
#' @param combined_inputs_list List of length N of combined_inputs from PMF_mark_CS (each has change_stats, in_mark, diffs, ...).
#' @param diffs_kernel_list List of length N of kernel diffs (t_i - times[times < t_i]).
#' @param mu_vec Length-N background rate at each event time.
#' @param times Length-N event times (used only for length check).
#' @return List of one function \code{f(params)} returning numeric vector of length N, or NULL on error.
#' @noRd
build_combined_intensity_funcs <- function(combined_inputs_list, diffs_kernel_list, mu_vec, times) {
  N <- length(times)
  if (N == 0L || length(combined_inputs_list) != N || length(diffs_kernel_list) != N || length(mu_vec) != N) {
    return(NULL)
  }
  n_cs <- NA_integer_
  for (i in seq_len(N)) {
    inp <- combined_inputs_list[[i]]
    if (!is.null(inp) && !isTRUE(inp$degenerate_edges) && nrow(inp$change_stats) > 0L) {
      n_cs <- ncol(inp$change_stats)
      break
    }
  }
  if (!is.finite(n_cs) || n_cs < 1L) return(NULL)

  # --- Pre-stack everything for fully vectorized evaluation ---
  nrows <- vapply(combined_inputs_list, function(x) {
    if (is.null(x) || isTRUE(x$degenerate_edges)) 0L else nrow(x$change_stats)
  }, 0L)
  row_bounds <- c(0L, cumsum(nrows))    # length N+1; segment i = (row_bounds[i]+1):row_bounds[i+1]
  total_rows <- row_bounds[N + 1L]
  if (total_rows == 0L) return(NULL)

  # Stacked change_stats matrix (already existed)
  change_stats_stacked <- do.call(rbind, lapply(seq_len(N), function(i) {
    inp <- combined_inputs_list[[i]]
    if (nrows[i] == 0L) matrix(0, 0, n_cs) else inp$change_stats
  }))

  # Stacked edge-decay diffs and in_mark logical
  diffs_stacked <- unlist(lapply(seq_len(N), function(i) {
    inp <- combined_inputs_list[[i]]
    if (nrows[i] == 0L) numeric(0) else inp$diffs
  }), use.names = FALSE)

  in_mark_stacked <- unlist(lapply(seq_len(N), function(i) {
    inp <- combined_inputs_list[[i]]
    if (nrows[i] == 0L) logical(0) else inp$in_mark
  }), use.names = FALSE)

  # Per-event constant vectors
  degenerate <- nrows == 0L
  new_minus_old <- vapply(seq_len(N), function(i) {
    inp <- combined_inputs_list[[i]]
    if (is.null(inp) || degenerate[i]) 0L else as.integer(inp$new_nodes - inp$old_nodes)
  }, 0L)
  past_max_node_time <- vapply(seq_len(N), function(i) {
    inp <- combined_inputs_list[[i]]
    if (is.null(inp) || is.null(inp$max_node_time)) FALSE
    else !is.null(inp$time) && inp$time > inp$max_node_time
  }, FALSE)

  # Event times for kernel recurrence
  event_times <- times

  # Segment boundary indices for cumsum segment-sum trick
  seg_end   <- row_bounds[-1]        # row_bounds[2:(N+1)]
  seg_start <- row_bounds[-(N + 1L)] # row_bounds[1:N]

  # --- Vertex categorical: pre-stack observed level indices ---
  obs_cat_attr_name  <- NULL
  obs_cat_stacked    <- NULL
  obs_cat_seg_bounds <- NULL   # length N+1

  # Find first event with observed_categorical to get attribute name and level names
  for (i in seq_len(N)) {
    inp <- combined_inputs_list[[i]]
    if (!is.null(inp) && length(inp$observed_categorical) > 0L) {
      obs_cat_attr_name <- names(inp$observed_categorical)[1L]
      break
    }
  }
  if (!is.null(obs_cat_attr_name)) {
    # Find level names
    level_names_cat <- NULL
    for (i in seq_len(N)) {
      ln <- combined_inputs_list[[i]]$level_names_by_attr[[obs_cat_attr_name]]
      if (!is.null(ln)) { level_names_cat <- ln; break }
    }
    if (!is.null(level_names_cat)) {
      unknown_idx <- match("unknown", level_names_cat)
      if (is.na(unknown_idx)) unknown_idx <- 1L
      # Pre-compute level indices per event
      idx_list <- vector("list", N)
      obs_lengths <- integer(N)
      for (i in seq_len(N)) {
        inp <- combined_inputs_list[[i]]
        if (!is.null(inp) && !is.null(inp$observed_categorical) &&
            obs_cat_attr_name %in% names(inp$observed_categorical)) {
          obs_vals <- inp$observed_categorical[[obs_cat_attr_name]]
          idx <- match(obs_vals, level_names_cat)
          idx[is.na(idx)] <- unknown_idx
          idx_list[[i]] <- idx
          obs_lengths[i] <- length(idx)
        } else {
          idx_list[[i]] <- integer(0)
          obs_lengths[i] <- 0L
        }
      }
      obs_cat_stacked <- unlist(idx_list, use.names = FALSE)
      obs_cat_seg_bounds <- c(0L, cumsum(obs_lengths))
    }
  }

  # --- Free per-event data that is now stacked ---
  for (i in seq_len(N)) {
    combined_inputs_list[[i]]$change_stats <- NULL
    combined_inputs_list[[i]]$diffs <- NULL
    combined_inputs_list[[i]]$in_mark <- NULL
    combined_inputs_list[[i]]$observed_categorical <- NULL
    combined_inputs_list[[i]]$level_names_by_attr <- NULL
  }
  rm(combined_inputs_list, diffs_kernel_list)

  eps <- 1e-10

  # ====================================================================
  # Fully vectorized closure: O(N) kernel + O(total_rows) mark density
  # ====================================================================
  combined_closure <- function(params) {
    # 1. Kernel sums via O(N) recurrence (Hawkes trick):
    #    R[i] = sum_{j<i} exp(-beta * (t_i - t_j))
    #         = (R[i-1] + 1) * exp(-beta * dt[i])
    beta <- params$beta_overall
    K_par <- params$K
    R <- numeric(N)
    if (N > 1L) {
      for (i in 2L:N) {
        R[i] <- (R[i - 1L] + 1) * exp(-beta * (event_times[i] - event_times[i - 1L]))
      }
    }
    kernel_sums <- K_par * R

    # 2. Vectorized edge probabilities: one matmul + one plogis + one exp
    eta_all <- as.vector(change_stats_stacked %*% params$CS_params)
    p_all   <- stats::plogis(eta_all) * exp(-params$beta_edges * diffs_stacked)
    p_all   <- pmin(pmax(p_all, eps), 1 - eps)

    # 3. Segment log-sums via cumsum (replaces per-event loop)
    log_contrib <- ifelse(in_mark_stacked, log(p_all), log1p(-p_all))
    # Guard: a single NaN would poison cumsum for all subsequent events.
    # Replace non-finite entries with 0 (neutral for sum) so only that event's
    # segment is affected, matching the old per-event anyNA guard.
    bad_lc <- !is.finite(log_contrib)
    if (any(bad_lc)) log_contrib[bad_lc] <- 0
    cs_log <- c(0, cumsum(log_contrib))
    log_edge_sums <- cs_log[seg_end + 1L] - cs_log[seg_start + 1L]  # length N

    # 4. Node density: vectorized dpois (log = TRUE for stability)
    node_dens <- rep(0, N)
    needs_node <- !past_max_node_time & !degenerate
    if (any(needs_node)) {
      node_dens[needs_node] <- stats::dpois(new_minus_old[needs_node],
                                            params$node_lambda, log = TRUE)
      bad <- !is.finite(node_dens)
      if (any(bad)) node_dens[bad] <- -1e10
    }

    # 5. Vertex categorical: vectorized via pre-stacked level indices
    if (!is.null(obs_cat_stacked) && length(obs_cat_stacked) > 0L) {
      vcat <- params$vertex_categorical
      if (!is.null(vcat) && is.list(vcat) && obs_cat_attr_name %in% names(vcat)) {
        p_n1   <- as.numeric(vcat[[obs_cat_attr_name]])
        p_full <- pmax(c(p_n1, 1 - sum(p_n1)), eps)
        log_p  <- log(p_full)
        log_p_obs <- log_p[obs_cat_stacked]
        cs_cat <- c(0, cumsum(log_p_obs))
        cat_sums <- cs_cat[obs_cat_seg_bounds[-1] + 1L] -
                    cs_cat[obs_cat_seg_bounds[-(N + 1L)] + 1L]
        node_dens <- node_dens + cat_sums
      }
    }

    # 6. Combine: intensity[i] = exp(log_mark_density[i]) * (mu[i] + kernel[i])
    log_dens <- log_edge_sums + node_dens
    out <- exp(log_dens) * (mu_vec + kernel_sums)
    out[degenerate] <- mu_vec[degenerate] + kernel_sums[degenerate]
    pmax(out, eps)
  }

  list(combined_closure)
}


#' Fit HawkesGrowthNet with inhomogeneous (KDE) background
#'
#' Backward-compatible wrapper around \code{\link{fit_hawkesNet}}.
#' Equivalent to calling
#' \code{fit_hawkesNet(..., mu_vec = mu_vec, integral_bg = integral_bg)}.
#'
#' @inheritParams fit_hawkesNet
#' @param mu_vec Background rate at each event time (from
#'   \code{\link{prepare_inhomogeneous_background}}).
#' @param integral_bg Integral of background over time window.
#' @param ... Passed to \code{\link{fit_hawkesNet}} / \code{PMF_mark}
#'   (e.g. \code{formula_RHS}, \code{truncation}, \code{cores}).
#' @return Same as \code{\link{fit_hawkesNet}}.
#' @seealso \code{\link{fit_hawkesNet}}, \code{\link{prepare_inhomogeneous_background}}
#' @export
fit_hawkesNet_inhom <- function(params_init,
                                      time_window,
                                      mark_filtration,
                                      PMF_mark,
                                      mu_vec,
                                      integral_bg,
                                      maxit,
                                      trace = 0,
                                      reltol = 1e-8,
                                      parscale = NULL,
                                      fixed_params = NULL,
                                      cache_intensity = TRUE,
                                      combine_intensity = TRUE,
                                      method = "Nelder-Mead",
                                      verbose = TRUE,
                                      parallel_type = "auto",
                                      ...) {
  fit_hawkesNet(
    params_init = params_init,
    time_window = time_window,
    mark_filtration = mark_filtration,
    PMF_mark = PMF_mark,
    maxit = maxit,
    trace = trace,
    reltol = reltol,
    parscale = parscale,
    fixed_params = fixed_params,
    cache_intensity = cache_intensity,
    combine_intensity = combine_intensity,
    method = method,
    verbose = verbose,
    parallel_type = parallel_type,
    mu_vec = mu_vec,
    integral_bg = integral_bg,
    ...
  )
}


#' Build mu_vec and integral_bg from mark_filtration (network) using KDE
#'
#' Returns mu_vec aligned with get_times(mark_filtration)$times for use in
#' \code{\link{fit_hawkesNet}} (with \code{mu_vec} and \code{integral_bg} arguments).
#'
#' @param mark_filtration Network with vertex and edge time attributes
#' @param time_attr Name of the time attribute (default "time")
#' @param bw Bandwidth for KDE; NULL = default
#' @param grid_n Number of grid points for KDE
#' @return List with mu_vec, integral_bg, times, mu_fit, Lambda_fun
#' @export
prepare_inhomogeneous_background <- function(mark_filtration, time_attr = "time", bw = NULL, grid_n = 2048) {
  times_obj <- get_times(mark_filtration, time_name = time_attr)
  times <- times_obj$times
  times <- times[!is.na(times)]
  if (length(times) < 5) stop("Need more event times for KDE")

  windowT <- c(min(times), max(times))
  mu_fit <- estimate_mu_kde(times, windowT = windowT, bw = bw, grid_n = grid_n)
  Lambda_obj <- make_cumhaz_fun(mu_fit)
  mu_vec_at_events <- mu_fit$mu_fun(times)
  mu_vec_at_events <- pmax(mu_vec_at_events, 1e-12)
  integral_bg <- Lambda_obj$Lambda_fun(windowT[2])

  list(
    mu_vec = mu_vec_at_events,
    integral_bg = integral_bg,
    times = times,
    mu_fit = mu_fit,
    Lambda_fun = Lambda_obj$Lambda_fun
  )
}
