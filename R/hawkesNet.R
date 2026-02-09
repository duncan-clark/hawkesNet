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
    # Debug: print what we found
    if (length(cs_indices) > 0) {
      message("DEBUG: Found CS_params indices: ", paste(cs_indices, collapse=", "))
      message("DEBUG: Found CS_params names: ", paste(fit_table$parameter[cs_indices], collapse=", "))
    }
    message("DEBUG: Expected ", n_cs, " statistics: ", paste(stat_names, collapse=", "))
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
    # --- START DEBUG ---
    # 1. Get the environment where variables like 'diffs_local' should live
    #    (This is the parent of the current execution environment)
    enclosure <- parent.env(environment())
    # --- END DEBUG ---
    # Explicitly using the variables expected in e_tiny
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

# According to CONOR we can't have a branching process due to the background rate issues:
# Now we want the to do a thinning approach, i.e. propose a bunch of points and the accept or reject them
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
                                hashed_edges = F,
                                verbose = F,
                                mu_multiplier = 10,
                                joint_accept = F,
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

  # Initialize the output list of events
  events = list()
  events$n = 0
  events$t = c()
  events$mark_density <- c()

  accept_probs <- c()

  # propose points to be thinned:
  n_bg = rpois(1, lambda * (time_window[2] - time_window[1]))
  event_queue <- data.table(time = sort(runif(n_bg, min=0, max=time_window[2])))
  # maintain order so no need to sort
  setkey(event_queue, time)
  # Initialize the list to store new events
  new_events_list <- list()
  list_index <- 1
  tot <- 0
  tot_attempt <- 0
  current_net <- network::network(matrix(1),directed = F)
  delete.vertices(current_net,1)
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
      # debug(PMF_mark)
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
      new_nodes <- (net %n% 'n') - (current_net %n% 'n')
      if(hashed_edges && length(net$mel)!=0){
        # hash the network edge list for fast lookup:
        edges <- network::as.edgelist(net)
        keys_vec <- paste(edges[,1], edges[,2], sep = "-")
        edge_hash <- hash::hash(keys = keys_vec, values = rep(TRUE, length(keys_vec)))
      }else{
        edge_hash <- NULL
      }
      # Get mu_at_t for inhomogeneous case
      mu_at_t <- if (use_inhom) {
        mu_val <- mu_fun(current_event$time)
        pmax(mu_val, 1e-12)  # Ensure positive
      } else {
        NULL
      }
      
      if(joint_accept){
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
            new_nodes <- (net %n% 'n') - (current_net %n% 'n')
            if(hashed_edges){
              # hash the network edge list for fast lookup:
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
          if (use_inhom) {
            tmp <- cond_intensity_inhom(new_net = net,
                                       t = current_event$time,
                                       mark_filtration = current_net,
                                       PMF_mark = PMF_mark,
                                       params = params,
                                       mu_at_t = mu_at_t,
                                       new_edge_hash = edge_hash,
                                       ...
            )
            intensity <- tmp$lambda + tmp$kernel_sum
          } else {
            tmp <- cond_intensity(new_net = net,
                                 t = current_event$time,
                                 mark_filtration = current_net,
                                 PMF_mark = PMF_mark,
                                 params = params,
                                 new_edge_hash = edge_hash,
                                 ...
            )
            intensity <- tmp$lambda + tmp$kernel_sum
          }
        }
        # Use the same proposed mark for acceptance and for updating (do not resample)
      }

      accept <- intensity/lambda
    }
    accept_probs <- c(accept_probs,accept)

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
      events$t[length(events$t)+1] <- current_event$time
      if(length(events$t) >2){
        events$mark_density <- c(events$mark_density,mark_sample$mark_density)
      }

    }else{
      # do nothing since we rejected the point
    }
    if(verbose){
      print(paste0("time is ",current_event$time, " size of net is ",current_net %n% 'n',' number of edges is ',length(current_net$mel)))
      print(paste0("time is ",current_event$time, " this iteration of while loop took ", round((proc.time()-t)[3],2)," seconds"))
    }
    # Concatenate new events to event_queue only if we have only one event left to go
    old_n <- dim(event_queue)[1]
    event_queue <- rbindlist(list(event_queue, rbindlist(new_events_list, use.names = TRUE)))
    new_events_list <- vector("list", length(new_events_list))  # Reset the list
    list_index <- 1
  }
  t1 <- proc.time() - t1
  print(paste0("simulation took ",round(t1[3],2)," seconds"))
  return(list(events = events,
              net = current_net,
              accept_probs = accept_probs))
}

#' Log-likelihood for the Hawkes network growth model
#'
#' Computes the log-likelihood for a given parameter vector and observed network.
#'
#' @param params List of parameters.
#' @param time_window Numeric \code{c(t0, t1)}.
#' @param mark_filtration Observed network (filtration).
#' @param PMF_mark Mark PMF function.
#' @param edge_hash_list Optional list of edge hashes per event (for internal use).
#' @param verbose Print timing (default \code{FALSE}).
#' @param intens_funcs Precomputed intensity functions (for fitting).
#' @param ... Passed to \code{cond_intensity} / \code{PMF_mark} (e.g. \code{truncation}, \code{formula_RHS}).
#' @return List with \code{loglik} and \code{intens_funcs}.
#' @seealso \code{\link[network]{as.edgelist}}, \code{\link[hash]{hash}}
#' @rdname loglik_hawkesNet
#' @export
loglik_hawkesNet = function(params,
                                  time_window,
                                  mark_filtration,
                                  PMF_mark,
                                  edge_hash_list = NULL,
                                  verbose = FALSE,
                                  intens_funcs = NULL,
                                  ...
){
  t<-proc.time()
  # don't allow negative parameters in first 2
  if(any(sapply(params[1:min(length(params),4)],function(x){x<0}))){
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
  
  if(is.null(intens_funcs)){
    # do the sum of the intensities:
    intens_sum <- 0
    times_precalc <- get_times(mark_filtration)
    
    # if parallelize do that here with PSOCK for simplicity:
    intens_func <- function(i){
      # need to do this on the fly otherwise too storage intensive
      current_net <- filtration_to_net(mark_filtration,times[i],equal = TRUE)
      if("formula_RHS" %in% names(list(...))){
        model = createCppModel(as.formula(paste("current_net ~ ",list(...)$formula_RHS)))
      }else{
        model <- NULL
      }
      intensity <- cond_intensity(new_net = current_net,
                                  t = times[i],
                                  mark_filtration = current_net,
                                  PMF_mark = PMF_mark,
                                  params = params,
                                  model = model,
                                  times = times_precalc,
                                  ...)
      return(list(result = intensity$result,
                  func = intensity$func
                  ))
    }
    cores <- list(...)$cores
    if(!is.null(cores) && is.numeric(cores) && cores > 1){
      if(verbose){
        print(paste0("using ",cores," cores on ",length(times), " objects for cond intensity list first calculation"))
      }
      print("starting intens list calculation")
      t <- proc.time()
      intens_list <- pbmcapply::pbmclapply(rev(seq_along(times)),
                                intens_func,
                                mc.cores = cores,
                                mc.preschedule=FALSE)
      print(paste0("intens list ", round((proc.time()-t)[3],2)," seconds"))
      intens_vec <- sapply(intens_list,function(x){x$result})
      intens_funcs <- sapply(intens_list,function(x){x$func})
    }else{
      t1 <- proc.time()
      intens_list <- lapply(seq_along(times),function(i){
        intens_func(i)})
      intens_vec <- sapply(intens_list,function(x){x$result})
      intens_funcs <- sapply(intens_list,function(x){x$func})
      print(paste0("intens list took ", round((proc.time()-t1)[3],2)," seconds"))
    }
  }else{
    t1 <- proc.time()
    intens_vec <- numeric(length(intens_funcs))
    for (i in seq_along(intens_funcs)) intens_vec[i] <- intens_funcs[[i]](params)
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

  # Integral due to kernel being density:

  max_t <- max(times)
  pieces <- 1 - exp(-params$beta_overall * (tval - times))
  # Guard against division by zero or very small beta_overall
  if (!is.finite(params$beta_overall) || params$beta_overall <= 0 || params$beta_overall < 1e-10) {
    return(list(loglik = -1e10, intens_funcs = intens_funcs))
  }
  integral <- params$mu * tval + (1/params$beta_overall)*params$K*sum(pieces)
  loglik <- intens_sum - integral  # If you hit Browse[] here, clear RStudio breakpoints (Debug -> Clear All) or run undebug(loglik_hawkesNet)
  
  # Guard: ensure loglik is finite for L-BFGS-B
  if (!is.finite(loglik) || !is.finite(intens_sum) || !is.finite(integral)) {
    if (verbose) {
      warning("loglik_hawkesNet: non-finite loglik detected. intens_sum=", intens_sum, 
              ", integral=", integral, ", returning -1e10")
    }
    return(list(loglik = -1e10, intens_funcs = intens_funcs))
  }
  
  # print(paste0("integral is ",integral))
  # print(paste0("intens_sum is ",intens_sum))
  # print(paste0("trigger part of integral is  ",(1/params$beta_overall)*params$K*sum(pieces)))
  t<-proc.time() - t
  if(verbose){
    print(paste0("this iteration of loglik took ", round(t[3],2)," seconds"))
  }
  
  return(list(loglik = loglik,
              intens_funcs = intens_funcs))
}

#' Fit the Hawkes network growth model by maximum likelihood
#'
#' @param params_init List of initial parameter values.
#' @param time_window Numeric \code{c(t0, t1)}.
#' @param mark_filtration Observed network.
#' @param PMF_mark Mark PMF function.
#' @param maxit Maximum number of iterations for the optimizer.
#' @param trace Trace level (default 0).
#' @param REPORT Reporting interval for \code{optim} (default 10).
#' @param reltol Relative convergence tolerance (default 1e-8).
#' @param parscale Scale vector for parameters (default all 1).
#' @param fixed_params Character vector of parameter names to hold fixed.
#' @param cache_intensity If \code{TRUE}, cache intensity functions (default \code{TRUE}).
#' @param method Optimization method: \code{"Nelder-Mead"} (default, derivative-free) or
#'   \code{"L-BFGS-B"} (gradient-based with box constraints; constrains mu, beta_overall,
#'   beta_edges, K, node_lambda > 0 automatically).
#' @param ... Passed to \code{loglik_hawkesNet} (e.g. \code{truncation}, \code{formula_RHS}).
#' @return List with \code{fit} (output of \code{optim}), \code{intens_funcs}, \code{fit_table} (parameter estimates and standard errors), and \code{hessian} (numerical Hessian of negative log-likelihood at MLE, if numDeriv available).
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
                                method = "Nelder-Mead",
                                ...){
  params_init_old <- params_init
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
  
  t<-proc.time()
  
  # pre-calculate the param -> conditonal intensity mapping
  # since the observation never changes - not need to do expensive network processes every iteration
  # then param -> likelihood should be very fast
  # 1. Pre-calculate ONLY if cache_intensity is TRUE
  if(cache_intensity){
    print("Pre-calculating intensity closures (Fast Mode)...")
    init_lik <- loglik_hawkesNet(params = params_init_old,
                                     time_window = time_window,
                                     mark_filtration = mark_filtration,
                                     PMF_mark = PMF_mark,
                                     ...)
    cached_funcs <- init_lik$intens_funcs
  } else {
    print("Caching disabled (Safe Mode) ...")
    cached_funcs <- NULL
    init_lik <- NULL
  }
  dot_args <- list(...)
  optim_func <- function(params){
    params_curr <- relist(params, skeleton = params_init)
    # Restore vertex_categorical_levels from the original params (not optimized)
    params_curr$vertex_categorical_levels <- params_init_old$vertex_categorical_levels
    # Re-inject the fixed parameters from the backup (params_init_old)
    if (!is.null(fixed_params)) {
      for (k in fixed_params) {
        params_curr[[k]] <- params_init_old[[k]]
      }
    }
    if (!point_process_params_valid(params_curr)) return(-1e10)
    
    # 2. Pass NULL to intens_funcs if caching is disabled
    # This forces loglik_hawkesNet to rebuild the density from scratch
    result <- tryCatch({
      do.call(loglik_hawkesNet, c(
        list(params = params_curr,
             time_window = time_window,
             mark_filtration = mark_filtration,
             PMF_mark = PMF_mark,
             intens_funcs = cached_funcs),
        dot_args
      ))
    }, error = function(e) {
      return(list(loglik = -1e10, intens_funcs = cached_funcs))
    })
    ll <- result$loglik
    # Final guard: ensure return value is finite for L-BFGS-B
    if (!is.finite(ll)) return(-1e10)
    return(ll)
  }
  flat_par <- unlist(params_init)
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
  fit <- do.call(optim, optim_args)
  print(paste0("fitting took ",round((proc.time()-t)[3],2)," seconds"))

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
  # optim returns hessian of fn (loglik), which is negative definite at a maximum.
  # The variance-covariance matrix is the inverse of the *negative* hessian (i.e. observed information).
  vcov <- tryCatch(solve(-hessian), error = function(e) NULL)
  if (!is.null(vcov)) {
      se <- sqrt(pmax(diag(vcov), 0))
      fit_table$std.error <- se
  }
  # Replace CS_params1, CS_params2, ... with actual ERNM statistic names
  fit_table <- rename_CS_params_in_table(fit_table, mark_filtration, list(...))
  message("Hawkes growth fit results:")
  message("Total parameters in fit table: ", nrow(fit_table))
  # Print all rows explicitly (max = NULL means print all rows)
  print(fit_table, max = NULL, row.names = TRUE)
  if (all(is.na(fit_table$std.error))) {
    message("(Standard errors not available; install numDeriv for SEs.)")
  }
  
  return(list(fit = fit,
              intens_funcs = if (cache_intensity) init_lik$intens_funcs else cached_funcs,
              fit_table = fit_table,
              hessian = hessian))
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
    integral <- params$mu * times[i] + (1/params$beta_overall)*params$K*sum(pieces[1:i])
    integral <- params$mu * times[i] + params$K*sum(pieces[1:i])
    return(integral)
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


