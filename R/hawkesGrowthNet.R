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
#' @return List with \code{result} (intensity value), \code{func} (function for gradient), and optional debug fields.
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
#'
#' @param params List of parameters (\code{mu}, \code{beta_overall}, \code{K}, \code{beta_edges}, and any mark-specific).
#' @param time_window Numeric vector \code{c(t0, t1)}.
#' @param PMF_mark Mark PMF function (e.g. \code{PMF_mark_BA} or \code{PMF_mark_CS}).
#' @param cond_intensity Conditional intensity function (e.g. \code{cond_intensity}).
#' @param hashed_edges If \code{TRUE}, use a hash for edge lookup (default \code{FALSE}).
#' @param verbose Print progress (default \code{FALSE}).
#' @param mu_multiplier Multiplier for thinning upper bound (default 10).
#' @param joint_accept Logical (default \code{FALSE}); joint acceptance for mark and time.
#' @param n_mark_sample Optional number of mark samples per proposal.
#' @param ... Passed to \code{PMF_mark} or \code{cond_intensity} (e.g. \code{truncation}, \code{formula_RHS}).
#' @return List with \code{events}, \code{net}, \code{accept_probs}.
#' @seealso \code{\link[network]{as.edgelist}}, \code{\link[hash]{hash}}
#' @rdname sim_hawkesGrowthNet
#' @export
#' @importFrom network as.edgelist
#' @importFrom hash hash
sim_hawkesGrowthNet <- function(params,
                                time_window,
                                PMF_mark, # function that both generates new mark and calculates the density of existing mark
                                cond_intensity, # function to calcualte condiational_intensity, takes in a kernel_func
                                hashed_edges = F,
                                verbose = F,
                                mu_multiplier = 10,
                                joint_accept = F,
                                n_mark_sample = NULL,
                                ... # to be past to PMF_mark


){
  t1 <- proc.time()
  # simulate the background points (can only simulate their times right now)
  mu <- params$mu
  theta <- params$theta
  beta <- params$beta
  K <- params$K

  # poisson in time lambda
  lambda <- mu_multiplier*mu

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
      if(joint_accept){
        intensity <- cond_intensity(new_net = net,
                                    t = current_event$time,
                                    mark_filtration = current_net,
                                    PMF_mark = PMF_mark,
                                    params = params,
                                    new_edge_hash = edge_hash,
                                    ...
        )$result
      }else{
        if(!is.null(n_mark_sample)){
          imp_sample <- sapply(1:n_mark_sample,function(i){
            mark_sample <- PMF_mark(time = current_event$time,
                                    params = params,
                                    mark_filtration = current_net,
                                    mark = NULL,
                                    generate_mark = TRUE,
                                    new_edge_hash = TRUE,
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
            intensity <- cond_intensity(new_net = net,
                                        t = current_event$time,
                                        mark_filtration = current_net,
                                        PMF_mark = PMF_mark,
                                        params = params,
                                        new_edge_hash = edge_hash,
                                        ...
            )
            return(intensity$result/mark_sample$mark_sample_density)
          })
          intensity <- mean(imp_sample)
        }else{
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

        mark_sample <- PMF_mark(time = current_event$time,
                                params = params,
                                mark_filtration = current_net,
                                mark = NULL,
                                generate_mark = TRUE,
                                new_edge_hash = NULL,
                                ...
        )
        net <- mark_sample$mark_sample
      }

      accept <- intensity/lambda
    }
    accept_probs <- c(accept_probs,accept)

    # if we accept the point add it in
    if(verbose){
      print(paste0("Number of edges proposed is ",length(net$mel)))
      print(paste0("Number of nodes proposed is ",length(net$oel)))
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
      print(paste0("time is ",current_event$t, " size of net is ",current_net %n% 'n',' number of edges is ',length(current_net$mel)))
      print(paste0("time is ",current_event$t, " this iteration of while loop took ", round((proc.time()-t)[3],2)," seconds"))
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
#' Computes the log-likelihood and optional gradients for a given parameter vector and observed network.
#'
#' @param params List of parameters.
#' @param time_window Numeric \code{c(t0, t1)}.
#' @param mark_filtration Observed network (filtration).
#' @param PMF_mark Mark PMF function.
#' @param edge_hash_list Optional list of edge hashes per event (for internal use).
#' @param verbose Print timing (default \code{FALSE}).
#' @param do_grad Compute gradients (default \code{FALSE}).
#' @param intens_funcs Precomputed intensity functions (for fitting).
#' @param ... Passed to \code{cond_intensity} / \code{PMF_mark} (e.g. \code{truncation}, \code{formula_RHS}).
#' @return List with \code{loglik} and optionally \code{grads}, \code{intens_funcs}.
#' @seealso \code{\link[network]{as.edgelist}}, \code{\link[hash]{hash}}
#' @rdname loglik_hawkesGrowthNet
#' @export
#' @importFrom network as.edgelist
#' @importFrom hash hash
loglik_hawkesGrowthNet = function(params,
                                  time_window,
                                  mark_filtration,
                                  PMF_mark,
                                  edge_hash_list = NULL,
                                  verbose = FALSE,
                                  do_grad = FALSE,
                                  intens_funcs = NULL,
                                  ...
){
  t<-proc.time()
  # don't allow negative parameters in first 2
  if(any(sapply(params[1:min(length(params),4)],function(x){x<0}))){
    return(list(loglik = -(10**(100)),
                grads = rep(0,length(params)))
    )
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
    if("cores" %in% names(list(...))){
      if(verbose){
        print(paste0("using ",list(...)$cores," cores on ",length(times), " objects for cond intensity list first calculation"))
      }
      print("starting intens list calculation")
      cores <- list(...)$cores
      t <- proc.time()
      intens_list <- pbmclapply(rev(seq_along(times)),
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
  integral <- params$mu * tval + (1/params$beta_overall)*params$K*sum(pieces)
  loglik <- intens_sum - integral
  # print(paste0("integral is ",integral))
  # print(paste0("intens_sum is ",intens_sum))
  # print(paste0("trigger part of integral is  ",(1/params$beta_overall)*params$K*sum(pieces)))
  # print(paste0("result is :",loglik))
  # print(paste0("this iteration of loglik took ", round((proc.time()-t)[3],2)," seconds"))

  # calcualte the gradients:

  # kernel_sum <- sapply(intens_list,function(x){x$kernel_sum})
  # decays <- lapply(intens_list,function(x){x$decays})
  # diffs <- lapply(intens_list,function(x){x$diffs})
  
  grads <- list()
  if(do_grad){
    # hawkes_mu
    grads$mu <- sum(1/(params$mu + params$K*kernel_sum)) - tval
    # hawkes K
    grads$K <-  sum(kernel_sum/(params$mu + params$K*kernel_sum)) -  (1/params$beta_overall)*sum(pieces)
    # hawkes beta
    grads$beta_overall <- sum(sapply(seq_along(decays), function(i){
      (params$K * sum(-diffs[[i]]*decays[[i]])) / (params$mu + kernel_sum[i])
    })) +
      (-params$K)/(params$beta_overall**2) * ( length(kernel_sum) + sum((tval - times - 1/(params$beta_overall^2))*exp(-params$beta_overall*(tval-times))))
    
    # hawkes edge decay from mark generator:
    decay_grads <- sapply(intens_list,function(x){x$decay_grad})
    grads$beta_edges <- sum(decay_grads)
    
    # hawkes theta_params (from mark generator)
    # get the mark PMF grads:
    mark_grads <- lapply(intens_list,function(x){x$mark_grad})
    grads$CS_params <- colSums(do.call(rbind,mark_grads))
  }


  # TODO fix the mark grad issues with simplification

  t<-proc.time() - t
  if(verbose){
    print(paste0("this iteration of loglik took ", round(t[3],2)," seconds"))
  }
  
  return(list(loglik = loglik,
              intens_funcs = intens_funcs,
              grads = grads))
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
#' @param get_hessian If \code{TRUE}, return Hessian in the fit object.
#' @param fixed_params Character vector of parameter names to hold fixed.
#' @param cache_intensity If \code{TRUE}, cache intensity functions (default \code{TRUE}).
#' @param ... Passed to \code{loglik_hawkesGrowthNet} (e.g. \code{truncation}, \code{formula_RHS}).
#' @return List with \code{fit} (output of \code{optim}) and \code{intens_funcs}.
#' @rdname fit_hawkesGrowthNet
#' @export
fit_hawkesGrowthNet <- function(params_init,
                                time_window,
                                mark_filtration,
                                PMF_mark,
                                maxit,
                                trace = 0,
                                REPORT = 10,
                                reltol = 1e-8,
                                parscale = NULL,
                                get_hessian = FALSE,
                                fixed_params = NULL,
                                cache_intensity = TRUE,
                                ...){
  params_init_old <- params_init
  if(!is.null(fixed_params)){
    for(k in fixed_params){
      params_init[[k]] <- NULL
    }
  }
  
  if (is.null(parscale)) {
    flat_params <- unlist(params_init)
    parscale <- rep(1, length(flat_params))
  }
  
  # ==================
  # Deprecated 
  # ==================
  
  # optim_func <- function(params,...){
  #   param_vec <- params
  #   params <- relist(params,skeleton = params_init)
  #   params[fixed_params] <- params_init_old[fixed_params]
  #   
  #   result <- loglik_hawkesGrowthNet(params = params,
  #                                    time_window = time_window,
  #                                    mark_filtration = mark_filtration,
  #                                    PMF_mark = PMF_mark,
  #                                    ...)
  # 
  #   # =================
  #   # USE numDERIV TO DEBUG !
  #   # =================
  #   wrapper <- function(x,...){
  #     x <- relist(x, skeleton = params_init)
  #     return(loglik_hawkesGrowthNet(params = x,
  #                                   ...)$loglik)
  #   }
  # 
  #   # numgrad <- grad(func = wrapper,
  #   #                 x    = param_vec,
  #   #                 time_window = time_window,  # or whatever your data is
  #   #                 events      = events,
  #   #                 PMF_mark    = PMF_mark,
  #   #                 ...)
  #   #
  #   # print("Numeric gradient vs analytic gradient:")
  #   # print(numgrad)
  #   #print("analytic gradient")
  #   #print(unlist(result$grads))
  # 
  #   # print("params are:")
  #   # print(params)
  #   # print("Loglik is :")
  #   # print(result$loglik)
  #   # print("Grads are:")
  #   # print(result$grads)
  #   return(list(value = result$loglik,
  #               #grad = numgrad
  #               grad = unlist(result$grads)
  #   ))
  # }
  # 
  # fn_wrapper <- function(par,...) {
  #   res <- optim_func(par, ...)
  #   return(res$value)
  # }
  # gr_wrapper <- function(par, ...) {
  #   res <- optim_func(par, ...)
  # 
  #   return(res$grad)
  # }
  t<-proc.time()
  
  # pre-calculate the param -> conditonal intensity mapping
  # since the observation never changes - not need to do expensive network processes every iteration
  # then param -> likelihood should be very fast
  # 1. Pre-calculate ONLY if cache_intensity is TRUE
  if(cache_intensity){
    print("Pre-calculating intensity closures (Fast Mode)...")
    init_lik <- loglik_hawkesGrowthNet(params = params_init_old,
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
  optim_func <- function(params, ...){
    params_curr <- relist(params, skeleton = params_init)
    # 2. Re-inject the fixed parameters from the backup (params_init_old)
    #    to reconstruct the full parameter list required by the loglik function.
    if (!is.null(fixed_params)) {
      for (k in fixed_params) {
        params_curr[[k]] <- params_init_old[[k]]
      }
    }
    
    # 2. Pass NULL to intens_funcs if caching is disabled
    # This forces loglik_hawkesGrowthNet to rebuild the density from scratch
    result <- loglik_hawkesGrowthNet(params = params_curr,
                                     time_window = time_window,
                                     mark_filtration = mark_filtration,
                                     PMF_mark = PMF_mark,
                                     intens_funcs = cached_funcs, # Pass NULL if disabled
                                     ...)
    return(result$loglik)
  }
  fit <- optim(par = unlist(params_init),
               fn = optim_func,
               method = "Nelder-Mead",
               control = list(fnscale = -1,
                              trace=trace,
                              maxit=maxit,
                              reltol = reltol,
                              parscale = parscale,
                              abstol = NULL),
               hessian = get_hessian,
               ...)
  print(paste0("fitting took ",round((proc.time()-t)[3],2)," seconds"))
  
  # Numerically estimate Hessian at optimal params
  # THIS IS GONNA TAKE FOREVER - NEED TO CODE UP GRADIENT!
  # hessian_estimate <- numDeriv::hessian(
  #   func = function(p) {
  #     cat(sprintf("Parameters: %s\n",paste(round(p, 4), collapse = ", ")))
  #     p <- relist(p, skeleton = params_init)
  #     result <- -loglik_hawkesGrowthNet(params = p,
  #                                       time_window = time_window,
  #                                       mark_filtration = mark_filtration,
  #                                       PMF_mark = PMF_mark,
  #                                       ...)$loglik
  #     cat(sprintf("  --> loglik: %f\n", -result))
  #     
  #     
  #     return(result)
  #   },
  #   x = fit$par
  # )
  # 
  # # Fisher information approximation is the negative Hessian at the optimum
  # fisher_info <- hessian_estimate
  # vcov_matrix <- solve(fisher_info)
  

  
  return(list(fit=fit,
              intens_funcs = init_lik$intens_funcs))
}


#' Compensators (integrated intensity) for the Hawkes growth model
#'
#' @param params List of parameters.
#' @param time_window Numeric \code{c(t0, t1)}.
#' @param mark_filtration Observed network.
#' @return Numeric vector of compensator values at each event time.
#' @export
compensators_hawkesGrowthNet <- function(params,
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
ks_test_pval_hawkesGrowthNet <- function(params,
                                         time_window,
                                         mark_filtration){
  compensators <- compensators_hawkesGrowthNet(params = unlist(params),
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


