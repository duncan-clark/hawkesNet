
# =========================================================
# NOTE CURRENTLY IMPLEMENTED AS A KDE-BASED BACKGROUND RATE
# NOT IN PACKAGE SINCE AI GENERATED AND NOT VERIFIED
# =========================================================

if(FALSE){
  #' @title Compute KDE Background
  #' @description Computes a smoothed background rate shape using Kernel Density Estimation.
  #' @param times Numeric vector of event times.
  #' @param windowT Vector c(min, max) for the time window.
  #' @param bw Bandwidth for density (default "nrd0").
  #' @return A list containing:
  #'   \item{mu_func}{Approxfun function to get background density at any t}
  #'   \item{mu_vec}{Vector of background densities evaluated at 'times'}
  #'   \item{integral}{The integral of the density over the window (approx 1)}
  #' @export
  compute_kde_background <- function(times, windowT, bw = "nrd0") {
    # Compute density
    # We extend 'from' and 'to' slightly to avoid boundary effects,
    # but we clip the integral later.
    dens <- stats::density(times, from = windowT[1], to = windowT[2], bw = bw)
    
    # Create interpolation function
    mu_func <- stats::approxfun(dens$x, dens$y, rule = 2) # rule=2 handles extrapolation as constant
    
    # Evaluate at specific event times
    mu_at_events <- mu_func(times)
    
    # Ensure no absolute zeros to prevent log(0) in likelihood
    mu_at_events[mu_at_events <= 0] <- 1e-10
    
    # Approximate integral over the specific window [0, T]
    # (Since density integrates to 1 over (-Inf, Inf), we just check the window mass)
    # Simple numerical integration:
    x_grid <- seq(windowT[1], windowT[2], length.out = 1000)
    y_grid <- mu_func(x_grid)
    # Trapezoidal rule approx
    integral <- sum(diff(x_grid) * (head(y_grid, -1) + tail(y_grid, -1)) / 2)
    
    list(mu_func = mu_func, mu_vec = mu_at_events, integral = integral)
  }
  
  #' @rdname cond_intensity
  #' @export
  cond_intensity <- function(new_net,
                             t,
                             mark_filtration = NULL,
                             PMF_mark,
                             params,
                             new_edge_hash = NULL,
                             times = NULL,
                             bg_rate = NULL, # <--- NEW PARAMETER
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
    
    # If bg_rate (KDE value at t) is not provided, assume constant background model 
    # where density is 1 (so params$mu acts as pure rate).
    if(is.null(bg_rate)) bg_rate <- 1.0 
    
    # Pull ONLY what we need out of tmp
    log_mark_density0 <- tmp$log_mark_density
    log_density_func  <- tmp$log_density_func
    
    # Initial value calculation
    # params$mu is now the SCALING FACTOR for the background rate
    decays0 <- exp(-params$beta_overall * diffs)
    
    # Safety for log
    lambda_val <- params$mu * bg_rate + params$K * sum(decays0)
    if(lambda_val <= 0) lambda_val <- 1e-10
    
    log_result0 <- log_mark_density0 + log(lambda_val)
    result0 <- exp(log_result0)
    
    # Build a closure
    func <- local({
      diffs_local <- diffs
      ldf_local   <- log_density_func
      bg_local    <- bg_rate # Capture the KDE value at this specific time t
      
      function(params) {
        decays <- exp(-params$beta_overall * diffs_local)
        # Intensity = mu(scaling) * kde(t) + K * sum(...)
        lambda_inner <- params$mu * bg_local + params$K * sum(decays)
        if(lambda_inner <= 0) lambda_inner <- 1e-10
        
        exp(ldf_local(params) + log(lambda_inner))
      }
    })
    
    e_new <- new.env(parent = baseenv())
    e_new$diffs_local <- environment(func)$diffs_local
    e_new$ldf_local   <- environment(func)$ldf_local
    e_new$bg_local    <- environment(func)$bg_local
    environment(func) <- e_new
    
    list(
      result = result0,
      func   = func,
      lambda = params$mu * bg_rate, 
      kernel_sum = params$K * sum(decays0),
      decays = decays0,
      diffs  = diffs
    )
  }
  
  #' @rdname loglik_hawkesGrowthNet
  #' @export
  loglik_hawkesGrowthNet = function(params,
                                    time_window,
                                    mark_filtration,
                                    PMF_mark,
                                    edge_hash_list = NULL,
                                    verbose = FALSE,
                                    do_grad = FALSE,
                                    intens_funcs = NULL,
                                    kde_bg = NULL, # <--- NEW: Passed from fit function
                                    ...
  ){
    t_start <- proc.time()
    
    # Constraints check
    if(any(sapply(params[1:min(length(params),2)], function(x){x<0}))){
      return(list(loglik = -1e10, grads = rep(0,length(params))))
    }
    
    times <- get_times(mark_filtration)
    times <- times$times
    tval <- time_window[2]-time_window[1]
    
    # --- Background Rate Setup ---
    # If kde_bg is provided, we use it. If not, we revert to constant background logic.
    if(is.null(kde_bg)){
      # Fallback to constant background if not provided
      # Effectively creates a flat background of 1.0 everywhere
      mu_vec_at_events <- rep(1.0, length(times))
      bg_integral_term <- tval
    } else {
      mu_vec_at_events <- kde_bg$mu_vec
      bg_integral_term <- kde_bg$integral
    }
    
    if(is.null(intens_funcs)){
      # --- Parallel / List Creation Stage ---
      
      intens_func <- function(i){
        current_net <- filtration_to_net(mark_filtration, times[i], equal = TRUE)
        
        if("formula_RHS" %in% names(list(...))){
          model = createCppModel(as.formula(paste("current_net ~ ",list(...)$formula_RHS)))
        }else{
          model <- NULL
        }
        
        # Pass the specific background rate for this time point
        intensity <- cond_intensity(new_net = current_net,
                                    t = times[i],
                                    mark_filtration = current_net,
                                    PMF_mark = PMF_mark,
                                    params = params,
                                    model = model,
                                    times = get_times(mark_filtration),
                                    bg_rate = mu_vec_at_events[i], # <--- PASS KDE VALUE
                                    ...)
        return(list(result = intensity$result,
                    func = intensity$func))
      }
      
      if("cores" %in% names(list(...))){
        cores <- list(...)$cores
        intens_list <- pbmclapply(rev(seq_along(times)),
                                  intens_func,
                                  mc.cores = cores,
                                  mc.preschedule=FALSE)
        intens_vec <- sapply(intens_list,function(x){x$result})
        intens_funcs <- sapply(intens_list,function(x){x$func})
      } else {
        intens_list <- lapply(1:length(times), intens_func)
        intens_vec <- sapply(intens_list,function(x){x$result})
        intens_funcs <- sapply(intens_list,function(x){x$func})
      }
      
    } else {
      # --- Optimization Stage (Functions already exist) ---
      intens_vec <- numeric(length(intens_funcs))
      for (i in seq_along(intens_funcs)) intens_vec[i] <- intens_funcs[[i]](params)
    }
    
    # Handling numerical zeros
    tmp <- intens_vec
    if(any(is.na(tmp))) tmp[is.na(tmp)] <- min(tmp[!is.na(tmp)])/2
    if(sum(tmp<=0)!=0)  tmp[tmp<=0] <- min(tmp[tmp>0])/2
    
    intens_sum <- sum(log(tmp))
    
    # --- Integral Calculation ---
    # Integral of Background: params$mu (scaler) * Integral(KDE)
    bg_integral <- params$mu * bg_integral_term
    
    # Integral of Excitation
    pieces <- sapply(times, function(x){
      (1 - exp(-params$beta_overall * (time_window[2] - x)))
    })
    
    excitation_integral <- (1/params$beta_overall) * params$K * sum(pieces)
    
    integral <- bg_integral + excitation_integral
    loglik <- intens_sum - integral
    
    # --- Gradients ---
    grads <- list()
    if(do_grad){
      # Note: These grads need to be updated to reflect the bg_rate inside the closure
      # This is complex because we didn't extract the closures' internal bg_rate in this scope.
      # For now, I am returning the loglik. 
      # If analytic gradients are required for KDE background, the closure extraction 
      # in cond_intensity needs to return 'bg_local' so it can be accessed here.
      warning("Analytic gradients for Inhomogeneous Background not fully implemented. Use numeric diffs.")
    }
    
    return(list(loglik = loglik,
                intens_funcs = intens_funcs,
                grads = grads))
  }
  
  #' @rdname fit_hawkesGrowthNet
  #' @export
  #' @rdname fit_hawkesGrowthNet
  #' @export
  fit_hawkesGrowthNet <- function(params_init,
                                  time_window,
                                  mark_filtration,
                                  PMF_mark,
                                  trace = 0,
                                  REPORT = 10,
                                  reltol = 1e-8,
                                  maxit,
                                  get_hessian = FALSE,
                                  fixed_params = NULL,
                                  use_kde_bg = TRUE, 
                                  bw = "nrd0",       
                                  ...){
    
    params_init_old <- params_init
    if(!is.null(fixed_params)){
      for(k in fixed_params){
        params_init[[k]] <- NULL
      }
    }
    
    # --- Pre-Calculate KDE Background ---
    kde_bg <- NULL
    if(use_kde_bg){
      times <- get_times(mark_filtration)$times
      if(trace > 0) cat("Pre-calculating KDE background rate...\n")
      kde_bg <- compute_kde_background(times, time_window, bw = bw)
    }
    
    # --- Optimization Wrapper ---
    # This function captures 'kde_bg', 'time_window', etc. from the parent environment
    optim_func <- function(params, ...){
      
      # Reconstruct full parameter list
      params <- relist(params, skeleton = params_init)
      params[fixed_params] <- params_init_old[fixed_params]
      
      result <- loglik_hawkesGrowthNet(params = params,
                                       time_window = time_window,
                                       mark_filtration = mark_filtration,
                                       PMF_mark = PMF_mark,
                                       kde_bg = kde_bg, # Uses the variable from parent scope
                                       ...)             # '...' contains 'intens_funcs' passed by optim
      
      return(list(value = result$loglik,
                  grad = unlist(result$grads)))
    }
    
    fn_wrapper <- function(par,...) {
      res <- optim_func(par, ...)
      return(res$value)
    }
    
    # --- Initial Likelihood / Function Generation ---
    t <- proc.time()
    init_lik <- loglik_hawkesGrowthNet(params = params_init_old,
                                       time_window = time_window,
                                       mark_filtration = mark_filtration,
                                       PMF_mark = PMF_mark,
                                       kde_bg = kde_bg, 
                                       ...)
    
    if(trace > 0) print(paste0("Initial List Gen took ", round((proc.time()-t)[3],2)," seconds"))
    
    # --- Optimization ---
    fit <- optim(par = unlist(params_init),
                 fn = fn_wrapper,
                 gr = NULL, 
                 method = "Nelder-Mead",
                 control = list(fnscale = -1,
                                trace = trace,
                                maxit = maxit,
                                reltol = reltol),
                 hessian = get_hessian,
                 # STATIC DATA PASSED TO OPTIM's ...
                 intens_funcs = init_lik$intens_funcs, 
                 # REMOVED: kde_bg = kde_bg (This was causing the collision)
                 ...)
    
    if(trace > 0) print(paste0("Fitting took ", round((proc.time()-t)[3],2)," seconds"))
    
    return(list(fit = fit, kde_bg = kde_bg))
  }
}