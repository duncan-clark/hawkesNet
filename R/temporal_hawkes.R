## ============================================================
## Power-law Hawkes (temporal) fit — KDE Background
##   - Background: mu(t) = gamma * KDE(t)
##   - Kernel: g(t) = (p-1) c^(p-1) / (t + c)^p 
##   - Optimization: L-BFGS-B with native constraints
## ============================================================

## ---- Helper: Pre-compute KDE Background ----
# Call this before fitting to get the background 'mu' vector
compute_kde_background <- function(times, windowT, bw = "nrd0") {
  dens <- stats::density(times, from = windowT[1], to = windowT[2], bw = bw)
  # Interpolate to get exact density at event times
  mu_at_events <- stats::approx(dens$x, dens$y, xout = times)$y
  # Return both the values at events and the total integral (should be ~1)
  list(mu_vec = mu_at_events, total_int = 1) 
}

hawkes_kernel_density <- function(dt, kernel = c("exp", "powerlaw"), beta = NULL, c = NULL) {
  kernel <- match.arg(kernel)
  out <- matrix(0, nrow = nrow(dt), ncol = ncol(dt))
  if (kernel == "exp") {
    if (is.null(beta) || beta <= 0) return(out)
    out <- beta * exp(-beta * dt)
  } else {
    p <- beta
    if (is.null(p) || is.null(c) || p <= 1 || c <= 0) return(out)
    out <- (p - 1) * (c^(p - 1)) / ((dt + c)^p)
  }
  out[dt < 0] <- 0
  out
}

hawkes_kernel_cdf <- function(u, kernel = c("exp", "powerlaw"), beta = NULL, c = NULL) {
  kernel <- match.arg(kernel)
  u[u < 0] <- 0
  if (kernel == "exp") {
    if (is.null(beta) || beta <= 0) return(rep(0, length(u)))
    return(1 - exp(-beta * u))
  } else {
    p <- beta
    if (is.null(p) || is.null(c) || p <= 1 || c <= 0) return(rep(0, length(u)))
    return(1 - (c / (u + c))^(p - 1))
  }
}

## ---- Updated Log-likelihood with KDE Background ----
loglik_temporal_hawk <- function(params, realiz, windowT, dists = NULL,
                                 kde_bg = NULL, kernel = c("exp", "powerlaw")) {
  kernel <- match.arg(kernel)
  
  # Params: gamma (KDE scale), beta/p (decay), K (branching), [c (powerlaw)]
  gamma <- params[1]
  beta  <- params[2]
  K     <- params[3]
  c_pl  <- if (kernel == "powerlaw") params[4] else NULL
  
  # Native constraints for Nelder-Mead (though we'll use L-BFGS-B)
  if (gamma < 0 || K < 0 || K >= 1) return(-1e10)
  
  # Intensity at event times: lambda(tj) = gamma*KDE(tj) + K * sum g(tj-ti)
  # kde_bg$mu_vec must be pre-calculated for these exact event times
  mu_term <- gamma * kde_bg$mu_vec
  
  if (is.null(dists)) {
    time_dist <- outer(realiz$t, realiz$t, "-")
  } else {
    time_dist <- dists$time_dist
  }
  
  kern_mat <- hawkes_kernel_density(dt = time_dist, kernel = kernel, beta = beta, c = c_pl)
  kern_mat[upper.tri(kern_mat, diag = TRUE)] <- 0
  
  g_sum <- rowSums(kern_mat)
  lambdas <- mu_term + K * g_sum
  
  if (any(lambdas <= 0) || any(is.na(lambdas))) return(-1e10)
  
  # Integral term: gamma * ∫KDE + K * sum G(T-ti)
  # Since KDE integrates to 1 over the window:
  u <- windowT[2] - realiz$t
  G_vals <- hawkes_kernel_cdf(u, kernel = kernel, beta = beta, c = c_pl)
  
  int_lam <- gamma * kde_bg$total_int + K * sum(G_vals)
  
  loglik <- sum(log(lambdas)) - int_lam
  return(loglik)
}

fit_temporal_hawkes <- function(params_init,
                                realiz,
                                windowT,
                                method = "L-BFGS-B", 
                                maxit = 200,
                                kernel = c("exp", "powerlaw"),
                                trace = 0,
                                low = NULL,
                                upp = NULL) {
  kernel <- match.arg(kernel)
  realiz <- realiz[order(realiz$t),,drop = FALSE]
  
  # Pre-calculate KDE background rate
  bg <- compute_kde_background(realiz$t, windowT)
  dists <- list(time_dist = outer(realiz$t, realiz$t, "-"))
  
  # --- 1. Define Default Constraints ---
  # We define them regardless, but we only pass them to optim if using L-BFGS-B
  default_low <- if(kernel == "exp") c(1e-5, 1e-5, 1e-5) else c(1e-5, 1.001, 1e-5, 1e-5)
  default_upp <- if(kernel == "exp") c(Inf, Inf, 0.9999) else c(Inf, Inf, 0.9999, Inf)
  
  # Use user-provided bounds if they exist, otherwise use defaults
  if(is.null(low)) low <- default_low
  if(is.null(upp)) upp <- default_upp
  
  # --- 2. Handle Optimization Method ---
  if (method == "Nelder-Mead") {
    # Nelder-Mead ignores 'lower'/'upper' in optim().
    # We must pass NULL to avoid warnings, and rely on the internal 
    # 'return(-1e10)' check in the log-likelihood function to handle bounds.
    optim_lower <- -Inf
    optim_upper <- Inf
  } else {
    # L-BFGS-B uses these bounds strictly
    optim_lower <- low
    optim_upper <- upp
  }
  
  optim(
    par = unlist(params_init),
    fn = loglik_temporal_hawk,
    method = method,
    lower = optim_lower, # Pass NULL/-Inf if Nelder-Mead
    upper = optim_upper, # Pass NULL/Inf if Nelder-Mead
    control = list(fnscale = -1, trace = trace, maxit = maxit),
    realiz = realiz,
    windowT = windowT,
    dists = dists,
    kde_bg = bg,
    kernel = kernel,
    hessian = TRUE
  )
}

## ---- Updated Compensator ----
compensator_temporal_hawkes <- function(params,
                                        realiz,
                                        windowT,
                                        kernel = c("exp", "powerlaw")
                                        ) {
  kernel <- match.arg(kernel)
  realiz <- realiz[order(realiz$t),,drop = FALSE]
  bg <- compute_kde_background(realiz$t, windowT)
  
  gamma <- params[1]; beta <- params[2]; K <- params[3]
  c_pl <- if (kernel == "powerlaw") params[4] else NULL
  
  # KDE integral up to each ti (Cumulative Density Function of the KDE)
  # We use the empirical CDF of the background here
  dens <- stats::density(realiz$t, from = windowT[1], to = windowT[2])
  # Integrate the KDE density numerically for the background part
  bg_cdf_fun <- stats::approxfun(dens$x, cumsum(dens$y)/sum(dens$y))
  bg_integral <- gamma * bg_cdf_fun(realiz$t)
  
  dt <- outer(realiz$t, realiz$t, "-")
  dt[upper.tri(dt, diag = TRUE)] <- NA_real_
  Gmat <- matrix(hawkes_kernel_cdf(as.vector(dt), kernel, beta, c_pl), nrow=nrow(realiz))
  Gmat[is.na(Gmat)] <- 0
  
  incremental <- bg_integral + K * rowSums(Gmat)
  incremental
}


#' @title Simulate a univariate Hawkes process (branching structure)
#'
#' @description
#' Generate event times from a Hawkes process with exponential triggering kernel
#' using a branching (cluster) representation:
#' - \eqn{\mu} controls the background (immigrant) events.
#' - \eqn{K} is the mean number of offspring per event.
#' - \eqn{\beta} is the exponential decay rate for child arrival times after the parent.
#'
#' @param mu Numeric > 0. Background rate.
#' @param K Numeric >= 0. Mean number of children per event (branching ratio). 
#'          Typically K < 1 for a subcritical process.
#' @param beta Numeric > 0. Decay rate of the exponential triggering kernel.
#' @param T Numeric > 0. Maximum time horizon. Events are restricted to [0, T].
#' @param seed Optional. Set an integer random seed for reproducibility. Default: \code{NULL}.
#'
#' @return A numeric vector of sorted event times within [0, T].
#'
#' @details
#' **Algorithm**:
#' 1. Draw background events (immigrants) from a Poisson(\eqn{\mu \times T}) process 
#'    and place them uniformly in [0, T].
#' 2. For each event at time \eqn{t_p}, draw \eqn{N_p \sim \mathrm{Poisson}(K)} children.
#'    Each child's time is \eqn{t_c = t_p + \Delta}, where \eqn{\Delta \sim \mathrm{Exp}(\beta)}.
#'    Keep only those \eqn{t_c \le T}.
#' 3. Each child then serves as a parent to further offspring, recursively, until no new events fall in [0,T].
#'
#' This method gives the same distribution as a Hawkes process with intensity
#' \eqn{\lambda(t) = \mu + \sum_{t_i < t} K \beta e^{-\beta (t - t_i)}}, but may be faster
#' if you only need the final set of event times (particularly for subcritical K).
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' times <- simulate_hawkes_branching(mu = 0.2, K = 0.5, beta = 1, T = 10)
#' print(times)
#' }
#'
#' @export
simulate_hawkes_branching <- function(mu, K, beta, T, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(mu > 0, beta > 0, K >= 0, T > 0)
  
  # 1) Generate background (immigrant) events in [0, T]
  #    N0 ~ Poisson(mu * T), times ~ Uniform(0, T)
  N0 <- rpois(1, lambda = mu * T)
  if (N0 > 0) {
    bg_times <- sort(runif(N0, min = 0, max = T))
  } else {
    bg_times <- numeric(0)
  }
  
  # We'll store all events in a growing list (start with background).
  # BFS approach: each event can spawn child events.
  events <- bg_times
  i <- 1   # index of "parent" event we're branching from
  
  while (i <= length(events)) {
    parent_time <- events[i]
    # 2) Number of children from this parent ~ Poisson(K)
    num_children <- rpois(1, K)
    if (num_children > 0) {
      # 3) Offspring arrival times are Exp(beta) after parent_time
      offsets <- rexp(num_children, rate = beta)
      child_times <- parent_time + offsets
      # Keep only those children that occur before T
      child_times <- child_times[child_times <= T]
      
      if (length(child_times) > 0) {
        # Append new children to the event list (unordered for now)
        events <- c(events, child_times)
      }
    }
    i <- i + 1
  }
  
  # Sort all event times for a canonical representation
  events <- sort(events)
  return(events)
}


if(FALSE){

  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # =========================================================================
  # 1. GENERAL SIMULATOR (EXP & POWER-LAW)
  # =========================================================================
  # This replaces the specific 'simulate_hawkes_branching' to handle both kernels
  simulate_hawkes_general <- function(kernel = c("exp", "powerlaw"), 
                                      mu, K, 
                                      beta = NULL,      # For Exp
                                      p = NULL, c = NULL, # For Power Law
                                      T_max, seed = NULL) {
    
    if (!is.null(seed)) set.seed(seed)
    kernel <- match.arg(kernel)
    
    # 1. Background Events (Poisson Process)
    # Expected N = mu * T_max
    N0 <- rpois(1, lambda = mu * T_max)
    if (N0 > 0) {
      events <- sort(runif(N0, min = 0, max = T_max))
    } else {
      events <- numeric(0)
    }
    
    # 2. Branching Loop
    i <- 1
    while (i <= length(events)) {
      parent_time <- events[i]
      # Number of children ~ Poisson(K)
      num_children <- rpois(1, K)
      
      if (num_children > 0) {
        if (kernel == "exp") {
          # Exponential waiting time: Exp(beta)
          offsets <- rexp(num_children, rate = beta)
        } else {
          # Power Law waiting time: Inverse CDF of Lomax
          # u ~ U(0,1) -> dt = c * ((1-u)^(-1/(p-1)) - 1)
          u <- runif(num_children)
          offsets <- c * ((1 - u)^(-1 / (p - 1)) - 1)
        }
        
        child_times <- parent_time + offsets
        child_times <- child_times[child_times <= T_max]
        
        if (length(child_times) > 0) {
          events <- c(events, child_times)
        }
      }
      i <- i + 1
    }
    
    return(sort(events))
  }
  
  # =========================================================================
  # 2. SIMULATION STUDY A: EXPONENTIAL KERNEL
  # =========================================================================
  
  # --- Config A ---
  N_REPS_EXP <- 100
  T_WIN      <- 1000
  TRUE_MU    <- 0.2
  TRUE_K     <- 0.5
  TRUE_BETA  <- 1.5
  TRUE_GAMMA <- TRUE_MU * T_WIN # Expected background integral
  
  results_exp <- list()
  
  cat(sprintf("\n--- STARTING STUDY A: EXPONENTIAL (%d Reps) ---\n", N_REPS_EXP))
  pb <- txtProgressBar(min = 0, max = N_REPS_EXP, style = 3)
  
  for(i in 1:N_REPS_EXP) {
    # 1. Simulate Exp Data
    times <- simulate_hawkes_general("exp", mu=TRUE_MU, K=TRUE_K, beta=TRUE_BETA, T_max=T_WIN)
    
    if(length(times) > 20) {
      # 2. Fit Exp Model
      try({
        fit <- fit_temporal_hawkes(
          params_init = list(gamma = length(times)/2, beta = 1.0, K = 0.2),
          realiz = data.frame(t = times),
          windowT = c(0, T_WIN),
          kernel = "exp",
          trace = 0
        )
        
        results_exp[[length(results_exp)+1]] <- data.frame(
          Rep = i,
          Gamma = fit$par[1],
          Beta = fit$par[2], 
          K = fit$par[3]
        )
      }, silent=TRUE)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # --- Plot A: Exponential Recovery ---
  df_exp <- do.call(rbind, results_exp) %>%
    pivot_longer(cols = c(Gamma, Beta, K), names_to = "Parameter", values_to = "Estimate")
  
  truth_exp <- data.frame(
    Parameter = c("Gamma", "Beta", "K"),
    Value     = c(TRUE_GAMMA, TRUE_BETA, TRUE_K)
  )
  
  plot_a <- ggplot(df_exp, aes(x = "Exp Fit", y = Estimate)) +
    geom_boxplot(fill = "#00BFC4", alpha = 0.6, outlier.shape = 21) +
    geom_hline(data = truth_exp, aes(yintercept = Value), 
               color = "red", linetype = "dashed", size = 1) +
    facet_wrap(~Parameter, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Study A: Parameter Recovery (Exponential Kernel)",
      subtitle = "Simulated Exp -> Fitted Exp (100 Reps)",
      y = "MLE Estimate", x = ""
    )
  
  print(plot_a)
  
  
  # =========================================================================
  # 3. SIMULATION STUDY B: POWER-LAW KERNEL
  # =========================================================================
  
  # --- Config B ---
  N_REPS_PL <- 100
  TRUE_P    <- 2.5
  TRUE_C    <- 0.5
  # (Mu and K stay the same as above)
  
  results_pl <- list()
  
  cat(sprintf("\n\n--- STARTING STUDY B: POWER-LAW (%d Reps) ---\n", N_REPS_PL))
  pb <- txtProgressBar(min = 0, max = N_REPS_PL, style = 3)
  
  for(i in 1:N_REPS_PL) {
    # 1. Simulate Power-Law Data
    times <- simulate_hawkes_general("powerlaw", mu=TRUE_MU, K=TRUE_K, p=TRUE_P, c=TRUE_C, T_max=T_WIN)
    
    if(length(times) > 20) {
      # 2. Fit Power-Law Model
      try({
        fit <- fit_temporal_hawkes(
          params_init = list(gamma = length(times)/2, p = 2.0, K = 0.2, c = 0.1),
          realiz = data.frame(t = times),
          windowT = c(0, T_WIN),
          kernel = "powerlaw",
          trace = 0
        )
        
        results_pl[[length(results_pl)+1]] <- data.frame(
          Rep = i,
          Gamma = fit$par[1],
          P = fit$par[2], 
          K = fit$par[3],
          C = fit$par[4]
        )
      }, silent=TRUE)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # --- Plot B: Power-Law Recovery ---
  df_pl <- do.call(rbind, results_pl) %>%
    pivot_longer(cols = c(Gamma, P, K, C), names_to = "Parameter", values_to = "Estimate")
  
  truth_pl <- data.frame(
    Parameter = c("Gamma", "P", "K", "C"),
    Value     = c(TRUE_GAMMA, TRUE_P, TRUE_K, TRUE_C)
  )
  
  plot_b <- ggplot(df_pl, aes(x = "PowerLaw Fit", y = Estimate)) +
    geom_boxplot(fill = "#F8766D", alpha = 0.6, outlier.shape = 21) +
    geom_hline(data = truth_pl, aes(yintercept = Value), 
               color = "red", linetype = "dashed", size = 1) +
    facet_wrap(~Parameter, scales = "free_y", nrow = 1) +
    theme_bw() +
    labs(
      title = "Study B: Parameter Recovery (Power-Law Kernel)",
      subtitle = "Simulated PowerLaw -> Fitted PowerLaw (100 Reps)",
      y = "MLE Estimate", x = ""
    )
  
  print(plot_b)
}
