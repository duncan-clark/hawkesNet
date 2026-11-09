## ============================================================
## Power-law Hawkes (temporal) fit — KDE Background
##   - Background: mu(t) = gamma * KDE(t)
##   - Kernel: g(t) = (p-1) c^(p-1) / (t + c)^p 
##   - Optimization: L-BFGS-B with native constraints
## ============================================================

## ---- Helper: Pre-compute KDE Background ----
# Call this before fitting to get the background 'mu' vector
compute_kde_background <- function(times, windowT, bw = "nrd0") {
  .temporal_validate_data(data.frame(t = times), windowT)
  if (length(times) < 2L) stop("At least two events are needed to estimate a KDE background")
  dens <- density(times, from = windowT[1], to = windowT[2], bw = bw)
  # Normalize on the observation window: an unrestricted KDE also assigns
  # mass outside it. Integrate the same linear interpolant used for rates.
  dx <- diff(dens$x)
  mass <- sum(dx * (head(dens$y, -1L) + tail(dens$y, -1L)) / 2)
  if (!is.finite(mass) || mass <= 0) stop("KDE has no finite mass in the observation window")
  y <- dens$y / mass
  cumulative <- c(0, cumsum(dx * (head(y, -1L) + tail(y, -1L)) / 2))
  j <- findInterval(times, dens$x, all.inside = TRUE)
  offset <- times - dens$x[j]
  integral <- cumulative[j] + y[j] * offset +
    (y[j + 1L] - y[j]) * offset^2 / (2 * dx[j])
  list(mu_vec = approx(dens$x, y, xout = times)$y, total_int = 1,
       integral_at_events = integral, grid = dens$x, density = y)
}

.temporal_validate_data <- function(realiz, windowT) {
  if (!is.numeric(windowT) || length(windowT) != 2L ||
      any(!is.finite(windowT)) || windowT[2L] <= windowT[1L])
    stop("windowT must contain finite start and end with end > start")
  if (!is.data.frame(realiz) || !is.numeric(realiz$t) ||
      any(!is.finite(realiz$t)) || any(realiz$t < windowT[1L] | realiz$t > windowT[2L]))
    stop("realiz must have finite numeric event times t inside windowT")
  invisible(NULL)
}

.temporal_valid_params <- function(params, kernel) {
  is.numeric(params) && length(params) == (if (kernel == "exp") 3L else 4L) &&
    all(is.finite(params)) && params[1L] >= 0 && params[3L] >= 0 &&
    params[2L] > (if (kernel == "exp") 0 else 1) &&
    (kernel == "exp" || params[4L] > 0)
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
    out <- (p - 1) / c * exp(-p * log1p(pmax(dt, 0) / c))
  }
  out[dt <= 0] <- 0
  out
}

hawkes_kernel_cdf <- function(u, kernel = c("exp", "powerlaw"), beta = NULL, c = NULL) {
  kernel <- match.arg(kernel)
  u[u < 0] <- 0
  if (kernel == "exp") {
    if (is.null(beta) || beta <= 0) return(rep(0, length(u)))
    return(-expm1(-beta * u))
  } else {
    p <- beta
    if (is.null(p) || is.null(c) || p <= 1 || c <= 0) return(rep(0, length(u)))
    return(-expm1(-(p - 1) * log1p(u / c)))
  }
}

## ---- Updated Log-likelihood with KDE Background ----
#' @title Log-likelihood for Temporal Hawkes Process with KDE Background
#' @description Computes the log-likelihood of a Hawkes process with KDE background.
#' The built-in KDE is normalized over the observation window. Times tied on
#' the recording grid use the same pre-bin intensity and do not excite one
#' another; this convention does not fit an explicit discrete-time model.
#' @param params Numeric vector of parameters:
#'              - For "exp": c(gamma, beta, K)
#'              - For "powerlaw": c(gamma, p, K, c)
#'              gamma: scale of KDE background
#'              beta/p: decay parameter
#'              K: branching ratio
#'              c: scale parameter for power-law kernel
#' @param realiz Data frame of event times with column 't'
#' @param windowT Vector of `[start, end]` times
#' @param dists Pre-computed distance matrices (optional)
#' @param kde_bg Optional list with background values \code{mu_vec} in event-row
#'   order and its window integral \code{total_int}. If NULL, estimate the KDE.
#' @param kernel "exp" or "powerlaw"
#' @return Numeric log-likelihood value
#' @examples
#' \donttest{
#' realiz <- data.frame(t = sort(runif(50, 0, 10)))
#' # Compute KDE background
#' dens <- density(realiz$t, from = 0, to = 10)
#' mu_vec <- approx(dens$x, dens$y, xout = realiz$t)$y
#' bg <- list(mu_vec = mu_vec, total_int = 1)
#' # Log-likelihood for an exponential kernel
#' loglik_temporal_hawkes(c(1, 0.5, 0.2), realiz, c(0, 10), kde_bg = bg, kernel = "exp")
#' }
#' @export
loglik_temporal_hawkes <- function(params, realiz, windowT, dists = NULL,
                                 kde_bg = NULL, kernel = c("exp", "powerlaw")) {
  kernel <- match.arg(kernel)
  .temporal_validate_data(realiz, windowT)
  params <- unlist(params, use.names = FALSE)
  if (!.temporal_valid_params(params, kernel)) return(-Inf)
  if (is.null(kde_bg)) kde_bg <- compute_kde_background(realiz$t, windowT)
  if (!is.numeric(kde_bg$mu_vec) || length(kde_bg$mu_vec) != nrow(realiz) ||
      any(!is.finite(kde_bg$mu_vec)) || any(kde_bg$mu_vec < 0) ||
      length(kde_bg$total_int) != 1L || !is.finite(kde_bg$total_int) || kde_bg$total_int < 0)
    stop("kde_bg must contain nonnegative mu_vec matching event rows and a finite total_int")
  
  # Params: gamma (KDE scale), beta/p (decay), K (branching), [c (powerlaw)]
  gamma <- params[1]
  beta  <- params[2]
  K     <- params[3]
  c_pl  <- if (kernel == "powerlaw") params[4] else NULL
  
  # Intensity at event times: lambda(tj) = gamma*KDE(tj) + K * sum g(tj-ti)
  # kde_bg$mu_vec must be pre-calculated for these exact event times
  mu_term <- gamma * kde_bg$mu_vec
  
  if (is.null(dists)) {
    time_dist <- outer(realiz$t, realiz$t, "-")
  } else {
    time_dist <- dists$time_dist
    if (!is.matrix(time_dist) || !identical(dim(time_dist), c(nrow(realiz), nrow(realiz))) ||
        any(!is.finite(time_dist)))
      stop("dists$time_dist must be a finite square matrix matching the number of events")
  }
  
  kern_mat <- hawkes_kernel_density(dt = time_dist, kernel = kernel, beta = beta, c = c_pl)
  kern_mat[time_dist <= 0] <- 0
  
  g_sum <- rowSums(kern_mat)
  lambdas <- mu_term + K * g_sum
  
  if (any(lambdas <= 0) || any(!is.finite(lambdas))) return(-Inf)
  
  # Integral term: gamma * ∫KDE + K * sum G(T-ti)
  # The supplied integral matches the rate values; the built-in KDE is normalized.
  u <- windowT[2] - realiz$t
  G_vals <- hawkes_kernel_cdf(u, kernel = kernel, beta = beta, c = c_pl)
  
  int_lam <- gamma * kde_bg$total_int + K * sum(G_vals)
  
  loglik <- sum(log(lambdas)) - int_lam
  return(loglik)
}

#' @title Fit Temporal Hawkes Process
#' @description Fits a Hawkes process with KDE background and specified kernel.
#' Default bounds restrict the branching ratio to less than one. Custom bounds
#' can relax this restriction for finite observation windows.
#' @param params_init List of initial parameters
#' @param realiz Data frame of event times
#' @param windowT Vector of `[start, end]` times
#' @param method Optimization method (default "L-BFGS-B")
#' @param kernel "exp" or "powerlaw"
#' @param maxit Maximum number of iterations for the optimizer (default 200).
#' @param trace Non-negative integer controlling optimizer output (default 0, silent).
#' @param low Optional numeric vector of lower bounds for L-BFGS-B.
#' @param upp Optional numeric vector of upper bounds for L-BFGS-B.
#' @return An \code{optim} result with \code{par} and \code{value}, plus
#'   \code{hessian}, \code{standard_errors}, \code{boundary} and
#'   \code{hessian_status}. Boundary or unconverged fits have an NA Hessian
#'   and no Wald standard errors. Interior fits use central differences that
#'   remain inside the parameter bounds; standard errors are withheld if the
#'   observed information is not positive definite and well conditioned.
#' @examples
#' \donttest{
#' realiz <- data.frame(t = sort(runif(50, 0, 10)))
#' # Fit a temporal Hawkes process with exponential kernel
#' fit <- fit_temporal_hawkes(c(1, 0.5, 0.2), realiz, c(0, 10), kernel = "exp")
#' fit$par
#' }
#' @export
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
    # Invalid model parameters have log likelihood -Inf.
    optim_lower <- -Inf
    optim_upper <- Inf
  } else {
    # L-BFGS-B uses these bounds strictly
    optim_lower <- low
    optim_upper <- upp
  }
  
  fit <- optim(
    par = unlist(params_init),
    fn = loglik_temporal_hawkes,
    method = method,
    lower = optim_lower, # Pass NULL/-Inf if Nelder-Mead
    upper = optim_upper, # Pass NULL/Inf if Nelder-Mead
    control = list(fnscale = -1, trace = trace, maxit = maxit),
    realiz = realiz,
    windowT = windowT,
    dists = dists,
    kde_bg = bg,
    kernel = kernel,
    hessian = FALSE
  )
  # optim's unconstrained Hessian can probe a negative branching ratio even
  # when L-BFGS-B has correctly converged to its lower bound. Assess interior
  # information separately and do not turn invalid likelihoods into penalties.
  par <- fit$par
  labels <- names(par)
  if (is.null(labels)) labels <- if (kernel == "exp") c("gamma", "beta", "K") else c("gamma", "p", "K", "c")
  physical_lower <- if (kernel == "exp") c(0, 0, 0) else c(0, 1, 0, 0)
  info_lower <- pmax(physical_lower, optim_lower)
  info_upper <- rep_len(optim_upper, length(par))
  tolerance <- 1e-6 * pmax(1, abs(par))
  at_boundary <- par - info_lower <= tolerance | info_upper - par <= tolerance
  fit$boundary <- labels[at_boundary]
  fit$hessian <- matrix(NA_real_, length(par), length(par), dimnames = list(labels, labels))
  fit$standard_errors <- stats::setNames(rep(NA_real_, length(par)), labels)
  if (any(at_boundary)) {
    fit$hessian_status <- "boundary estimate; Wald information not reported"
  } else if (fit$convergence != 0L) {
    fit$hessian_status <- "optimizer did not converge; Wald information not reported"
  } else {
    step <- pmin(1e-4 * pmax(1, abs(par)), (par - info_lower) / 4, (info_upper - par) / 4)
    hessian <- tryCatch(stats::optimHess(par, loglik_temporal_hawkes,
      realiz = realiz, windowT = windowT, dists = dists, kde_bg = bg, kernel = kernel,
      control = list(ndeps = step)), error = function(e) NULL)
    if (is.null(hessian) || any(!is.finite(hessian))) {
      fit$hessian_status <- "finite-difference information unavailable"
    } else {
      fit$hessian <- hessian
      information <- -hessian
      eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
      if (min(eigenvalues) <= 0 || min(eigenvalues) / max(eigenvalues) < 1e-10) {
        fit$hessian_status <- "observed information is not positive definite or is ill conditioned"
      } else {
        fit$standard_errors <- stats::setNames(sqrt(diag(solve(information))), labels)
        fit$hessian_status <- "interior observed information"
      }
    }
  }
  fit
}

#' @title Compensator for Temporal Hawkes Process
#' @description Computes the compensator values at each event time.
#'              Optionally uses a KDE for the background rate, otherwise assumes Uniform.
#' @param params Numeric vector of parameters:
#'              - For "exp": c(gamma, beta, K)
#'              - For "powerlaw": c(gamma, p, K, c)
#'              gamma: Total background events (scale of KDE or Uniform mass)
#'              beta/p: decay parameter
#'              K: branching ratio
#'              c: scale parameter for power-law kernel
#' @param realiz Data frame of event times with column 't'
#' @param windowT Vector of `[start, end]` times
#' @param kernel "exp" or "powerlaw"
#' @param use_kde Logical or 0/1. If TRUE, uses KDE background. Default is 0 (FALSE).
#' @return Numeric vector of compensator values at each event time
#' @examples
#' \donttest{
#' realiz <- data.frame(t = sort(runif(50, 0, 10)))
#' # Compensators for a fitted model
#' comp <- compensator_temporal_hawkes(c(1, 0.5, 0.2), realiz, c(0, 10), kernel = "exp")
#' plot(comp, type = "s")
#' }
#' @export
compensator_temporal_hawkes <- function(params,
                                        realiz,
                                        windowT,
                                        kernel = c("exp", "powerlaw"),
                                        use_kde = FALSE
) {
  kernel <- match.arg(kernel)
  .temporal_validate_data(realiz, windowT)
  params <- unlist(params, use.names = FALSE)
  if (!.temporal_valid_params(params, kernel)) stop("Invalid temporal Hawkes parameters for ", kernel, " kernel")
  realiz <- realiz[order(realiz$t),,drop = FALSE]
  
  # Parse Parameters
  gamma <- params[1]
  beta  <- params[2]
  K     <- params[3]
  c_pl  <- if (kernel == "powerlaw") params[4] else NULL
  
  # --- 1. Calculate Background Integral ---
  # If use_kde is TRUE, we integrate the density.
  # If use_kde is FALSE, we integrate a Uniform(windowT) distribution.
  if (use_kde) {
    # KDE integral up to each ti (Cumulative Density Function of the KDE)
    bg_integral <- gamma * compute_kde_background(realiz$t, windowT)$integral_at_events
  } else {
    # Uniform Background: Rate = gamma / (T_end - T_start)
    # Integral(t) = gamma * (t - T_start) / (T_end - T_start)
    T_start <- windowT[1]
    T_end   <- windowT[2]
    bg_integral <- gamma * (realiz$t - T_start) / (T_end - T_start)
  }
  
  # --- 2. Calculate Triggering Kernel Integral ---
  dt <- outer(realiz$t, realiz$t, "-")
  dt[upper.tri(dt, diag = TRUE)] <- NA_real_
  
  # Calculate CDF of the kernel for every pair
  Gmat <- matrix(hawkes_kernel_cdf(as.vector(dt), kernel, beta, c_pl), nrow=nrow(realiz))
  Gmat[is.na(Gmat)] <- 0
  
  # Total compensator = Background Integral + Branching Integral
  incremental <- bg_integral + K * rowSums(Gmat)
  incremental
}

#' @title KS Test for Temporal Hawkes Process
#' @description Does KS test on pure temporal hawkes
#' @param realiz list with elements, n,lon,lat,t
#' @param windowT vector with elements, start,end
#' @param hawkes_par list with elements,   mu,alpha,beta,K
#' @param kernel "exp" or "powerlaw"
#' @param use_kde Logical/Numeric. Default 0 (FALSE).
#' @return p-value of KS test
#' @examples
#' \donttest{
#' realiz <- data.frame(t = sort(runif(50, 0, 10)))
#' # KS test for a fitted model
#' ks_test_pval_temporal(realiz, c(0, 10), c(1, 0.5, 0.2), kernel = "exp")
#' }
#' @export
ks_test_pval_temporal <- function(realiz,
                                  windowT,
                                  hawkes_par,
                                  kernel = c("exp", "powerlaw"),
                                  use_kde = FALSE
){
  compensators <- compensator_temporal_hawkes(
    params = unlist(hawkes_par),
    realiz = realiz,
    windowT = windowT,
    kernel = kernel,
    use_kde = use_kde
  )
  compensator_incs <- diff(c(0, compensators))
  
  # Transform to uniform using the time-rescaling theorem
  test_dist <- 1 - exp(-compensator_incs)
  
  test <- ks.test(test_dist, "punif")
  return(test$p.value)
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
#' @param T Numeric > 0. Maximum time horizon. Events are restricted to
#' @param seed Optional. Set an integer random seed for reproducibility. Default: \code{NULL}.
#'
#' @return A numeric vector of sorted event times within `[0, T]`.
#' @examples
#' # Simulate a Hawkes process with branching ratio 0.5
#' t <- simulate_hawkes_branching(mu = 1, K = 0.5, beta = 2, T = 10)
#' length(t)
#' @details
#' **Algorithm**:
#' 1. Draw background events (immigrants) from a Poisson(\eqn{\mu \times T}) process 
#'    and place them uniformly in `[0, T]`.
#' 2. For each event at time \eqn{t_p}, draw \eqn{N_p \sim \mathrm{Poisson}(K)} children.
#'    Each child's time is \eqn{t_c = t_p + \Delta}, where \eqn{\Delta \sim \mathrm{Exp}(\beta)}.
#'    Keep only those \eqn{t_c \le T}.
#' 3. Each child then serves as a parent to further offspring, recursively, until no new events fall in `[0,T]`.
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
  
  # 1) Generate background (immigrant) events in `[0, T]`
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
