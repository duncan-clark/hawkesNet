#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param params PARAM_DESCRIPTION
#' @param realiz PARAM_DESCRIPTION
#' @param windowT PARAM_DESCRIPTION
#' @param dists PARAM_DESCRIPTION, Default: NULL
#' @param density_approx PARAM_DESCRIPTION, Default: TRUE
#' @param numeric_integral PARAM_DESCRIPTION, Default: FALSE
#' @param cl PARAM_DESCRIPTION, Default: NULL
#' @param optimized PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[parallel]{clusterApply}}
#' @rdname loglik_temporal_hawk
#' @export
loglik_temporal_hawk = function(params,
                                realiz,
                                windowT,
                                dists = NULL,
                                density_approx = TRUE,
                                numeric_integral= FALSE,
                                cl = NULL,
                                optimized = T
){
  tval <- windowT[2]-windowT[1]
  max_t <- max(realiz$t)
  if(tval < max(realiz$t) - min(realiz$t)){
    stop("realization has points outside time window")
  }
  mu<-params[1]
  beta<-params[2]
  K<-params[3]

  # don't allow negative parameters
  if(min(mu,K,beta)<0) return(-999999)
  # don't allow explosive growth
  if(K>.99999) return(-999999)

  # Precompute distances if not supplied:
  if(is.null(dists)){
    realiz <- realiz[order(realiz$t),]
    time_dist <- outer(realiz$t, realiz$t, "-")
  }else{
    time_dist <- dists$time_dist
  }
  exp_decay_mat <- exp(-beta * time_dist)
  # for stability
  exp_decay_mat[upper.tri(exp_decay_mat)] <- 0

  # Adjust for areas with no background intensity
  adjust_factor <- 1
  if(density_approx){
    # For intlam to work according to
    # Facilitated estimation of ETAS. Bulletin of the Seismological Society of America, 103(1), 601-605
    # We need \int_{xy} \mu(x,y,t) dx,dy = \mu T
    # This means we need some rescaling of mu in later in the code
    # that is if the window is not of area 1 the "mu" in the likelihood is rescaled below
    # K is E # points, alpha and beta from exponential decay
    # note that we need alpha/pi to make sure the triggering function is a density
    mu_star <- mu
    intlam <- adjust_factor*mu*tval + K*dim(realiz)[1]
    const  <- K*beta
  }else{
    mu_star <- mu
    intlam <- adjust_factor*mu*tval
    const  <- K*beta
    if(numeric_integral){
      const  <- K*beta
      int_func <- function(x,realiz){
        # must take arg x, x is vector c(t,x,y)
        t<- x[1]
        t_diff <- t - realiz$t
        time_dist <- t_diff
        # cut down the dist matrices based on integral
        result <- sum(exp(-beta * t_diff))
        if((result==0)){
          return(1e-10)
        }
        return(K*result)
      }
      func <- function(i){
        max_t <- realiz$t[i]
        min_t <- realiz$t[i-1]
        hcubature(int_func,
                  realiz = realiz[realiz$t==min_t,],
                  lowerLimit = c(min_t, windowS$xrange[1], windowS$yrange[1]),
                  upperLimit = c(max_t, windowS$xrange[2], windowS$yrange[2])
        )$integral
      }
      if(!is.null(cl)){
        # pass the func to the cluster:
        parallel::clusterExport(cl, c("func",
                                      "int_func",
                                      "hcubature",
                                      "realiz",
                                      "x_diff",
                                      "y_diff",
                                      "t_diff",
                                      "space_dist",
                                      "time_dist"),
                                envir = environment())
        pieces <- parSapply(cl = cl,X = 2:dim(realiz)[1],FUN = func)
      }else{
        pieces <- sapply(2:dim(realiz)[1],func)
      }
    }else{
      func_3 <- function(x){
        t_comp <- (1-exp(-beta*(max_t-t)))
        return(t_comp)
      }
      pieces <- func_3(realiz$t)
    }
    if(K!=0 & (dim(realiz)[1] !=1)){
      intlam <- intlam + K*sum(pieces)
    }
  }
  # initialize log sum:
  sum_log <- log(mu_star)
  lamjs <- c()

  if(optimized){
    # Then
    gij_vec <- rowSums( # sum across columns up to j - 1
      lower.tri(exp_decay_mat, diag = FALSE) * exp_decay_mat
    )
    lamjs <- mu_star + const * gij_vec
    lamjs <- lamjs[2:length(lamjs)]
  }else{
    if(nrow(realiz) >= 2){
      for (j in 2:nrow(realiz)){
        # only points in the past can trigger the current point
        gij <- sum(exp_decay_mat[j, 1:(j-1)])
        lamjs[j-1] <- mu_star + const * gij
      }
    }else{
      sum_log <- 0
    }
  }
  # add the log lamjs in:
  if(any(is.na(lamjs)) || any(lamjs < 0)){
    return(-999999)
  }else{
    sum_log <- sum_log + sum(log(lamjs))
  }
  loglik <- sum_log - intlam
  if(loglik == -Inf){
    return(-999999)
  }
  return(loglik)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param params_init PARAM_DESCRIPTION
#' @param realiz PARAM_DESCRIPTION
#' @param windowT PARAM_DESCRIPTION
#' @param trace PARAM_DESCRIPTION, Default: 0
#' @param maxit PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname fit_temporal_hawkes
#' @export
fit_temporal_hawkes <- function(params_init,
                                realiz,
                                windowT,
                                trace = 0,
                                maxit,
                                ...){
  if(inherits(params_init, "list")){params_init <- unlist(params_init)}
  realiz <- realiz[order(realiz$t),]
  time_dist <- outer(realiz$t, realiz$t, "-")
  dists <- list(time_dist = time_dist)

  fit <- optim(par = params_init,
               fn = loglik_temporal_hawk,
               method = "Nelder-Mead", # since no hessian is available
               control = list(fnscale = -1,
                              trace=trace,
                              maxit=maxit),
               realiz = realiz,
               windowT = windowT,
               dists = dists,
               hessian = T,
               ...)
  return(fit)
}

# Compensator function of spatial temporal hawkes - for RCT theorem
# roxgen documentation
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param params PARAM_DESCRIPTION
#' @param realiz PARAM_DESCRIPTION
#' @param windowT PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname compensator_temporal_hawkes
#' @export
compensator_temporal_hawkes <- function(params,
                                        realiz,
                                        windowT){
  # make realiz
  realiz <- realiz[order(realiz$t),]
  # make realiz start at time 0
  realiz$t <- realiz$t - windowT[1]
  realiz <- realiz[realiz$t >= 0,]

  # get params
  mu<-params[1]
  beta<-params[2]
  K<-params[3]
  intlam <- mu

  # do analytic method:
  t_comps <- outer(
    X = realiz$t,
    Y = realiz$t,
    FUN = function(x,y){
      (x >= y)*(1 - exp(-beta * (x - y)))
    }
  )
  t_comps[is.na(t_comps)] <- 0
  pieces <- rowSums(t_comps)
  incremental <- realiz$t * intlam + K * pieces
  return(incremental)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param realiz PARAM_DESCRIPTION
#' @param windowT PARAM_DESCRIPTION
#' @param hawkes_par PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ks_test_pval_temporal
#' @export
ks_test_pval_temporal <- function(realiz,
                                  windowT,
                                  hawkes_par
){
  compensators <- compensator_temporal_hawkes(params = unlist(hawkes_par),
                                              realiz = realiz,
                                              windowT = windowT)
  compensator_incs <- diff(compensators)
  test_dist <- 1 - exp(-compensator_incs)
  # hist(test_dist)
  # print(test$p.value)
  test <- ks.test(test_dist,"punif")
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
#' @param T Numeric > 0. Maximum time horizon. Events are restricted to  \code{[0, T]}.
#' @param seed Optional. Set an integer random seed for reproducibility. Default: \code{NULL}.
#'
#' @return A numeric vector of sorted event times within \code{[0, T]}.
#'
#' @details
#' **Algorithm**:
#' 1. Draw background events (immigrants) from a Poisson(\eqn{\mu \times T}) process 
#'    and place them uniformly in \code{[0, T]}.
#' 2. For each event at time \eqn{t_p}, draw \eqn{N_p \sim \mathrm{Poisson}(K)} children.
#'    Each child's time is \eqn{t_c = t_p + \Delta}, where \eqn{\Delta \sim \mathrm{Exp}(\beta)}.
#'    Keep only those \eqn{t_c \le T}.
#' 3. Each child then serves as a parent to further offspring, recursively, until no new events fall in \code{[0, T]}.
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

