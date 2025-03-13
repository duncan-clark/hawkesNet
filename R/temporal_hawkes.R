

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
#' @importFrom parallel clusterExport
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
  print(loglik)
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
  if(class(params_init) == "list"){params_init <- unlist(params_init)}
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
  hist(test_dist)
  test <- ks.test(test_dist,"punif")
  print(test$p.value)
  return(test$p.value)
}
