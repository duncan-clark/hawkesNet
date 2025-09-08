
safe_exp <- function(x) {
  y <- exp(pmin(x, 700)) 
  y[!is.finite(y)] <- 0
  y
}

init_state <- function(c = 1, tau = 0, lambda_new = 0.5,
                       start_perps = 0, start_crimes = 0, t0 = 0,
                       perps_init = NULL) {
  stopifnot(c >= 0, c <= 1, tau >= 0, lambda_new >= 0,
            start_perps >= 0, start_crimes >= 0)

  if (!is.null(perps_init)) {
    stopifnot(all(c("id", "deg", "t_last") %in% names(perps_init)))
    perps <- perps_init[order(perps_init$id), c("id","deg","t_last")] 
  } else {
    perps <- data.frame(id = integer(0), deg = integer(0), t_last = numeric(0))
  }

  list(
    params = list(c = c, tau = tau, lambda_new = lambda_new),
    perps = perps,
    crimes = data.frame(id = integer(0), t = numeric(0)),
    next_ids = list(crime_id = start_crimes + 1L,
                    perp_id  = if (nrow(perps) > 0) max(perps$id) + 1L else start_perps + 1L),
    t0 = t0
  )
}

ba_probs <- function(state, t) {
  perps <- state$perps
  nV <- nrow(perps)
  if (nV == 0) return(data.frame(id = integer(0), p_ba = numeric(0),
                                 delta = numeric(0), deg = integer(0), t_last = numeric(0)))
  c_off <- state$params$c
  tau   <- state$params$tau
  ## delta_v = (c + d_v) * exp(-tau * (t - t_last_v))
  rec <- safe_exp(-tau * pmax(0, t - perps$t_last))
  delta <- (c_off + perps$deg) * rec
  S <- sum(delta)
  p_ba <- if (S > 0) delta / S else rep(0, nV)
  data.frame(id = perps$id, p_ba = p_ba, delta = delta,
             deg = perps$deg, t_last = perps$t_last, row.names = NULL)
}

ba_bipartite_step <- function(state, t) {
  stopifnot(is.list(state), is.finite(t))
  crime_id <- state$next_ids$crime_id
  state$crimes <- rbind(state$crimes, data.frame(id = crime_id, t = t))
  state$next_ids$crime_id <- state$next_ids$crime_id + 1L
  probs <- ba_probs(state, t)
  ##  e_v ~ Bern(p_v^{..})
  if (nrow(probs) > 0) {
    e_old <- rbinom(n = nrow(probs), size = 1L, prob = probs$p_ba)
    old_hits <- probs$id[e_old == 1L]
  } else {
    old_hits <- integer(0)
  }
  ## 4 K ~ Poisson(lambda_new) 
  lambda_new <- state$params$lambda_new
  K <- rpois(1L, lambda_new)
  new_ids <- if (K > 0) seq.int(state$next_ids$perp_id, length.out = K) else integer(0)
  if (K > 0) state$next_ids$perp_id <- state$next_ids$perp_id + K

  edges_event <- data.frame(
    crime_id = integer(0), perp_id = integer(0), is_new = logical(0)
  )
  if (length(old_hits) > 0) {
    edges_event <- rbind(edges_event,
                         data.frame(crime_id = crime_id, perp_id = old_hits, is_new = FALSE))
  }
  if (K > 0) {
    edges_event <- rbind(edges_event,
                         data.frame(crime_id = crime_id, perp_id = new_ids, is_new = TRUE))
  }
  if (length(old_hits) > 0) {
    idx <- match(old_hits, state$perps$id)
    state$perps$deg[idx] <- state$perps$deg[idx] + 1L
    state$perps$t_last[idx] <- t
  }
  ## New perps
  if (K > 0) {
    state$perps <- rbind(state$perps,
                         data.frame(id = new_ids, deg = rep(1L, K), t_last = rep(t, K)))
  }

  list(
    state = state,
    event = list(t = t, crime_id = crime_id, K = K, edges = edges_event,
                 probs = probs)
  )
}

ba_bipartite_simulate <- function(t_vec,
                                  c = 1, tau = 0, lambda_new = 0.5,
                                  start_perps = 0, start_crimes = 0, t0 = 0,
                                  perps_init = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(t_vec) == 0 || all(diff(t_vec) > 0))

  st <- init_state(c = c, tau = tau, lambda_new = lambda_new,
                   start_perps = start_perps, start_crimes = start_crimes,
                   t0 = t0, perps_init = perps_init)

  events <- vector("list", length(t_vec))
  all_edges <- list()

  for (i in seq_along(t_vec)) {
    res <- ba_bipartite_step(st, t = t_vec[i])
    st <- res$state
    events[[i]] <- list(t = res$event$t,
                        crime_id = res$event$crime_id,
                        K = res$event$K)
    if (nrow(res$event$edges) > 0) all_edges[[length(all_edges) + 1L]] <- res$event$edges
  }

  edges_df <- if (length(all_edges) > 0) do.call(rbind, all_edges) else
    data.frame(crime_id = integer(0), perp_id = integer(0), is_new = logical(0))

  list(
    state = st,
    events = do.call(rbind, lapply(events, as.data.frame)),
    edges = edges_df
  )
}


loglik <- function(t_vec, events_df, edges_df,
                             c = 0.5, tau = 0.1, lambda_new = 1,
                             start_perps = 0, start_crimes = 0) {

  stopifnot(all(c("crime_id","t","K") %in% names(events_df)))
  stopifnot(all(c("crime_id","perp_id","is_new") %in% names(edges_df)))
  events_df$crime_id <- as.integer(as.character(events_df$crime_id))
  edges_df$crime_id  <- as.integer(as.character(edges_df$crime_id))
  edges_df$perp_id   <- as.integer(as.character(edges_df$perp_id))
  perps <- data.frame(
    id = integer(0),
    deg = integer(0),
    t_last = numeric(0)
  )

  ll <- 0

  for (i in seq_along(t_vec)) {
    t <- t_vec[i]
    crime <- events_df$crime_id[i]
    
    sub_edges <- edges_df[edges_df$crime_id == crime, ]
    old_perps <- sub_edges$perp_id[sub_edges$is_new == FALSE]
    new_perps <- sub_edges$perp_id[sub_edges$is_new == TRUE]
    ##old perp edges log-likelihood
    if (nrow(perps) > 0) {
      delta <- (c + perps$deg) * exp(-tau * (t - perps$t_last))
      if (sum(delta) > 0) {
        p <- delta / sum(delta)
      } else {
        p <- rep(0, nrow(perps))
      }
      is_old <- perps$id %in% old_perps
      ll <- ll + sum(log(p[is_old] + 1e-12)) + sum(log(1 - p[!is_old] + 1e-12))
    }
    ## new perp count log-likelihood
    K <- length(new_perps)
    ll <- ll + dpois(K, lambda_new, log = TRUE)
    ## update
    if (length(old_perps) > 0 && nrow(perps) > 0) {
      perps$deg[perps$id %in% old_perps] <- perps$deg[perps$id %in% old_perps] + 1
      perps$t_last[perps$id %in% old_perps] <- t
    }
    if (K > 0) {
      max_id <- ifelse(nrow(perps) == 0, 0, max(perps$id))
      ## ensure sequential IDs
      new_ids <- (max_id + 1L):(max_id + K)
      perps <- rbind(perps, data.frame(id = new_ids, deg = rep(1,K), t_last = rep(t,K)))
    }
  }

  return(ll)
}

negloglik <- function(par, t_vec, events_df, edges_df,
                      start_perps = 0, start_crimes = 0) {
  c     <- par[1]
  tau   <- par[2]
  lam   <- par[3]
  ll <- loglik(t_vec, events_df, edges_df,
                        c = c, tau = tau, lambda_new = lam,
                        start_perps = start_perps,
                        start_crimes = start_crimes )
  return(-ll)
}


## t_vec <- cumsum(rexp(50, rate = 0.5))  
## sim <- ba_bipartite_simulate(t_vec, c = 0.2, tau = 0.1, lambda_new = 0.6,
##                              start_perps = 0, start_crimes = 0, seed = 42)
## str(sim$events)
## head(sim$edges)
## head(sim$state$perps)


## ll <- loglik(t_vec = t_vec,
##                       events_df = sim$events,
##                       edges_df  = sim$edges,
##                       c = 0.2, tau = 0.1, lambda_new = 0.6)
## ll

## fit <- optim(
##   par = c(c = 0.5, tau = 0.1, lambda_new = 0.5),  
##   fn  = negloglik,
##   method = "L-BFGS-B",
##   lower = c(0, 0, 0),  
##   upper = c(1, Inf, Inf),
##   t_vec = t_vec,
##   events_df = sim$events,
##   edges_df  = sim$edges
## )

## fit$par       
## fit$value  
