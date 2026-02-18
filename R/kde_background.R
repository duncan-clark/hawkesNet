# =============================================================================
# Inhomogeneous background rate estimation with KDE
# =============================================================================
# Estimates a time-varying baseline rate mu(t) via kernel density estimation,
# and provides integrated baseline Lambda(t). Used by fit_hawkesNet_inhom
# and prepare_inhomogeneous_background.
# =============================================================================

#' Ensure times are sorted numeric
#' @param t Numeric vector of times
#' @return Sorted numeric vector
#' @noRd
ensure_sorted <- function(t) sort(as.numeric(t))

#' Estimate inhomogeneous baseline rate mu(t) from event times using KDE
#'
#' @param t Numeric vector of event times
#' @param windowT Length-2 vector (t0, t1) defining the observation window. If NULL, uses range(t).
#' @param bw Bandwidth for density(). If NULL, uses default (Silverman).
#' @param grid_n Number of grid points for the density / rate estimate
#' @return List with: mu_fun (function mapping t -> rate), grid, mu_grid, bw, windowT.
#'         mu_grid is the rate at each grid point (density * n so that integral ≈ n).
#' @examples
#' t <- sort(c(rexp(50, 2), runif(20, 0, 10)))
#' mu_fit <- estimate_mu_kde(t, windowT = c(0, 10))
#' mu_fit$mu_fun(5)
#' @export
estimate_mu_kde <- function(t, windowT = NULL, bw = NULL, grid_n = 4096) {
  t <- ensure_sorted(t)
  if (is.null(windowT)) windowT <- c(min(t), max(t))
  t0 <- windowT[1]
  t1 <- windowT[2]
  t <- t[t >= t0 & t <= t1]
  n <- length(t)
  if (n < 5) stop("Need more events to estimate baseline reliably.")

  if (is.null(bw)) {
    bw <- 1.06 * sd(t) * n^(-1/5)
    if (!is.finite(bw) || bw <= 0) bw <- (t1 - t0) / 50
  }

  dens <- density(t, bw = bw, from = t0, to = t1, n = grid_n)
  # density integrates to 1; scale by n so that integral of rate ≈ n (total events)
  mu_grid <- dens$y * n
  mu_grid <- pmax(mu_grid, 1e-12)

  mu_fun <- approxfun(dens$x, mu_grid, rule = 2)

  list(
    mu_fun = mu_fun,
    grid = dens$x,
    mu_grid = mu_grid,
    bw = bw,
    windowT = windowT,
    n = n
  )
}

#' Integrated baseline Lambda(t) = integral_0^t mu(s) ds (cumulative hazard)
#'
#' @param mu_fit Output from estimate_mu_kde (list with grid, mu_grid)
#' @return List with Lambda_fun (function t -> Lambda(t)) and Lambda_grid (data.frame t, Lambda)
#' @examples
#' t <- sort(c(rexp(50, 2), runif(20, 0, 10)))
#' mu_fit <- estimate_mu_kde(t, windowT = c(0, 10))
#' ch <- make_cumhaz_fun(mu_fit)
#' ch$Lambda_fun(5)
#' @export
make_cumhaz_fun <- function(mu_fit) {
  x <- mu_fit$grid
  y <- mu_fit$mu_grid
  dx <- diff(x)
  area <- dx * (head(y, -1) + tail(y, -1)) / 2
  Lambda <- c(0, cumsum(area))
  Lambda_fun <- approxfun(x, Lambda, rule = 2)
  list(
    Lambda_fun = Lambda_fun,
    Lambda_grid = data.frame(t = x, Lambda = Lambda)
  )
}

#' Time-rescale event times so baseline becomes constant: tau_i = Lambda(t_i)
#'
#' Under the time-change, the rescaled process has approximately constant rate 1.
#'
#' @param t Numeric vector of event times
#' @param mu_fit Output from estimate_mu_kde
#' @return List with tau (rescaled times), Lambda_fun, Lambda_grid
#' @examples
#' t <- sort(c(rexp(50, 2), runif(20, 0, 10)))
#' mu_fit <- estimate_mu_kde(t, windowT = c(0, 10))
#' rescaled <- time_rescale_by_baseline(t, mu_fit)
#' head(rescaled$tau)
#' @export
time_rescale_by_baseline <- function(t, mu_fit) {
  t <- ensure_sorted(t)
  Lambda_obj <- make_cumhaz_fun(mu_fit)
  tau <- Lambda_obj$Lambda_fun(t)
  list(
    tau = tau,
    Lambda_fun = Lambda_obj$Lambda_fun,
    Lambda_grid = Lambda_obj$Lambda_grid
  )
}

#' Rescale vertex and edge times of a network using KDE baseline
#'
#' Gets all node and edge times from the network, estimates mu(t) via KDE,
#' computes Lambda(t), and replaces the time attribute with tau = Lambda(t).
#'
#' @param net A network object with vertex and edge time attributes
#' @param time_attr Name of the time attribute (default "time"). Must exist for vertices and edges.
#' @param bw Bandwidth for KDE; NULL = default
#' @param grid_n Grid size for KDE
#' @return List with: net_rescaled (network with time_attr set to rescaled times),
#'         time_window_rescaled = c(0, max_tau), mu_fit, Lambda_fun.
#' @examples
#' \donttest{
#' # Build a small network with time attributes
#' net <- network::network.initialize(5, directed = FALSE)
#' network::set.vertex.attribute(net, "time", c(0.1, 0.2, 0.5, 1.0, 2.0))
#' network::add.edge(net, 1, 2)
#' network::add.edge(net, 2, 3)
#' network::add.edge(net, 3, 4)
#' network::add.edge(net, 4, 5)
#' network::add.edge(net, 1, 5)
#' network::set.edge.attribute(net, "time", c(0.1, 0.3, 0.6, 1.2, 2.1))
#' res <- network_rescale_times_by_kde(net)
#' res$time_window_rescaled
#' }
#' @export
network_rescale_times_by_kde <- function(net, time_attr = "time", bw = NULL, grid_n = 4096) {
  node_times <- get.vertex.attribute(net, time_attr)
  edge_times <- get.edge.attribute(net, time_attr)
  t_all <- sort(unique(c(node_times, edge_times)))
  t_all <- t_all[!is.na(t_all)]
  if (length(t_all) < 5) stop("Need more event times for KDE.")

  windowT <- c(min(t_all), max(t_all))
  mu_fit <- estimate_mu_kde(t_all, windowT = windowT, bw = bw, grid_n = grid_n)
  Lambda_obj <- make_cumhaz_fun(mu_fit)
  Lambda_fun <- Lambda_obj$Lambda_fun

  tau_node <- Lambda_fun(node_times)
  tau_edge <- Lambda_fun(edge_times)

  net_rescaled <- network.copy(net)
  set.vertex.attribute(net_rescaled, time_attr, tau_node)
  set.edge.attribute(net_rescaled, time_attr, tau_edge)

  max_tau <- max(c(tau_node, tau_edge), na.rm = TRUE)

  list(
    net_rescaled = net_rescaled,
    time_window_rescaled = c(0, max_tau),
    mu_fit = mu_fit,
    Lambda_fun = Lambda_fun,
    Lambda_grid = Lambda_obj$Lambda_grid
  )
}

#' Plot estimated baseline rate and integrated baseline (for sanity checks)
#' @param mu_fit Output from estimate_mu_kde
#' @param rescale Optional list with Lambda_grid (for second panel)
#' @param main Plot title
#' @return Invisible \code{NULL}; called for its side effect (plotting).
#' @examples
#' t <- sort(c(rexp(50, 2), runif(20, 0, 10)))
#' mu_fit <- estimate_mu_kde(t, windowT = c(0, 10))
#' plot_kde_background(mu_fit)
#' @export
plot_kde_background <- function(mu_fit, rescale = NULL, main = "KDE baseline") {
  op <- par(mfrow = c(1, 2))
  on.exit(par(op))
  plot(mu_fit$grid, mu_fit$mu_grid, type = "l",
       xlab = "t", ylab = "mu_hat(t)", main = paste(main, "- rate"))
  if (!is.null(rescale)) {
    plot(rescale$Lambda_grid$t, rescale$Lambda_grid$Lambda, type = "l",
         xlab = "t", ylab = "Lambda(t)", main = paste(main, "- integrated"))
  }
  invisible(NULL)
}
