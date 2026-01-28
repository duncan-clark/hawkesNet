

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net_list PARAM_DESCRIPTION
#' @param times PARAM_DESCRIPTION
#' @param adjust PARAM_DESCRIPTION, Default: 10
#' @param file PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname make_network_growth_animation
#' @export
make_network_growth_animation <- function(net_list,
                                          times,
                                          adjust = 10,
                                          file = NULL){

  for(i in 1:length(net_list)){
    set.vertex.attribute(net_list[[i]],'vertex.names',1:(net_list[[i]] %n% 'n'))
  }

  animate <- networkDynamic(network.list=net_list,
                            onsets = adjust*times,
                            termini = adjust*c(times[2:length(times)],max(times)),
                            vertex.pid = "vertex.names"
  )
  render.animation(animate,
                   render.par = list(tween.frames = 1,
                                     show.time = TRUE,
                                     show.stats = NULL,
                                     extraPlotCmds=NULL,
                                     initial.coords=0),
                   displaylabels = FALSE,
                   displayisolates = TRUE
  )
  if(!is.null(file)){
    saveVideo(ani.replay(),
              video.name=file,
              other.opts="-b 5000k",
              clean=TRUE)

  }
  return(animate)
}

# Helper function for edge hash maps
# Roxygen docs:
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param i PARAM_DESCRIPTION
#' @param j PARAM_DESCRIPTION
#' @param edge_hash PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname has_edge
#' @export
has_edge <- function(i, j,edge_hash) {
  existing_keys <- keys(edge_hash)
  keys_to_check <- paste(i, j, sep = "-")
  keys_to_check %in% existing_keys
}


# roxgen documentation
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param events_list PARAM_DESCRIPTION
#' @param net PARAM_DESCRIPTION, Default: NULL
#' @param directed PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[network]{network}}, \code{\link[network]{attribute.methods}}, \code{\link[network]{add.vertices}}
#' @rdname events_to_net
#' @export
#' @importFrom network network set.vertex.attribute add.vertices set.edge.attribute
events_to_net <- function(events_list,
                          net = NULL,
                          directed = F){

  for(k in 1:length(events_list$i)){
    if(k==1 & is.null(net)){
      net <- network::network(matrix(c(events_list$i[k],events_list$j[k]),nrow = 1),directed = directed)
      network::set.vertex.attribute(net,"time",events_list$t[k],v=events_list$i[k])
      network::set.vertex.attribute(net,"time",events_list$t[k],v=events_list$j[k])
    }else{
      N <- net %n% 'n'
      over_i <- events_list$i[k] - N
      over_j <- events_list$j[k] - N
      if(over_i > 0){
        net <- network::add.vertices(net,over_i)
        network::set.vertex.attribute(net,"time",events_list$t[k],v=events_list$i[k])
        # since network is now bigger amend the over j
        over_j <- over_j - 1
      }
      if(over_j > 0){
        net <- network::add.vertices(net,over_j)
        network::set.vertex.attribute(net,"time",events_list$t[k],v=events_list$j[k])
      }
      add.edge(net, events_list$i[k], events_list$j[k])
    }
    e <- get.dyads.eids(net,
                        events_list$i[k],
                        events_list$j[k])
    network::set.edge.attribute(net,
                                "time",
                                events_list$t[k],
                                e = e[[1]]
    )
  }
  set.network.attribute(net,'n',max(c(events_list$i,events_list$j)))
  return(net)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net PARAM_DESCRIPTION
#' @param t PARAM_DESCRIPTION
#' @param equals PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname filtration_to_net
#' @export
filtration_to_net <- function(net,
                              t,
                              equals = FALSE){
  # make sure to leave one less edge or vertex that if equals
  delete.edges(net, which(get.edge.attribute(net,"time")>t))
  delete.vertices(net,which(get.vertex.attribute(net,"time")>t))
  if(!equals){
    e_times <- get.edge.attribute(net,"time")
    n_times <- get.vertex.attribute(net,"time")
    t_to_delete <- max(c(e_times, n_times))
    delete.edges(net,which(e_times == t_to_delete))
    delete.vertices(net,which(n_times == t_to_delete))
  }
  # no need for vertex names
  delete.vertex.attribute(net,'vertex.names')
  return(net)
}

filtration_to_net_both <- function(net,t){
  e_times <- get.edge.attribute(net,"time")
  n_times <- get.vertex.attribute(net,"time")
  
  times <- c(e_times,n_times)
  t_to_delete <- max(times[times <= t])
  net_less <- net
  
  del_e_less <- which(e_times>t)
  del_n_less <- which(n_times>t)
  
  del_e_eq <- setdiff(which(e_times>=t_to_delete),del_e_less) 
  del_n_eq <- setdiff(which(n_times>=t_to_delete),del_n_less)
  
  
  # make sure to leave one less edge or vertex that if equals
  delete.edges(net_less, which(e_times>t))
  delete.vertices(net_less,which(n_times>t))
  
  net_eq <- net_less
  delete.edges(net_eq, del_e_eq)
  delete.vertices(net_eq, del_n_eq)
  # no need for vertex names
  delete.vertex.attribute(net,'vertex.names')
  return(list(eq = net_eq,
              less = net_less))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname get_times
#' @export
get_times <- function(net,time_name = 'time'){
  node_times <- network::get.vertex.attribute(net,time_name)
  edge_times <- network::get.edge.attribute(net,time_name)
  return(list(node_times = node_times,
              edge_times = edge_times,
              times = sort(unique(c(node_times,edge_times)))
  ))
}

# function to plot pp on line:
pp_line_plot <- function(t,title=NULL){
  plot(c(min(t),max(t)), c(-1, 1), type = "n", yaxt = "n",
     xlab = "Value", ylab = "", main = paste0("Vector on a Number Line: ",title))
abline(h = 0, col = "gray", lwd = 2)
points(t, rep(0, length(t)), pch = 19, col = "blue", cex = 1.5)
}

plot_kde_intensity <- function(event_times,
                               bw_adjust = 0.5) {
  # event_times: numeric vector (0 to 1)
  df <- data.frame(time = event_times)
  
  ggplot(df, aes(x = time)) +
    # 1. The KDE Line (Raw Intensity)
    stat_density(
      aes(y = after_stat(density)), 
      geom = "line", 
      color = "#2c3e50", 
      size = 1.2, 
      adjust = bw_adjust
    ) +
    # 2. Add rug marks to see the actual event locations
    geom_rug(alpha = 0.4, color = "firebrick") +
    # 3. Formatting
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Estimated Event Intensity over Time",
      subtitle = paste0("KDE Line Plot (Bandwidth Adjust = ", bw_adjust, ")"),
      x = "Normalized Time (0 = Oldest, 1 = Newest)",
      y = "Intensity (Event Density)"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold")
    )
}

# --- Execution ---
# times <- network_data %v% "time_scaled"
# plot_kde_intensity(times, bw_adjust = 0.3)

get_latest_times <- function(nw) {
  el <- as.matrix.network.edgelist(nw, names = FALSE)
  times <- get.edge.attribute(nw, "time")
  
  n <- network.size(nw)
  best <- rep(-Inf, n)
  
  # Treat NA times as -Inf so they don't win
  t2 <- times
  t2[is.na(t2)] <- -Inf
  
  # One vectorized pass
  best[el[, 1]] <- pmax(best[el[, 1]], t2)
  best[el[, 2]] <- pmax(best[el[, 2]], t2)
  
  # Convert untouched nodes to NA
  best[is.infinite(best)] <- NA_real_
  
  # Fallback to node times
  node_times <- get_times(nw)$node_times
  idx <- is.na(best)
  best[idx] <- node_times[idx]
  
  best
}

#' Normalize node + edge "time" attributes to [0, 1]
#'
#' Uses get_times() to compute the global min/max across BOTH node and edge times,
#' then rescales vertex and edge "time" attributes into [0, 1].
#'
#' @param net A `network` object.
#' @param attr Name of the time attribute (default "time").
#' @param keep_na If TRUE, leave NA times as NA (default TRUE). If FALSE, error on NA.
#' @param constant_value Value to assign when all non-NA times are identical (default 0).
#' @return The network with normalized times.
#' @export
normalize_times_01 <- function(net, attr = "time", keep_na = TRUE, constant_value = 0) {
  # Pull times using your existing helper
  times_obj <- get_times(net)
  
  # Global range across node+edge times (unique & sorted already)
  all_times <- times_obj$times
  
  if (length(all_times) == 0L) return(net)
  
  if (!keep_na && anyNA(all_times)) {
    stop("NA times found; set keep_na = TRUE to keep them as NA.")
  }
  
  rng <- range(all_times, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(net)
  
  # Helper
  norm01 <- function(x) (x - rng[1]) / (rng[2] - rng[1])
  
  # Vertex times
  v_times <- network::get.vertex.attribute(net, attr)
  if (!is.null(v_times)) {
    v_times <- as.numeric(v_times)
    if (rng[2] == rng[1]) {
      v_times[!is.na(v_times)] <- constant_value
    } else {
      idx <- !is.na(v_times)
      v_times[idx] <- norm01(v_times[idx])
    }
    network::set.vertex.attribute(net, attr, v_times)
  }
  
  # Edge times
  e_times <- network::get.edge.attribute(net, attr)
  if (!is.null(e_times)) {
    e_times <- as.numeric(e_times)
    if (rng[2] == rng[1]) {
      e_times[!is.na(e_times)] <- constant_value
    } else {
      idx <- !is.na(e_times)
      e_times[idx] <- norm01(e_times[idx])
    }
    network::set.edge.attribute(net, attr, e_times)
  }
  
  net
}


# intensity_ggplot:
pp_intensity_ggplot <- function(times,
                                line_multiplier = 1,
                                base_multiplier = 1,
                                tlim = NULL,
                                dt = NULL,
                                mu = 0,
                                K = 1,
                                beta = 1,
                                spikes = TRUE,
                                smooth = TRUE,
                                title = "Point process intensity-style plot") {
  
  times <- sort(as.numeric(times))
  times <- times[is.finite(times)]
  if (!length(times)) stop("times is empty after removing non-finite values.")
  
  if (is.null(tlim)) tlim <- range(times)
  t0 <- tlim[1]; t1 <- tlim[2]
  
  if (is.null(dt)) dt <- (t1 - t0) / 2000
  dt <- max(dt, .Machine$double.eps)
  
  df_events <- data.frame(time = times)
  
  plot_layers <- list()
  
  # ---- Smooth Hawkes-style intensity ----
  if (smooth) {
    grid <- seq(t0, t1, by = dt)
    
    # Efficient recursion
    decay <- exp(-beta * dt)
    s <- numeric(length(grid))
    
    idx <- findInterval(times, grid, rightmost.closed = TRUE)
    counts <- tabulate(pmax(1L, idx), nbins = length(grid))
    
    for (k in 2:length(grid)) {
      s[k] <- decay * s[k - 1] + counts[k - 1]
    }
    
    lambda <- mu + K * s
    df_lambda <- data.frame(time = grid, intensity = lambda)
    
    plot_layers <- c(
      plot_layers,
      list(
        geom_line(
          data = df_lambda,
          aes(time, intensity),
          linewidth = 1*line_multiplier
        )
      )
    )
    
    ymax <- max(lambda)
  } else {
    ymax <- 1
  }
  
  # ---- Event spikes ----
  if (spikes) {
    plot_layers <- c(
      plot_layers,
      list(
        geom_linerange(
          data = df_events,
          aes(x = time, ymin = 0, ymax = ymax),
          alpha = 0.6
        )
      )
    )
  }
  
  ggplot() +
    plot_layers +
    coord_cartesian(xlim = c(t0, t1), ylim = c(0, ymax)) +
    labs(
      x = "time",
      y = "intensity",
      title = title
    ) +
    theme_minimal(base_size = 1*base_multiplier)
}

# Compensator fuunction of spatial temporal hawkes - for RCT theorem
# roxgen documentation
#' @description
#' @param params list with elements,   mu,alpha,beta,K
#' @param realiz list with elements, n,lon,lat,t
#' @param windowT vector with elements, start,end
#' @param windowS os.win
#' @output SOME OUTPUT
compensator_hawkes <- function(params,
                               realiz,
                               windowT,
                               windowS = NULL,
                               zero_background_region = NULL,
                               optimized = T){
  # make realiz
  realiz <- realiz[order(realiz$t),]
  # make realiz start at time 0
  realiz$t <- realiz$t - windowT[1]
  realiz <- realiz[realiz$t >= 0,]
  tval <- windowT[2]-windowT[1]
  max_t <- max(realiz$t)
  
  # get params
  mu<-params[which(names(params)=="mu")]
  alpha<-params[which(names(params)=="alpha")]
  beta<-params[which(names(params)=="beta")]
  K<-params[which(names(params)=="K")]
  
  
  
  if(is.null(zero_background_region) | is.null(windowS)){
    adjust_factor = 1
    in_zero_background <- rep(0,dim(realiz)[1])
  }else{
    area_zero_background <- spatstat.geom::area(owin(zero_background_region))
    adjust_factor <- (spatstat.geom::area.owin(windowS) - area_zero_background) / spatstat.geom::area.owin(windowS)
    in_zero_background <- 1*inside.owin(realiz[,c("x","y")],w=zero_background_region)
  }
  if(is.null(windowS)){
    mu_star <- mu
  }else{
    windowS <- as.owin(windowS)
    mu_star <- mu/(spatstat.geom::area(windowS))
  }
 
  intlam <- adjust_factor*mu
  
  # do analytic method:
  # each piece consider the integral only caring about the k-th time point
  func_2 <- function(i,max_t){
    t_comp <- (1-exp(-beta*(max_t-realiz$t[i])))
    if(!is.null(windowS)){
      s_1 <- (pnorm(windowS$xrange[2],mean = realiz$x[i],sd = 1/sqrt(2*alpha)) -
                pnorm(0,mean = realiz$x[i],sd = 1/sqrt(2*alpha)))
      s_2 <- (pnorm(windowS$yrange[2],mean = realiz$y[i],sd = 1/sqrt(2*alpha)) -
                pnorm(0,mean = realiz$y[i],sd = 1/sqrt(2*alpha)))
      return(t_comp*s_1*s_2)
    }else{
      return(t_comp)
    }
  }
  
  if(optimized){
    t_comps <- outer(
      X = realiz$t,
      Y = realiz$t,
      FUN = function(x,y){
        (x >= y)*(1 - exp(-beta * (x - y)))
      }
    )
    
    if(!is.null(windowS)){
      s_1 <- (pnorm(windowS$xrange[2],mean = realiz$x,sd = 1/sqrt(2*alpha)) -
                pnorm(0,mean = realiz$x,sd = 1/sqrt(2*alpha)))
      s_2 <- (pnorm(windowS$yrange[2],mean = realiz$y,sd = 1/sqrt(2*alpha)) -
                pnorm(0,mean = realiz$y,sd = 1/sqrt(2*alpha)))
      #s_1 and s_2 should act on columns:
      t_comps <- t_comps * matrix(s_1,
                                  byrow = T,
                                  nrow = dim(t_comps)[1],
                                  ncol = dim(t_comps)[2])
      t_comps <- t_comps * matrix(s_2,
                                  byrow = T,
                                  nrow = dim(t_comps)[1],
                                  ncol = dim(t_comps)[2])
    }
    pieces <- rowSums(t_comps)
    incremental <- realiz$t * intlam + K * pieces
  }else{
    # Unoptimized - for information
    pieces <- lapply(1:dim(realiz)[1],function(i){
      sapply(1:i,function(j){func_2(j,realiz$t[i])
      }
      )
    }
    )
    incremental <- sapply(1:length(pieces),function(i){
      realiz$t[i]*(intlam) + K*sum(pieces[[i]])
    })
  }
  return(incremental)
}

# Some documentaton:
#' @description
#' @param realiz list with elements, n,lon,lat,t
#' @param windowT vector with elements, start,end
#' @param windowS os.win
#' @param hawkes_par list with elements,   mu,alpha,beta,K
#' @param zero_background_region os.win
ks_test_pval <- function(realiz,
                         windowT,
                         hawkes_par,
                         zero_background_region,
                         windowS = NULL
){
  compensators <- compensator_hawkes(params = unlist(hawkes_par),
                                     realiz = realiz,
                                     windowT = windowT,
                                     windowS = windowS,
                                     zero_background_region = zero_background_region,
                                     optimized = T)
  compensator_incs <- diff(compensators)
  test_dist <- 1 - exp(-compensator_incs)
  test <- ks.test(test_dist,"punif")
  hist(compensator_incs)
  hist(test_dist)
  return(test$p.value)
}
