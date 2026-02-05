#' Mark PMF for Barabási–Albert-style (degree-weighted) attachment
#'
#' Probability mass function for the mark (new edge) given the filtration; uses degree-weighted attachment with exponential time decay.
#'
#' @param time Current event time.
#' @param params List with \code{beta_edges} (and optionally \code{beta_overall}, \code{K}, \code{mu}).
#' @param mark_filtration Observed network up to \code{time}.
#' @param mark Optional network state at \code{time}; if \code{NULL}, derived from \code{mark_filtration}.
#' @param generate_mark If \code{TRUE}, also sample a new edge (default \code{FALSE}).
#' @param new_edge_hash Optional hash of existing edges for fast lookup.
#' @return List with \code{log_mark_density}, \code{log_density_func}, and optionally sampled edge / probabilities.
#' @seealso \code{\link[network]{network}}, \code{\link[network]{add.vertices}}
#' @rdname PMF_mark_BA
#' @export
#' @importFrom network network add.vertices
PMF_mark_BA <- function(time,
                        params,
                        mark_filtration,
                        mark = NULL,
                        generate_mark = FALSE,
                        generate_density = TRUE,
                        grad = FALSE,
                        new_edge_hash = NULL,
                        truncation = NULL,
                        ...){
  
  if(is.null(mark)){
    mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
  }
  last_net <- filtration_to_net(mark_filtration, time, equals = FALSE)
  new_net <- last_net

  if(last_net %n% 'n' != 0){
    new_nodes <- mark %n% 'n'
    old_nodes <- last_net %n% 'n'
    network::add.vertices(new_net, nv = new_nodes)
    set.vertex.attribute(new_net, "time", c(get.vertex.attribute(last_net, "time"), rep(time, new_nodes)))
  } else {
    last_net <- NULL
    new_net <- network(matrix(0, 1, 1), directed = FALSE)
    set.vertex.attribute(new_net, "time", time)
    old_nodes <- 0
    new_nodes <- 1
  }

  # get the possible edges for the given truncation:
  # if no nodes have been added then :
  poss_tails <- ((old_nodes+1) : new_nodes)
  poss_tails <- poss_tails[poss_tails>0]
  poss_heads <- 1:old_nodes
  poss_heads <- poss_heads[poss_heads>0]
  poss_edges <- expand.grid(poss_tails,poss_heads)
  poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
  tails <- poss_edges[,1]
  heads <- poss_edges[,2]
  
  # only consider edges that were not already in the old network
  if(!is.null(last_net)){
    in_old_net <- sapply(seq_along(heads), function(i) {
      length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
    })
    tails <- tails[!in_old_net]
    heads <- heads[!in_old_net]
  }
  
  if(!is.null(last_net) && (last_net %n% 'n' > 2)){
    times <- get.vertex.attribute(last_net, "time")
    node_degrees <- degree(last_net)
    degs <- node_degrees * exp(-params$beta_edges * (time - times))
    total_deg <- sum(degs)
    
    if(total_deg == 0 || is.na(total_deg)){
      probs <- rep(1, length(heads))
    } else {
      probs <- degs[heads] / total_deg
    }
  } else {
    node_degrees <- NULL
    probs <- rep(1, length(heads))
  }
  
  if(!is.null(mark) && length(probs) !=0 && (last_net %n% 'n' > 2)){
    if(is.null(new_edge_hash)){
      in_mark <- sapply(seq_along(heads), function(i) {
        length(get.edgeIDs(mark, heads[i], tails[i])) != 0
      })
    } else {
      in_mark <- has_edge(heads, tails, new_edge_hash)
    }
    
    log_mark_density <- sum(log(probs[in_mark])) + sum(log(1 - probs[!in_mark]))
    mark_density <- exp(log_mark_density)
  } else {
    in_mark <- rep(1, length(heads))
    log_mark_density <- 0
    mark_density <- 1
  }
  
  # 1. Define the lightweight log-density function for BA
  log_density_func_light <- function(params) {
    # BA Logic:
    # If using BA, the edge probabilities are often just based on degrees 
    # and the beta_edges decay, calculated inside the main function.
    # RE-CALCULATING this efficiently inside a closure is tricky because 
    # BA depends on the full graph history (degrees), not just change stats.
    
    # IF you stored the calculated 'probs' and 'in_mark' relative to the specific 
    # edge set at 'time', you can use them here *assuming params don't change degrees*.
    # But usually, 'beta_edges' changes the degrees.
    
    # CRITICAL NOTE: Implementing a fast recalculation for BA inside a closure 
    # is harder than CS because you need the degrees of the whole network.
    # For now, if you just need it to run with *fixed* probs (approximate) 
    # or if you are only optimizing params that don't change structure:
    if(!is.null(node_degrees)){
      degs <- node_degrees * exp(-params$beta_edges * (time - times))
      total_deg <- sum(degs)
      if(is.na(total_deg) || total_deg == 0){
        return(0)
      }else {
        probs <- degs[heads] / total_deg
      }
      return(sum(log(probs[in_mark])) + sum(log(1 - probs[!in_mark])))
    }else{
      return(0)
    }
  }
  
  # 2. Capture the environment 
  # (You need 'probs' and 'in_mark' to be captured)
  environment(log_density_func_light) <- list2env(
    list(
      heads = heads,
      time = time,
      times = times,
      node_degrees = node_degrees,
      in_mark = in_mark
    ), 
    parent = baseenv()
  )
  
  # 3. Define the wrapper
  density_func_light <- function(params) {
    exp(log_density_func_light(params))
  }
  
  # 4. Link wrapper to the log function
  # CRITICAL: You must include 'log_density_func_light' in the environment!
  environment(density_func_light) <- list2env(
    list(log_density_func_light = log_density_func_light), 
    parent = baseenv()
  )

  if(generate_mark){
    # use latest mark as baseline:
    last_net <- mark
    times <- get.vertex.attribute(last_net, "time")
    if(!is.null(last_net) && (last_net %n% 'n') > 2){
      mark_sample <- last_net
      old_nodes <- last_net %n% 'n'
      new_nodes <- 1
      mark_sample <- network::add.vertices(mark_sample,new_nodes)
      set.vertex.attribute(mark_sample,"time",c((last_net %v% 'time'),rep(time,new_nodes)))
      new_nodes <- mark_sample %n% 'n'
      
      poss_tails <- (old_nodes+1) : (new_nodes)
      poss_tails <- poss_tails[poss_tails>0]
      poss_heads <- 1:old_nodes
      poss_heads <- poss_heads[poss_heads>0]
      poss_edges <- expand.grid(poss_tails,poss_heads)
      poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
      tails <- poss_edges[,1]
      heads <- poss_edges[,2]
      
      # only consider edges that are not in the old net
      if(!is.null(last_net)){
        in_old_net <- sapply(seq_along(heads), function(i) {
          length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
        })
        tails <- tails[!in_old_net]
        heads <- heads[!in_old_net]
      }
      degs <- degree(last_net) * exp(-params$beta_edges * (time - times))
      total_deg <- sum(degs)
      
      if(total_deg == 0){
        probs <- rep(1, length(heads))
      } else {
        probs <- degs[heads] / total_deg
      }
      
      if(any(is.na(probs))){
        browser()
      }
      
      add <- runif(length(probs)) < probs
      network::add.edges(mark_sample, heads[add], tails[add])
      
      log_mark_sample_density <- sum(log(probs[add])) + sum(log(1 - probs[!add]))
      mark_sample_density <- exp(log_mark_sample_density)
    }else{
      if(is.null(last_net)){
        mark_sample <- network::network(matrix(1),directed = F)
        set.vertex.attribute(mark_sample,"time",time)
      }else{
        mark_sample <- last_net
      }
      times <- mark_sample %v% 'time'
      mark_sample <- network::add.vertices(mark_sample,1)
      if(mark_sample %n% 'n' == 2){
        network::add.edges(mark_sample, 2, 1) # Force the first edge to create a seed
      }
      set.vertex.attribute(mark_sample,"time",c(times,time))
      
      mark_sample_density <- 1
      log_mark_sample_density <- 0
    }
  } else {
    mark_sample <- new_net
    log_mark_sample_density <- 0
    mark_sample_density <- 1
  }
  
  return(list(
    mark_density = mark_density,
    log_mark_density = log_mark_density,
    density_func = density_func_light,
    log_density_func = log_density_func_light,
    edge_probs = probs,
    # mark sample:
    mark_sample = mark_sample,
    mark_sample_density = exp(log_mark_sample_density),
    log_mark_sample_density = log_mark_sample_density
  ))
}

#' Mark PMF for change statistic (ERGM-style) attachment
#'
#' Probability mass function for the mark using an ERNM/ERGM-style model with change statistics and optional truncation.
#'
#' @param time Current event time.
#' @param params List including \code{node_lambda}, \code{CS_params}, and optionally \code{beta_edges}, \code{K}, etc.
#' @param mark_filtration Observed network up to \code{time}.
#' @param mark Optional network at \code{time}; if \code{NULL}, derived from \code{mark_filtration}.
#' @param generate_mark If \code{TRUE}, sample a new edge (default \code{FALSE}).
#' @param new_edge_hash Optional hash of existing edges for fast lookup.
#' @param formula_RHS Character RHS of the ERNM formula (e.g. \code{"edges + triangles() + star(c(2,3))"}).
#' @param grad If \code{TRUE}, compute gradient (default \code{FALSE}).
#' @param truncation Truncation window: 1 = only new-to-old edges; k = edges from k steps before new nodes.
#' @return List with \code{log_mark_density}, \code{log_density_func}, and optionally sampled edge / probabilities.
#' @seealso \code{\link[network]{network}}, \code{\link[network]{add.vertices}}, \code{\link[ernm]{as.BinaryNet}}
#' @rdname PMF_mark_CS
#' @export
#' @importFrom network network add.vertices
#' @importFrom ernm as.BinaryNet
PMF_mark_CS <- function(time,
                        params,
                        mark_filtration,
                        mark = NULL,
                        generate_mark = FALSE,
                        generate_density = TRUE,
                        new_edge_hash = NULL,
                        formula_RHS,
                        grad = FALSE,
                        truncation = 1,
                        mark_decay = 'node_entrance',
                        model = NULL,
                        max_node_time = NULL,
                        ...
){
  if(is.null(mark)){
    mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
  }
  last_net <- filtration_to_net(mark_filtration, time, equals = FALSE)
  new_net <- last_net
  if(is.null(max_node_time)){
    max_node_time <- Inf
  }
  
  if(last_net %n% 'n' != 0){
    new_nodes <- mark %n% 'n'
    old_nodes <- last_net %n% 'n'
    network::add.vertices(new_net,new_nodes-old_nodes)
    set.vertex.attribute(new_net,"time",c(last_net %v% 'time',rep(time,new_nodes)))
  }else{
    last_net <- NULL
    new_net <- network::network(matrix(1),directed = F)
    delete.vertex.attribute(new_net,'na')
    set.vertex.attribute(new_net,"time",time)
    old_nodes <- 0
    new_nodes <- 1
  }
  # get the possible edges for the given truncation:
  # if no nodes have been added then :
  poss_tails <- (old_nodes - truncation) : (new_nodes)
  poss_tails <- poss_tails[poss_tails>0]
  poss_heads <- (new_nodes - truncation-1):new_nodes
  poss_heads <- poss_heads[poss_heads>0]
  poss_edges <- expand.grid(poss_tails,poss_heads)
  poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
  tails <- poss_edges[,1]
  heads <- poss_edges[,2]

  # only consider edges that were not already in the old network
  if(!is.null(last_net) & length(heads) != 0){
    in_old_net <- sapply(seq_along(heads), function(i) {
      length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
    })
    tails <- tails[!in_old_net]
    heads <- heads[!in_old_net]
  }

  if(!is.null(last_net) & generate_density){
    if(last_net %n% 'n' > 0){
      # if new net has less than 4 nodes add some:
      if(new_net %n% 'n' < 4){
        old_new_net <- new_net
        new_net <- network::add.vertices(new_net,4 - (new_net %n% 'n'))
      }else{
        old_new_net <- new_net
      }
      
      if(max(tails)>new_net %n% 'n'){
        stop("accidently adding a edge into the network that doesn't have that node yet")
      }
      
      # delete NAs to prevent C++ using them
      delete.vertex.attribute(new_net,'na')
      if(is.null(model)){
        model <- createCppModel(as.formula(paste("new_net ~ ",formula_RHS)))
      }else{
        model$setNetwork(as.BinaryNet(new_net))
      }
      new_net <- old_new_net
      model$calculate()
      stat <- model$statistics()
      change_stats <- model$computeChangeStats(tails, heads)
      
      probs <- apply(change_stats, 1, function(c){
        1/(1+exp(-sum(c*params$CS_params)))
      })
      # use either node times or last node activity:
      if(mark_decay == 'activity'){
        node_times <- get_latest_times(new_net)
      }
      if(mark_decay == 'node_entrance'){
        node_times <- new_net %v% 'time'
      }
      diffs <- time - node_times[heads]
      factor <- exp(-params$beta_edges*(diffs))
      probs <- probs * factor
      
      if(length(probs)==0){
        in_mark <- logical(0)
        probs <- NULL
      }
      
    }else{
      change_stats <- matrix(0, nrow = 0, ncol = length(params$CS_params))
      in_mark <- logical(0)
      times <- last_net %v% 'time'
      probs <- c(1)
    }
  }else{
    change_stats <- matrix(0, nrow = 0, ncol = length(params$CS_params))
    in_mark <- logical(0)
    probs <- NULL
  }
  
  if(!is.null(mark) & !is.null(probs)){
    if(is.null(new_edge_hash)){
      in_mark <- sapply(seq_along(heads), function(i) {
        length(get.edgeIDs(mark, heads[i], tails[i])) != 0
      })
    }else{
      in_mark <- has_edge(heads,tails,new_edge_hash)
    }
    
    if(length(probs)==1){
      log_mark_density <- 0
      mark_density <-1
      mark_grad <- 0
      decay_grad <- 0
    }else{
      if(time >max_node_time){
        node_dens <- 0
      }else{
        node_dens <- log(stats::dpois(new_nodes-old_nodes,params$node_lambda))
      }
      log_mark_density <- sum(log(probs[in_mark])) + sum(log(1-probs[!in_mark])) + node_dens
      mark_density <- exp(log_mark_density)
      
      if(grad){
        # get the mark grad:
        probs_grads <- lapply(change_stats,function(c){
          (-c)*exp(-sum(c*params$CS_params))/((1+exp(-sum(c*params$CS_params)))^2)}
        )
        e <- in_mark*1
        derivs <- mapply(probs_grads,probs,in_mark,FUN = function(dp,p,e){
          tmp <- e*p^(e-1)*(1-p)^(1-e) - (1-e)*p^e*(1-p)^(-e)
          tmp <- dp*tmp
          tmp <- sign(tmp)*exp(log(abs(tmp)) - log(p))
          return(tmp)
        },SIMPLIFY =F)
        # get devided by the right prob:
        derivs <- do.call(rbind,derivs)
        mark_grad <- colSums(derivs)
        
        # get the mark grad:
        decay_grads <- -(time-times)*exp(-params$beta_edges*(time - times))
        e <- in_mark*1
        derivs <- mapply(decay_grads,probs,in_mark,FUN = function(dp,p,e){
          tmp <- e*p^(e-1)*(1-p)^(1-e) - (1-e)*p^e*(1-p)^(-e)
          tmp <- dp*tmp
          tmp <- sign(tmp)*exp(log(abs(tmp)) - log(p))
          return(tmp)
        },SIMPLIFY =F)
        # get devided by the right prob:
        derivs <- do.call(rbind,derivs)
        decay_grad <- colSums(derivs)
      }else{
        mark_grad <- 0
        decay_grad <- 0
      }
      
      
    }
  }else{
    mark_density <- 1
    mark_grad <- 0
    decay_grad <- 0
    mark_density_normalized <- NULL
    log_mark_density <- 0
  }
  
  # if node_dens doesn't exist set it to 0
  if(!exists("node_dens")){
    node_dens <- 0
  }
  
  log_density_func_light <- function(params) {
    eta <- change_stats %*% params$CS_params
    p   <- stats::plogis(eta)
    # p <- 1/(1+exp(-eta))
    log_edge_part <- sum(log(p[in_mark])) + sum(log1p(-p[!in_mark]))
    log_edge_part + node_dens
  }

  
  environment(log_density_func_light) <- list2env(
    list(
      change_stats = change_stats,
      in_mark      = in_mark,
      node_dens    = node_dens,
      plogis = stats::plogis
    ),
    parent = baseenv()
  )
  
  # define the function (will be rebound to a minimal env right after)
  log_density_func_light <- function(params) {
    # Everything it needs will come from its environment:
    # change_stats, in_mark, diffs, new_nodes, old_nodes, time, max_node_time, degenerate_edges
    
    if (degenerate_edges) {
      return(0)
    }
    
    eta    <- as.vector(change_stats %*% params$CS_params)
    p_base <- stats::plogis(eta)
    
    # same decay factor as direct
    p <- p_base * exp(-params$beta_edges * diffs)
    
    # 2. SAFETY CLAMP
    # Ensure p is never exactly 0 or 1. 
    # This prevents log(0) and log(1-1) errors.
    epsilon <- 1e-10
    p[p > (1 - epsilon)] <- 1 - epsilon
    p[p < epsilon] <- epsilon
    
    if (anyNA(p)) return(NA_real_)
    log_edge_part <- sum(log(p[in_mark])) + sum(log1p(-p[!in_mark]))
    node_dens <- if (!is.null(max_node_time) && time > max_node_time) {
      0
    } else {
      log(stats::dpois(new_nodes - old_nodes, params$node_lambda))
    }
    log_edge_part + node_dens
  }
  
  # Decide if you’re in the same degenerate branch as the direct computation
  degenerate_edges <- is.null(probs) || length(probs) == 1L
  
  # Now *force* a tiny environment (no local needed)
  environment(log_density_func_light) <- list2env(
    list(
      change_stats     = change_stats,
      in_mark          = in_mark,
      diffs            = if (exists("diffs", inherits = FALSE)) diffs else numeric(0),
      new_nodes        = new_nodes,
      old_nodes        = old_nodes,
      time             = time,
      max_node_time    = max_node_time,
      degenerate_edges = degenerate_edges
    ),
    parent = baseenv()
  )
  
  density_func_light <- function(params){
    log_density <- log_density_func_light(params)
    exp(log_density)
  }
  
  environment(density_func_light) <- list2env(
    list(
      change_stats = change_stats,
      in_mark      = in_mark,
      node_dens    = node_dens,
      plogis = stats::plogis
    ),
    parent = baseenv()
  )

  # =============
  # generate mark
  # =============
  if(generate_mark){
    # use latest mark as baseline:
    last_net <- mark
    # use function sample a new mark
    # since we have poisson number of nodes added  we need to redo the probabilities
    if(!is.null(last_net) && (last_net %n% 'n') >= 1){
      mark_sample <- last_net
      old_nodes <- last_net %n% 'n'
      
      if(mark_sample %n% 'n' < 4){
        old_new_net <- mark_sample
        new_nodes <- 4 - (mark_sample %n% 'n')
      }else{
        if(time > max_node_time){
          new_nodes <- 0
        }else{
          new_nodes <- rpois(1,params$node_lambda)
        }
      }
      mark_sample <- network::add.vertices(mark_sample,new_nodes)
      # if(mark_sample %n% 'n' > 4){
      #   browser()
      # }
      set.vertex.attribute(mark_sample,"time",c((last_net %v% 'time'),rep(time,new_nodes)))
      new_size <- mark_sample %n% 'n'
      
      # get new poss edges
      poss_tails <- (old_nodes - truncation) : (new_size)
      poss_tails <- poss_tails[poss_tails>0]
      poss_heads <- (new_size - truncation-1):new_size
      poss_heads <- poss_heads[poss_heads>0]
      poss_edges <- expand.grid(poss_tails,poss_heads)
      poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
      tails <- poss_edges[,1]
      heads <- poss_edges[,2]

      # only consider edges that are not in the old net
      if(!is.null(last_net)){
        in_old_net <- sapply(seq_along(heads), function(i) {
          length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
        })
        tails <- tails[!in_old_net]
        heads <- heads[!in_old_net]
      }


      delete.vertex.attribute(mark_sample,'na')
      model <- createCppModel(as.formula(paste("mark_sample ~ ",formula_RHS)))
      model$calculate()
      
      change_stats <- model$computeChangeStats(tails, heads)
      probs <- apply(change_stats, 1, function(c){
        1/(1+exp(-sum(c*params$CS_params)))
      })
      
      # reset to when we did not add more edges
      #mark_sample <- old_new_net
      # logistic regression on change stats:
      if(length(change_stats) == 0){
      stop("these parameters result ixn full networks - you probably don't want this")
      }

      if(mark_decay == 'activity'){
        node_times <- get_latest_times(mark_sample)
      }
      if(mark_decay == 'node_entrance'){
        node_times <- mark_sample %v% 'time'
      }
      diffs <- sapply(seq_along(tails), function(i) {
        node_times[tails[i]] - node_times[heads[i]]
      })
      # factor <- params$eta + (1-params$eta)*exp(-params$beta_edges*(diffs))
      factor <- exp(-params$beta_edges*(diffs))
      probs <- factor * probs

      add <- runif(length(probs)) < probs
      add.edges(mark_sample,
                heads[add],
                tails[add]
      )
      set.edge.attribute(mark_sample,"time",c(mark_sample %e% 'time',rep(time,sum(add))))
      mark_sample_density = prod(probs[add])*prod(1-probs[!add])*stats::dpois(new_nodes-old_nodes,params$node_lambda)
      log_mark_sample_density <- sum(log(probs[add])) +
                                 sum(log(1-probs[!add])) +
                                 log(stats::dpois(new_nodes-old_nodes,params$node_lambda))
      }else{
        if(is.null(last_net)){
          mark_sample <- network::network(matrix(1),directed = F)
          set.vertex.attribute(mark_sample,"time",time)
        }else{
          mark_sample <- last_net
        }
        times <- mark_sample %v% 'time'
        mark_sample <- network::add.vertices(mark_sample,1)
        set.vertex.attribute(mark_sample,
                             "time",
                             c(times,time))
        mark_sample_density <- 1
        log_mark_sample_density <- 0
      }
    }else{
      mark_sample <- new_net
      mark_sample_density <- 1
      log_mark_sample_density <- 0
    }

  return(list(
    # density of provided marks
    mark_density = mark_density,
    log_mark_density = log_mark_density,
    density_func = density_func_light,
    log_density_func = log_density_func_light,
    edge_probs = probs,
    mark_grad = mark_grad,
    decay_grad = decay_grad,
    # mark_sample
    mark_sample = mark_sample,
    mark_sample_density = mark_sample_density,
    log_mark_sample_density = log_mark_sample_density
  ))
}
