# roxgen documentation
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param time PARAM_DESCRIPTION
#' @param params PARAM_DESCRIPTION
#' @param mark_filtration PARAM_DESCRIPTION
#' @param mark PARAM_DESCRIPTION, Default: NULL
#' @param generate_mark PARAM_DESCRIPTION, Default: FALSE
#' @param new_edge_hash PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[network]{network}}, \code{\link[network]{add.vertices}}
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
    in_old_net <- sapply(1:length(heads),function(i){
      length(get.edgeIDs(last_net, heads[i],tails[i])) !=0
    })
    tails <- tails[!in_old_net]
    heads <- heads[!in_old_net]
  }
  
  if(!is.null(last_net) && (last_net %n% 'n' > 2)){
    times <- get.vertex.attribute(last_net, "time")
    degs <- degree(last_net) * exp(-params$beta_edges * (time - times))
    total_deg <- sum(degs)
    
    if(total_deg == 0){
      probs <- rep(1, length(heads))
    } else {
      probs <- degs[heads] / total_deg
    }
  } else {
    probs <- rep(1, length(heads))
  }
  
  if(!is.null(mark) && length(probs) !=0 && (last_net %n% 'n' > 2)){
    if(is.null(new_edge_hash)){
      in_mark <- sapply(1:length(heads),function(i){
        length(get.edgeIDs(mark, heads[i], tails[i])) != 0
      })
    } else {
      in_mark <- has_edge(heads, tails, new_edge_hash)
    }
    
    log_mark_density <- sum(log(probs[in_mark])) + sum(log(1 - probs[!in_mark]))
    mark_density <- exp(log_mark_density)
  } else {
    log_mark_density <- 0
    mark_density <- 1
  }
  
  if(generate_mark){
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
        in_old_net <- sapply(1:length(heads),function(i){
          length(get.edgeIDs(last_net, heads[i],tails[i])) !=0
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
    edge_probs = probs,
    mark_sample = mark_sample,
    mark_sample_density = exp(log_mark_sample_density),
    log_mark_sample_density = log_mark_sample_density
  ))
}

# roxgen documentation
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param time PARAM_DESCRIPTION
#' @param params PARAM_DESCRIPTION
#' @param mark_filtration PARAM_DESCRIPTION
#' @param mark PARAM_DESCRIPTION, Default: NULL
#' @param generate_mark PARAM_DESCRIPTION, Default: FALSE
#' @param new_edge_hash PARAM_DESCRIPTION, Default: NULL
#' @param formula_RHS PARAM_DESCRIPTION
#' @param grad PARAM_DESCRIPTION, Default: FALSE
#' @param truncation if truncation = 1, only consider edges from new nodes to old nodes,
#' truncation = k considers edges from k time steps before the new nodes to the old nodes:
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[network]{network}}, \code{\link[network]{add.vertices}}
#'  \code{\link[ernm]{as.BinaryNet}}
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
                        ...
){
  # if we are not starting from nothing:
  last_net <- filtration_to_net(mark_filtration,time,equals = FALSE)
  # set equal to last net - then add things
  new_net <- last_net
  if(is.null(mark)){
    mark <- filtration_to_net(mark_filtration,time,equals = TRUE)
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
  if(!is.null(last_net)){
    in_old_net <- sapply(1:length(heads),function(i){
      length(get.edgeIDs(last_net, heads[i],tails[i])) !=0
    })
    tails <- tails[!in_old_net]
    heads <- heads[!in_old_net]
  }

  # =============
  # mark density
  # =============
  if(!is.null(last_net) & generate_density){
    if(last_net %n% 'n' > 2){
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
      # print("Max tail")
      # print(max(tails))
      # print("max head")
      # print(max(heads))
      # print("new network")
      # print(summary(new_net,print.adj = F))
      # print("vertex names")
      # print(new_net %v% 'vertex.names')
      # print(as.factor(new_net %v% 'vertex.names'))
      # print(summary(new_net %v% 'time'))
      # print(summary(new_net %e% 'time'))
      # print("making model")
      
      # delete NAs to prevent C++ using them
      delete.vertex.attribute(new_net,'na')
      model <- createCppModel(as.formula(paste("new_net ~ ",formula_RHS)))
      # print("model made")
      # model$setNetwork(ernm::as.BinaryNet(new_net))
      new_net <- old_new_net
      model$calculate()
      stat <- model$statistics()
      # print("Model statistics")
      # print(stat)
      # print("old new network - used in change stats")
      # print(summary(old_new_net,print.adj = F))
      # print("doing change stats")
      change_stats <- lapply(1:length(tails),FUN=function(i){
        # update - note no need to update just need to take away  old stat
        old_stat <- model$statistics()
        model$dyadUpdate(tails[i],heads[i])
        new_stat <- model$statistics()
        return(new_stat - old_stat)
      })
      # print("done with change stats")
      # logistic regression on change stats:
      probs <- 1/(1+exp(-sapply(change_stats,function(c){sum(c*params$CS_params)})))
      # use either node times or last node activity:
      if(mark_decay == 'activity'){
        node_times <- get_latest_times(new_net)
      }
      if(mark_decay == 'node_entranace'){
        node_times <- new_net %v% 'time'
      }
      diffs <- time - node_times[heads]
      #factor <- params$eta + (1-params$eta)*exp(-params$beta_edges*(diffs))
      factor <- exp(-params$beta_edges*(diffs))
      probs <- probs * factor

    }else{
      times <- last_net %v% 'time'
      probs <- c(1)
    }
  }else{
    probs <- NULL
  }

  if(!is.null(mark) & !is.null(probs)){
    if(is.null(new_edge_hash)){
      in_mark <- sapply(1:length(heads),function(i){
        length(get.edgeIDs(mark, heads[i],tails[i])) !=0
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
      log_mark_density <- sum(log(probs[in_mark])) + sum(log(1-probs[!in_mark])) + log(dpois(new_nodes-old_nodes,params$node_lambda))
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

  # =============
  # generate mark
  # =============
  if(generate_mark){
    # use function sample a new mark
    # since we have poisson number of nodes added  we need to redo the probabilities
    if(!is.null(last_net) && (last_net %n% 'n') > 2){
      mark_sample <- last_net
      old_nodes <- last_net %n% 'n'
      new_nodes <- rpois(1,params$node_lambda)
      #browser()
      mark_sample <- network::add.vertices(mark_sample,new_nodes)
      set.vertex.attribute(mark_sample,"time",c((last_net %v% 'time'),rep(time,new_nodes)))
      new_nodes <- mark_sample %n% 'n'

      # get new poss edges
      poss_tails <- (old_nodes - truncation) : (new_nodes)
      poss_tails <- poss_tails[poss_tails>0]
      poss_heads <- (new_nodes - truncation-1):new_nodes
      poss_heads <- poss_heads[poss_heads>0]
      poss_edges <- expand.grid(poss_tails,poss_heads)
      poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
      tails <- poss_edges[,1]
      heads <- poss_edges[,2]

      # only consider edges that are not in the old net
      if(!is.null(last_net)){
        in_old_net <- sapply(1:length(heads),function(i){
          length(get.edgeIDs(last_net, heads[i],tails[i])) !=0
        })
        tails <- tails[!in_old_net]
        heads <- heads[!in_old_net]
      }

      # get the probs:
      if(mark_sample %n% 'n' < 4){
        old_new_net <- mark_sample
        mark_sample <- network::add.vertices(mark_sample,4 - (mark_sample %n% 'n'))
      }else{
        old_new_net <- mark_sample
      }
      delete.vertex.attribute(mark_sample,'na')
      model <- createCppModel(as.formula(paste("mark_sample ~ ",formula_RHS)))
      # model$setNetwork(ernm::as.BinaryNet(new_net))
      model$calculate()
      change_stats <- lapply(1:length(tails),FUN=function(i){
        # update - note no need to update just need to take away  old stat
        old_stat <- model$statistics()
        model$dyadUpdate(tails[i],heads[i])
        new_stat <- model$statistics()
        return(new_stat - old_stat)
      })
      # reset to when we did not add more edges
      mark_sample <- old_new_net
      # logistic regression on change stats:
      if(length(change_stats) == 0){
      stop("these parameters result ixn full networks - you probalby don't want this")
      }
      probs <- 1/(1+exp(-sapply(change_stats,function(c){sum(c*params$CS_params)})))
      if(mark_decay == 'activity'){
        node_times <- get_latest_times(mark_sample)
      }
      if(mark_decay == 'node_entranace'){
        node_times <- mark_sample %v% 'time'
      }
      diffs <- sapply(1:length(tails),function(i){
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
      mark_sample_density = prod(probs[add])*prod(1-probs[!add])*dpois(new_nodes-old_nodes,params$node_lambda)
      log_mark_sample_density <- sum(log(probs[add])) +
                                 sum(log(1-probs[!add])) +
                                 log(dpois(new_nodes-old_nodes,params$node_lambda))
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
    edge_probs = probs,
    mark_grad = mark_grad,
    decay_grad = decay_grad,
    # mark_sample
    mark_sample = mark_sample,
    mark_sample_density = mark_sample_density,
    log_mark_sample_density = log_mark_sample_density
  ))
}
