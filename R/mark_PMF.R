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
                        new_edge_hash = NULL
){
  if(length(mark_filtration$times) !=0){
    times <- mark_filtration$times
    last_net <- mark_filtration$marks[[length(mark_filtration$marks)]]
    new_nodes <- 1
    old_nodes <- last_net %n% 'n'
    new_net <- last_net
    current_times <- get.vertex.attribute(new_net,"time")
    add.vertices(new_net,nv=new_nodes)
    set.vertex.attribute(new_net,"time",c(current_times,time))
  }else{
    times  <- c()
    last_net <- NULL
    new_net <- network::network(matrix(1),directed = F)
    set.vertex.attribute(new_net,"time",time)
    old_nodes <- 0
    new_nodes <- 1
  }

  # get the possible edges
  tails <- rep((old_nodes+1):(old_nodes + new_nodes),times = old_nodes)
  heads <- rep(1:old_nodes,each = new_nodes)
  keep <- tails != heads
  heads <- heads[keep]
  tails <- tails[keep]

  # get the edge probabilities
  if(!is.null(last_net)){
    degs <- c(degree(last_net,gmode = 'graph'),rep(0,new_nodes))
    times <- c(last_net %v% 'time',rep(time,new_nodes))
    if(length(degs) != length(times)){
      stop("length of degrees and times do not match")
    }
    degs <- degs*exp(-params$beta_BA_edges*(time - times))
    total_deg <- sum(degs)
    if(total_deg == 0){
      probs <- rep(1,length(heads))
      in_mark <- rep(TRUE,length(heads))
    }else{
      probs <- degs[heads] / total_deg
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
    }else{
      log_mark_density <- sum(log(probs[in_mark])) + sum(log(1-probs[!in_mark]))
      mark_density <- exp(log_mark_density)
    }
  }else{
    mark_density <- 1
    mark_density_normalized <- NULL
    log_mark_density <- 0
  }

  if(generate_mark & !is.null(probs)){
    # use function sample a new mark
    if(!is.null(last_net)){
      mark_sample <- last_net
      mark_sample <- network::add.vertices(mark_sample,new_nodes)
      set.vertex.attribute(mark_sample,"time",c((last_net %v% 'time'),rep(time,new_nodes)))
    }else{
      mark_sample <- network::network(matrix(1),directed = F)
      set.vertex.attribute(mark_sample,"time",rep(time,new_nodes))
    }
    add <- runif(length(probs)) < probs
    add.edges(mark_sample,
              heads[add],
              tails[add]
    )
    mark_sample_density = prod(probs[add])*prod(1-probs[!add])
    log_mark_sample_density <- sum(log(probs[add])) + sum(log(1-probs[!add]))
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
    # mark_sample
    mark_sample = mark_sample,
    mark_sample_density = mark_sample_density,
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
                        new_edge_hash = NULL,
                        formula_RHS,
                        grad = FALSE
){
  if(length(mark_filtration$times) !=0){
    times <- mark_filtration$times
    last_net <- mark_filtration$marks[[length(mark_filtration$marks)]]
    new_nodes <- 1
    old_nodes <- last_net %n% 'n'
    new_net <- last_net
    add.vertices(new_net,nv=new_nodes)
  }else{
    times  <- c()
    last_net <- NULL
    new_net <- network::network(matrix(1),directed = F)
    set.vertex.attribute(new_net,"time",time)
    old_nodes <- 0
    new_nodes <- 1
  }

  # get the possible edges
  tails <- rep((old_nodes+1):(old_nodes + new_nodes),times = old_nodes)
  heads <- rep(1:old_nodes,each = new_nodes)
  keep <- tails != heads
  heads <- heads[keep]
  tails <- tails[keep]

  # get the edge probabilities
  if(!is.null(last_net)){
    if(last_net %n% 'n' > 2){
      # if new net has less than 4 nodes add some:
      if(new_net %n% 'n' < 4){
        old_new_net <- new_net
        new_net <- network::add.vertices(new_net,4 - (new_net %n% 'n'))
      }else{
        old_new_net <- new_net
      }
      model <- ernm(as.formula(paste("new_net ~ ",formula_RHS)),
                    tapered = FALSE,
                    maxIter = 3,
                    mcmcBurnIn = 100,
                    mcmcInterval = 10,
                    mcmcSampleSize = 100,
                    verbose = 0)
      model <- model$m$sampler$getModel()
      model$setNetwork(ernm::as.BinaryNet(new_net))
      change_stats <- lapply(1:length(tails),FUN=function(i){
        #reset:
        model$calculate()
        stat <- model$statistics()
        # update
        model$dyadUpdate(tails[i],heads[i])
        new_stat <- model$statistics()
        return(new_stat - stat)
      })
      new_net <- old_new_net
      # logistic regression on change stats:
      probs <- 1/(1+exp(-sapply(change_stats,function(c){sum(c*params$CS_params)})))
      times <- last_net %v% 'time'
      probs <- probs*exp(-params$beta_edges*(time - times))
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
      log_mark_density <- sum(log(probs[in_mark])) + sum(log(1-probs[!in_mark]))
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

  if(generate_mark & !is.null(probs)){
    # use function sample a new mark
    if(!is.null(last_net)){
      mark_sample <- last_net
      mark_sample <- network::add.vertices(mark_sample,new_nodes)
      set.vertex.attribute(mark_sample,"time",c((last_net %v% 'time'),rep(time,new_nodes)))
    }else{
      mark_sample <- network::network(matrix(1),directed = F)
      set.vertex.attribute(mark_sample,"time",rep(time,new_nodes))
    }
    add <- runif(length(probs)) < probs
    add.edges(mark_sample,
              heads[add],
              tails[add]
    )
    mark_sample_density = prod(probs[add])*prod(1-probs[!add])
    log_mark_sample_density <- sum(log(probs[add])) + sum(log(1-probs[!add]))
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
