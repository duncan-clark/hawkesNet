#' Mark probability mass function for the network generation process
#'
#' Calculates the probability mass function (PMF) for marks (network structures) at a given time, supporting multiple network growth models (Barabasi–Albert, Change Statistic Hawkes, and BA-bipartite). 
#'
#' @references
#' Barabási, A.-L. & Albert, R. (1999). Emergence of scaling in random networks. *Science*, 286, 509–512. \doi{10.1126/science.286.5439.509}
#' @param time Numeric. The time at which to evaluate the PMF.
#' @param params Named list. Model parameter values required for the chosen \code{type}:
#'   - For \code{type = "BA"}: \code{beta_edges} (numeric).
#'   - For \code{type = "CS"}: \code{beta_edges} (numeric), \code{node_lambda} (numeric), \code{CS_params} (numeric vector of coefficients for \code{formula_RHS}).
#'   - For \code{type = "BA-bip"}: similar to BA, plus any relevant bipartite parameters.
#' @param mark_filtration Network or compatible object. The network history (filtration) up to the current time.
#' @param type Character. Model type: one of \code{"BA"}, \code{"CS"}, or \code{"BA-bip"}. Default: \code{"BA"}.
#' @param mark Network or NULL. The current mark/network structure. If \code{NULL}, derived from \code{mark_filtration}.
#' @param generate_mark Logical. If \code{TRUE}, generates a new mark/sample; otherwise computes density for the supplied mark. Default: \code{FALSE}.
#' @param generate_density Logical. If \code{TRUE}, computes the density for the provided mark. Default: \code{TRUE}.
#' @param grad Logical. If \code{TRUE}, also computes gradients of the mark density. Default: \code{FALSE}.
#' @param new_edge_hash Hash or NULL. Optional hashed edge list for fast lookup. Default: \code{NULL}.
#' @param truncation Integer or NULL. For \code{"CS"} models, controls which edges are considered (e.g., \code{1} means only new-to-old). Default: \code{NULL}.
#' @param formula_RHS Character or formula. For \code{"CS"} models, specifies the right-hand-side for change statistics calculation. Default: \code{NULL}.
#' @param mark_decay Character or NULL. How edge decay is modeled, e.g. \code{"node_entrance"}, \code{"activity"}, etc. Default: \code{NULL}.
#' @param model Object or NULL. Preconstructed model object (for efficiency); if \code{NULL}, will be built internally. Default: \code{NULL}.
#' @param max_node_time Numeric. The last time at which a node can enter the network (CS model). Default: \code{10}.
#' @param ... Additional arguments, passed to model-specific PMF functions.
#' @return Named list containing:
#'   \item{mark_density}{Numeric. Density of the provided or generated mark.}
#'   \item{log_mark_density}{Numeric. Log-density of the mark.}
#'   \item{edge_probs}{Numeric vector. Probabilities for each possible edge.}
#'   \item{mark_grad}{Numeric vector. Gradient of the mark density (if \code{grad = TRUE}).}
#'   \item{decay_grad}{Numeric vector. Gradient with respect to decay (if \code{grad = TRUE}).}
#'   \item{mark_sample}{Network. The sampled mark/network object.}
#'   \item{mark_sample_density}{Numeric. Density of the sampled mark.}
#'   \item{log_mark_sample_density}{Numeric. Log-density of the sampled mark.}
#'
#' @details Computes the mark PMF, \eqn{q(m\vert t,\mathcal{H}_{t})} (see \code{\link{cond_intensity}}).
#' Currently three options: \code{type = "BA"}, \code{type = "CS"}, and \code{type = "BA-bip"}. 
#'
#' For \code{type = "BA"} the Barabasi Albert (BA) preferential attachment model is used where the mark distribution is defined as
#' \deqn{
#'   q(m \mid t, \mathcal{H}_t) =
#'   \prod_{i=1}^{N_{t-}} \left(p_i^{BA}\right)^{e_i} \cdot
#'   \left(1 - p_i^{BA}\right)^{1 - e_i}.
#' }
#' Here, the attachment probability, \eqn{p_i^{BA}}, is defined as
#' \deqn{
#'   p_i^{BA} = \frac{\delta_i}{\sum_{k=1}^{N} \delta_k}
#' }
#' where \eqn{\delta_{i}^{t} = \exp(\tau \cdot (t - t_i)) \cdot d_{i}^{t}}, 
#' and \eqn{d_{i}^{t}} is the sna::degree of node \eqn{i} just before time \eqn{t}.
#' 
#' For \code{type = "CS"} the change statistic (CS) model is used where the mark distribution is defined (similar to above) as
#' \deqn{
#'   q(m \mid t, \mathcal{H}_t) =
#'   \prod_{i=1}^{N_{t-}} \left(p_i^{CS}\right)^{e_i} \cdot
#'   \left(1 - p_i^{CS}\right)^{1 - e_i}.
#' }
#' where the attachment probability, \eqn{p_i^{CS}}, is defined as
#' \deqn{
#' p_i^{CS} = \left(\nu + \exp(\tau \cdot(t - t_i))\right)\cdot\frac{1}{1 + \exp(-\theta^{\top} \cdot C_{i,N_t})}
#' }
#' @examples
#' \dontrun{
#' if(interactive()){
#'  data(net, package = "hawkesGrowthNet")
#' time <- get_times(net)$times
#' ## BA
#' params_ba <-  list( beta_edges = 0.1)
#' mark_filtration <-  filtration_to_net(net,10)
#' pmf_ba <- PMF_mark(time[10],  params_ba, mark_filtration)
#' ## CS
#' devtools::install_github("duncan-clark/ernm", ref = "R_change_stats")
#' require(ernm)
#' params_cs <-  list(beta_edges = 0.1,node_lambda = 1,CS_params =  c(-10,0,0,0))
#' pmf_cs <- PMF_mark(time = time[10],  params = params_cs,
#' mark_filtration = mark_filtration, type = "CS",  truncation = 1,
#' formula_RHS = "edges + triangles + star(c(2,3))",
#' max_node_time = 1)
#'  }
#' }
#' @seealso
#' \code{\link{cond_intensity}}
#' \code{\link{PMF_mark_BA}}
#' \code{\link{PMF_mark_CS}}
#' \code{\link[network]{network}}, \code{\link[network]{add.vertices}}
#' \code{\link[ernm]{as.BinaryNet}}
#' @rdname PMF_mark
#' @export
PMF_mark <- function(time,
                     params,
                     mark_filtration,
                     type = c("BA", "CS"),
                     mark = NULL,
                     generate_mark = FALSE,
                     generate_density = TRUE,
                     grad = FALSE,
                     new_edge_hash = NULL,
                     truncation = NULL,
                     formula_RHS,
                     mark_decay = 'node_entrance',
                     model = NULL,
                     max_node_time = 10,
                     ...){
    type <- type[1]
    if (!(type %in% c("BA", "CS", "BA-bip"))) {
        stop("type can only be one of `BA` for Barabasi–Albert, `CS` for change statistic Hawkes, or `BA-bip` for bipartite Barabási–Albert.")
    }
    if(type == "BA"){
        pmf <- PMF_mark_BA(time, params, mark_filtration, mark,
                           generate_mark, generate_density, grad,
                           new_edge_hash, truncation, ...)
    }else{
        if(type == "BA-bip"){
            pmf <- PMF_mark_BA_bipartite(time, params,  mark_filtration,
                                      mark, generate_mark, generate_density,
                                      grad, new_edge_hash, truncation, ...)
        }else{
            if(type == "CS"){
                pmf <- PMF_mark_CS(time, params, mark_filtration,
                                   mark, generate_mark, generate_density,
                                   grad,  new_edge_hash, truncation, formula_RHS,
                                   mark_decay, model, max_node_time, ...)
            }
        }
    }
    return(pmf)
}


#' Internal function to prepare for mark PMFs
#' @inheritParams PMF_mark
#' @noRd
mark_setup <- function(mark, mark_filtration, time){
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
    ## get the possible edges for the given truncation:
    ## if no nodes have been added then :
    poss_tails <- seq.int((old_nodes + 1), new_nodes)
    poss_tails <- poss_tails[poss_tails>0]
    poss_heads <- 1:old_nodes
    poss_heads <- poss_heads[poss_heads>0]
    poss_edges <- expand.grid(poss_tails,poss_heads)
    poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
    tails <- poss_edges[,1]
    heads <- poss_edges[,2]
    ## only consider edges that were not already in the old network
    if(!is.null(last_net)){
        in_old_net <- sapply(1:length(heads),function(i){
            length(get.edgeIDs(last_net, heads[i],tails[i])) !=0
        })
        tails <- tails[!in_old_net]
        heads <- heads[!in_old_net]
    }
    return(list(mark = mark, new_net = new_net, last_net = last_net,
                poss_tails = poss_tails, poss_heads = poss_heads,
                poss_edges = poss_edges,
                tails = tails, heads = heads, new_nodes = new_nodes,
                old_nodes = old_nodes))
}
#' Internal function to prepare for bipartite mark PMFs
#' @inheritParams PMF_mark
#' @noRd
mark_setup_bipartite <- function(mark = NULL, mark_filtration, time){
    if(is.null(mark)){
        mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
    }
    last_net <- filtration_to_net(mark_filtration, time, equals = FALSE)
    new_net <- last_net
    if(is.null(last_net) || (last_net %n% "n") == 0){
        last_net <- NULL
        new_net <- network::network(matrix(0, 1, 1), directed = FALSE, bipartite = 0)
        set.vertex.attribute(new_net, "time", time)
        set.vertex.attribute(new_net, "role", "event")  
        old_nodes <- 0
        new_nodes <- 1
    } else {
        new_nodes <- mark %n% "n"
        old_nodes <- last_net %n% "n"
        network::add.vertices(new_net, nv = new_nodes - old_nodes)
        set.vertex.attribute(new_net, "time",
                             c(get.vertex.attribute(last_net, "time"),
                               rep(time, new_nodes - old_nodes)))
        roles <- get.vertex.attribute(last_net, "role")
        if(is.null(roles)) roles <- rep("perp", old_nodes)
        new_roles <- c(roles, rep("event", new_nodes - old_nodes))
        set.vertex.attribute(new_net, "role", new_roles)
    }
    if(is.null(last_net)){
        last_net <- network::network(matrix(0, 0, 0), directed = FALSE, bipartite = 0)
        set.vertex.attribute(last_net, "time", numeric(0))
        set.vertex.attribute(last_net, "role", character(0))
    }
    perp_nodes <- which(get.vertex.attribute(new_net, "role") == "perp")
    event_nodes <- which(get.vertex.attribute(new_net, "role") == "event")
    poss_tails <- seq.int(from = old_nodes + 1, to = new_nodes)
    poss_tails <- poss_tails[poss_tails > 0]
    poss_heads <- perp_nodes
    poss_heads <- poss_heads[poss_heads > 0]
    poss_edges <- expand.grid(poss_tails, poss_heads)
    colnames(poss_edges) <- c("tail", "head")
    tails <- poss_edges[, "tail"]
    heads <- poss_edges[, "head"]
    if(!is.null(last_net) && length(heads) > 0){
        in_old_net <- sapply(1:length(heads), function(i){
            length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
        })
        tails <- tails[!in_old_net]
        heads <- heads[!in_old_net]
    }
    return(list(mark = mark,
                new_net = new_net,
                last_net = last_net,
                poss_tails = poss_tails,
                poss_heads = poss_heads,
                poss_edges = poss_edges,
                tails = tails,
                heads = heads,
                new_nodes = new_nodes,
                old_nodes = old_nodes))
}
#' Function for Barabási–Albert (BA) probability mass function
#' @rdname PMF_mark
#' @export
PMF_mark_BA <- function(time,
                        params,
                        mark_filtration,
                        mark,
                        generate_mark,
                        generate_density,
                        grad,
                        new_edge_hash,
                        truncation, ...){
    ## shared setup
    setup <- mark_setup(mark, mark_filtration, time)
    mark <- setup$mark
    new_net <- setup$new_net
    last_net <- setup$last_net
    poss_tails <- setup$poss_tails
    poss_heads <- setup$poss_heads
    poss_edges <- setup$poss_edges
    tails <- setup$tails
    heads <- setup$heads
    
    if(!is.null(last_net) && (last_net %n% 'n' > 2)){
        times <- get.vertex.attribute(last_net, "time")
        degs <- sna::degree(last_net) * exp(-params$beta_edges * (time - times))
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
        last_net <- mark
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
            ## only consider edges that are not in the old net
            if(!is.null(last_net)){
                in_old_net <- sapply(1:length(heads),function(i){
                    length(get.edgeIDs(last_net, heads[i],tails[i])) !=0
                })
                tails <- tails[!in_old_net]
                heads <- heads[!in_old_net]
            }
            times <- get.vertex.attribute(last_net, "time")
            degs <- sna::degree(last_net) * exp(-params$beta_edges * (time - times))
            total_deg <- sum(degs, na.rm = TRUE)
            
            if(total_deg == 0){
                probs <- rep(1, length(heads))
            } else {
                probs <- degs[heads] / total_deg
            }
            
            if(any(is.na(probs))){
                browser()
            }
            add <- runif(length(probs)) <  probs
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
#' Function for BA-style bipartite probablibity mass function
#' @rdname PMF_mark
#' @export
PMF_mark_BA_bipartite <- function(time,
                                  params,
                                  mark_filtration,
                                  mark = NULL,
                                  generate_mark = TRUE,
                                  generate_density = TRUE,
                                  grad = NULL,
                                  new_edge_hash = NULL,
                                  truncation = NULL, ...) {
  
  ## --- Shared setup ---
  setup <- mark_setup_bipartite(mark, mark_filtration, time)
  mark <- setup$mark
  new_net <- setup$new_net
  last_net <- setup$last_net
  poss_tails <- setup$poss_tails
  poss_heads <- setup$poss_heads
  poss_edges <- setup$poss_edges
  tails <- setup$tails
  heads <- setup$heads
  new_nodes <- setup$new_nodes
  old_nodes <- setup$old_nodes
  perp_nodes <- which(get.vertex.attribute(new_net, "role") == "perp")
    if(length(perp_nodes) > 0){
        times <- get.vertex.attribute(new_net, "time")
        c_offset <- 0.01 ## small sna::degree fix
        degs <- c_offset + sna::degree(new_net)[perp_nodes] * exp(-params$beta_edges * (time - times[perp_nodes]))
        total_deg <- sum(degs, na.rm = TRUE)
        if(total_deg == 0){
            probs <- rep(1, length(perp_nodes))
        } else {
            probs <- degs / total_deg
        }
        heads <- perp_nodes  
    } else {
        probs <- numeric(0)
        heads <- integer(0)
    }
  if(!is.null(mark) && length(probs) != 0){
    if(is.null(new_edge_hash)){
      in_mark <- sapply(1:length(heads), function(i){
        length(get.edgeIDs(mark, heads[i], poss_tails)) != 0
      })
    } else {
      in_mark <- has_edge(heads, poss_tails, new_edge_hash)
    }
    log_mark_density <- sum(log(probs[in_mark])) + sum(log(1 - probs[!in_mark]))
    mark_density <- exp(log_mark_density)
  } else {
    log_mark_density <- 0
    mark_density <- 1
  }
  if(generate_mark){
    mark_sample <- new_net
    new_event_id <- as.integer(mark_sample %n% "n") + 1
    mark_sample <- network::add.vertices(mark_sample, 1)
    set.vertex.attribute(mark_sample, "time", c(get.vertex.attribute(mark_sample, "time"), time))
    set.vertex.attribute(mark_sample, "role", c(get.vertex.attribute(mark_sample, "role"), "event"))
    perp_nodes <- which(get.vertex.attribute(mark_sample, "role") == "perp")
    if(length(perp_nodes) > 0){
      times <- get.vertex.attribute(mark_sample, "time")
      degs <- sna::degree(mark_sample)[perp_nodes] * exp(-params$beta_edges * (time - times[perp_nodes]))
      total_deg <- sum(degs)
      if(total_deg == 0){
        probs <- rep(1, length(perp_nodes))
      } else {
        probs <- degs / total_deg
      }
      add <- runif(length(probs)) < probs
      network::add.edges(mark_sample, perp_nodes[add], rep(new_event_id, sum(add)))
      log_mark_sample_density <- sum(log(probs[add])) + sum(log(1 - probs[!add]))
      mark_sample_density <- exp(log_mark_sample_density)
    } else {
      log_mark_sample_density <- 0
      mark_sample_density <- 1
    }
    K <- rpois(1, lambda = params$lambda_new)
    if(K > 0){
      new_perp_ids <- (mark_sample %n% "n") + seq_len(K)
      mark_sample <- network::add.vertices(mark_sample, K)
      set.vertex.attribute(mark_sample, "time",
                           c(get.vertex.attribute(mark_sample, "time"), rep(time, K)))
      set.vertex.attribute(mark_sample, "role",
                           c(get.vertex.attribute(mark_sample, "role"), rep("perp", K)))
      network::add.edges(mark_sample, new_perp_ids, rep(new_event_id, K))
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
    mark_sample_density = mark_sample_density,
    log_mark_sample_density = log_mark_sample_density
  ))
}

#' Function for change statistic (CS) HawkesNet probability mass function
#' @rdname PMF_mark
#' @export
PMF_mark_CS <- function(time,
                        params,
                        mark_filtration,
                        mark,
                        generate_mark,
                        generate_density,
                        grad,
                        new_edge_hash,
                        truncation,
                        formula_RHS,
                        mark_decay,
                        model,
                        max_node_time, ...){
 
    ## shared setup
    setup <- mark_setup(mark, mark_filtration, time)
    mark <- setup$mark
    new_net <- setup$new_net
    last_net <- setup$last_net
    poss_tails <- setup$poss_tails
    poss_heads <- setup$poss_heads
    poss_edges <- setup$poss_edges
    tails <- setup$tails
    heads <- setup$heads
    ## CS specific
    new_nodes <- setup$new_nodes
    old_nodes <-  setup$old_nodes
    ## =============
    ## mark density
    ## =============
    if(!is.null(last_net) & generate_density){
        if(last_net %n% 'n' > 2){
            ## if new net has less than 4 nodes add some:
            if(new_net %n% 'n' < 4){
                old_new_net <- new_net
                new_net <- network::add.vertices(new_net,4 - (new_net %n% 'n'))
            }else{
                old_new_net <- new_net
            }

            if(max(tails)>new_net %n% 'n'){
                stop("accidently adding a edge into the network that doesn't have that node yet")
            }
            ## delete NAs to prevent C++ using them
            delete.vertex.attribute(new_net,'na')
            if(is.null(model)){
                model <- ernm::createCppModel(as.formula(paste("new_net ~ ",formula_RHS)))
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
            ## use either node times or last node activity:
            if(mark_decay == 'activity'){
                node_times <- get_latest_times(new_net)
            }
            if(mark_decay == 'node_entrance'){
                node_times <- new_net %v% 'time'
            }
            diffs <- time - node_times[heads]
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
            if(time >max_node_time){
                node_dens <- 0
            }else{
                node_dens <- log(dpois(new_nodes-old_nodes,params$node_lambda))
            }
            log_mark_density <- sum(log(probs[in_mark])) + sum(log(1-probs[!in_mark])) + node_dens
            mark_density <- exp(log_mark_density)

            if(grad){
                ## get the mark grad:
                probs_grads <- lapply(change_stats,function(c){
                    (-c)*exp(-sum(c*params$CS_params))/((1+exp(-sum(c*params$CS_params)))^2)}
                    )
                e <- in_mark*1
                derivs <- mapply(probs_grads,probs,in_mark,FUN = function(dp,p,e){
                    tmp <- e*p^(e-1)*(1-p)^(1-e) - (1-e)*p^e*(1-p)^(-e)
                    tmp <- dp*tmp
                    tmp <- sign(tmp)*exp(log(abs(tmp)) - log(p))
                    return(tmp)
                },SIMPLIFY = FALSE)
                ## get devided by the right prob:
                derivs <- do.call(rbind,derivs)
                mark_grad <- colSums(derivs)

                ## get the mark grad:
                decay_grads <- -(time-times)*exp(-params$beta_edges*(time - times))
                e <- in_mark*1
                derivs <- mapply(decay_grads,probs,in_mark,FUN = function(dp,p,e){
                    tmp <- e*p^(e-1)*(1-p)^(1-e) - (1-e)*p^e*(1-p)^(-e)
                    tmp <- dp*tmp
                    tmp <- sign(tmp)*exp(log(abs(tmp)) - log(p))
                    return(tmp)
                },SIMPLIFY = FALSE)
                ## get devided by the right prob:
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

    ## =============
    ## generate mark
    ## =============
    if(generate_mark){
        ## use latest mark as baseline:
        last_net <- mark
        ## use function sample a new mark
        ## since we have poisson number of nodes added  we need to redo the probabilities
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
            set.vertex.attribute(mark_sample,"time",c((last_net %v% 'time'),rep(time,new_nodes)))
            new_size <- mark_sample %n% 'n'
            
            ## get new poss edges
            poss_tails <- (old_nodes - truncation) : (new_size)
            poss_tails <- poss_tails[poss_tails>0]
            poss_heads <- (new_size - truncation-1):new_size
            poss_heads <- poss_heads[poss_heads>0]
            poss_edges <- expand.grid(poss_tails,poss_heads)
            poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2],]
            tails <- poss_edges[,1]
            heads <- poss_edges[,2]

            ## only consider edges that are not in the old net
            if(!is.null(last_net)){
                in_old_net <- sapply(1:length(heads),function(i){
                    length(get.edgeIDs(last_net, heads[i],tails[i])) !=0
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
            diffs <- sapply(1:length(tails),function(i){
                node_times[tails[i]] - node_times[heads[i]]
            })
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
        ## density of provided marks
        mark_density = mark_density,
        log_mark_density = log_mark_density,
        edge_probs = probs,
        mark_grad = mark_grad,
        decay_grad = decay_grad,
        ## mark_sample
        mark_sample = mark_sample,
        mark_sample_density = mark_sample_density,
        log_mark_sample_density = log_mark_sample_density
    ))
}
