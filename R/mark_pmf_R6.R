
#' @details Base Class for markKernel - that adds network updates to the network
#' note that the mark space is dynamic so the kernel depends on the filtration
#' so as the filtration grows, the markKernel changes
#' we let the mark kernel be stateful - i.e. it stores the most up to date state - but also refelcts the dynamic mark space:
#' 
#' Base class for mark-update kernels on dynamic mark spaces
#' Shared plumbing; subclasses supply edge-probability logic.
MarkKernel <- R6::R6Class(
  "MarkKernel",
  public = list(
    params     = NULL,   # list of numeric parameters
    filtration = NULL,   # filtration object your helpers work with
    opts       = NULL,   # list(generate_mark, generate_density, truncation, max_node_time, ...)
    bipartite  = FALSE,  # hint for setup
    model      = NULL,   # optional compiled model / adapter (e.g., ernm)
    decay_fun  = NULL,   # function(node_times, time, params, heads_or_tails) -> factor
    initialize = function(params, filtration, opts = list(),
                          decay_fun = NULL, model = NULL, bipartite = FALSE) {
      self$params     <- params
      self$filtration <- filtration
      self$model      <- model
      if(is.null(decay_fun)){
        self$decay_fun  <- function(node_times, time, params, idx) {
          # default: exp(-beta * (time - node_time)) on heads
          exp(-params$beta_edges * (time - node_times[idx]))
        }
      }else{
        self$decay_fun  <- decay_fun
      }

      self$bipartite <- isTRUE(bipartite)
      defaults <- list(
        grad             = FALSE,
        truncation       = NULL,
        max_node_time    = 10,
        new_edge_hash    = NULL,
        formula_RHS      = NULL,
        mark_decay       = "node_entrance" # or "activity"
      )
      self$opts <- utils::modifyList(defaults, opts, keep.null = TRUE)
    },
    
    #' Evaluate only: returns log pmf etc. (no sampling)
    pmf = function(time, mark = NULL) {
      setup <- private$setup_common(time, mark)
      private$eval_density(time, setup)
    },
    
    #' Sample only: returns sampled mark & its log density (no provided mark density)
    sample = function(time, mark = NULL) {
      setup <- private$setup_common(time, mark)
      private$generate_mark(time, setup)
    },
    
    #' Evaluate and (optionally) sample in one go (back-compat API)
    compute = function(time, mark = NULL) {
      setup <- private$setup_common(time, mark)
      out_d <- private$eval_density(time, setup)
      out_s <- private$generate_mark(time, setup)
      c(out_d, out_s)
    }
  ),
  
  private = list(
    #' -------- Shared setup (non-bipartite & bipartite) --------
    setup_common = function(time, mark) {
      if (self$bipartite) return(private$setup_bipartite(time, mark))
      # ---- non-bipartite setup (adapted from your mark_setup) ----
      
      # if no mark is supplied take the lastest update 
      if (is.null(mark)) {
        mark <- filtration_to_net(self$filtration, time, equals = TRUE)
        last_net <- filtration_to_net(self$filtration, time, equals = FALSE)
      }else{
        last_net <- filtration_to_net(self$filtration, time, equals = TRUE)
      }
      new_net  <- last_net
      if (!is.null(last_net) && (last_net %n% "n") != 0) {
        new_nodes <- mark %n% "n"
        old_nodes <- last_net %n% "n"
        network::add.vertices(new_net, nv = new_nodes)
        set.vertex.attribute(new_net,
                             "time",
                             c(get.vertex.attribute(last_net, "time"),rep(time, new_nodes)))
      } else {
        last_net <- NULL
        new_net  <- network::network(matrix(0, 1, 1), directed = FALSE)
        set.vertex.attribute(new_net, "time", time)
        old_nodes <- 0
        new_nodes <- 1
      }
      cand <- private$candidate_edges(last_net, old_nodes, new_nodes)
      list(
        mark = mark,
        new_net = new_net,
        last_net = last_net,
        poss_tails = cand$poss_tails,
        poss_heads = cand$poss_heads,
        poss_edges = cand$poss_edges,
        tails = cand$tails,
        heads = cand$heads,
        new_nodes = new_nodes,
        old_nodes = old_nodes
      )
    },
    
    #============================
    # NOT IMPLEMENTED FOR NOW !!!
    # ===========================
    setup_bipartite = function(time, mark) {
      stop("bipartite not supported in the base class yet.")
    },
    
    candidate_edges = function(last_net, old_nodes, new_nodes) {
      poss_tails <- seq.int((old_nodes + 1), new_nodes)
      poss_tails <- poss_tails[poss_tails > 0]
      poss_heads <- seq_len(old_nodes)
      poss_heads <- poss_heads[poss_heads > 0]
      poss_edges <- expand.grid(poss_tails, poss_heads)
      poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2], , drop = FALSE]
      tails <- poss_edges[,1]
      heads <- poss_edges[,2]
      if (!is.null(last_net)) {
        in_old <- vapply(seq_along(heads), function(i)
          length(get.edgeIDs(last_net, heads[i], tails[i])) != 0, logical(1))
        tails <- tails[!in_old]; heads <- heads[!in_old]
      }
      # do the truncation:
      if (!is.null(self$opts$truncation)) {
        if (length(tails) > 0) {
          valid <- (tails > max(0, new_nodes - self$opts$truncation)) |
            (heads > max(0, old_nodes - self$opts$truncation))
          tails <- tails[valid]
          heads <- heads[valid]
          poss_edges <- poss_edges[valid, , drop = FALSE]
        }
      }
      list(poss_tails=poss_tails,
           poss_heads=poss_heads,
           poss_edges=poss_edges,
           tails=tails,
           heads=heads)
    },
    
    #' -------- Density evaluation pipeline --------
    eval_density = function(time, setup) {
      if (!isTRUE(self$opts$generate_density))
        return(list(mark_density=1, log_mark_density=0, edge_probs=NULL))
      # subclasses compute probs aligned to heads/tails:
      edge_probs <- private$edge_probs_from_setup(setup, time)
      
      mark <- setup$mark
      if (is.null(mark) || length(edge_probs) == 0)
        return(list(mark_density=1, log_mark_density=0, edge_probs=edge_probs))
      
      in_mark <- if (is.null(self$opts$new_edge_hash)) {
        vapply(seq_along(setup$heads), function(i)
          length(get.edgeIDs(mark, setup$heads[i], setup$tails[i])) != 0, logical(1))
      } else {
        has_edge(setup$heads, setup$tails, self$opts$new_edge_hash)
      }
      
      # Optional node count term (CS uses this)
      node_dens <- private$node_growth_density(mark, setup$last_net)
      
      if (length(edge_probs) == 1) {
        log_mark_density <- node_dens
      } else {
        log_mark_density <- sum(log(edge_probs[in_mark])) + sum(log(1 - edge_probs[!in_mark])) + log(node_dens)
      }
      list(mark_density = exp(log_mark_density),
           log_mark_density = log_mark_density,
           edge_probs = edge_probs)
    },
    
    generate_mark = function(time, setup) {
      # Default: add a single node (or Poisson many if node_lambda present),
      # recompute candidates/probs, Bernoulli add edges.
      last_net <- setup$mark
      if (!is.null(last_net) && (last_net %n% "n") >= 1) {
        mark_sample <- last_net
        old_nodes <- last_net %n% "n"
        # growth policy
        new_nodes <- private$node_growth_sample(time, last_net)
        if (new_nodes > 0) {
          network::add.vertices(mark_sample, new_nodes)
          set.vertex.attribute(mark_sample, "time",
                               c((last_net %v% "time"), rep(time, new_nodes)))
        }
        new_size <- mark_sample %n% "n"
        poss <- private$candidate_edges(last_net, old_nodes, new_size)
        tails <- poss$tails
        heads <- poss$heads
        
        if (length(tails) == 0) {
          return(list(mark_sample = mark_sample,
                      mark_sample_density = 1,
                      log_mark_sample_density = 0))
        }
        delete.vertex.attribute(mark_sample, "na")
        # Compute probabilities under the grown graph
        #grown_setup <- modifyList(setup, list(new_net = mark_sample, tails = tails, heads = heads))
        grown_setup <- private$setup_common(time, mark = mark_sample)
        # need to update possible edges too?
        #grown_setup$poss_edges <- cbind(tails, heads)
        
        probs <- private$edge_probs_from_setup(grown_setup, time, generation = TRUE)
        add <- stats::runif(length(probs)) < probs
        
        if(any(heads[add] > mark_sample %n% 'n') | any(tails[add] > mark_sample %n% 'n')){
          browser()
        }
        network::add.edges(mark_sample, heads[add], tails[add])
        if (!is.null(mark_sample %e% "time")) {
          set.edge.attribute(mark_sample, "time",
                             c(mark_sample %e% "time", rep(time, sum(add))))
        }
        logp <- sum(log(probs[add])) + sum(log(1 - probs[!add]))
        if (!is.null(self$params$node_lambda))
          logp <- logp + log(stats::dpois(new_nodes - old_nodes, self$params$node_lambda))
        list(mark_sample = mark_sample,
             mark_sample_density = exp(logp),
             log_mark_sample_density = logp)
      } else {
        # initialize singleton
        mark_sample <- if (is.null(last_net)) network::network(matrix(1), directed = FALSE) else last_net
        set.vertex.attribute(mark_sample, "time",
                             c(mark_sample %v% "time", time))
        list(mark_sample = mark_sample,
             mark_sample_density = 1,
             log_mark_sample_density = 0)
      }
    },
    
    #' -------- Hook for subclasses: compute probs vector aligned to heads/tails --------
    edge_probs_from_setup = function(setup, time, generation = FALSE) {
      stop("edge_probs_from_setup() must be implemented in subclass.")
    },
    
    node_growth_sample = function(time, last_net){
      stop("node_growth_sample() must be implemented in subclass.")
    },
    node_growth_density = function(mark, last_net){
      stop("node_growth_density() must be implemented in subclass.")
    }
  )
)

#' @details BA Kernel
BAKernel <- R6::R6Class(
  "BAKernel",
  inherit = MarkKernel,
  private = list(
    edge_probs_from_setup = function(setup, time, generation = FALSE) {
      # degree-based with exponential aging on heads
      last_net <- setup$last_net %||% setup$new_net
      if (is.null(last_net) || (last_net %n% "n") <= 2) return(rep(1, length(setup$heads)))
      times <- get.vertex.attribute(last_net, "time")
      degs  <- sna::degree(last_net)
      factor <- self$decay_fun(times, time, self$params, idx = seq_along(times))
      degs_w <- degs * factor
      tot <- sum(degs_w, na.rm = TRUE)
      if (length(setup$heads) == 0) return(numeric(0))
      if (tot==0) rep(1, length(setup$heads)) else degs_w[setup$heads] / tot
    },
    
    node_growth_sample = function(time, last_net){
      if (is.null(last_net) || (last_net %n% "n") < 4) {
        return(4)
      } else {
        if (time > self$opts$max_node_time) {
          return(0)
        } else if (!is.null(self$params$node_lambda)) {
          return(rpois(1, self$params$node_lambda))
        } else {
          return(1)
        }
      }
    },
    
    node_growth_density = function(mark, last_net){
      new_nodes <- mark %n% "n" - (last_net %n% "n" %||% 0)
      if (!is.null(self$params$node_lambda)) {
        return(stats::dpois(new_nodes, self$params$node_lambda))
      } else {
        return(1)
      }
    }
  )
)

#' @details CS Kernel



#' @details BA_bipartite Kernel



