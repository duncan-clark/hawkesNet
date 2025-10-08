
#' @details Base Class for markKernel - that adds network updates to the network
#' note that the mark space is dynamic so the kernel depends on the filtration
#' so as the filtration grows, the markKernel changes
#' Currently works by supplyig bernoulli edge probabilities and then storing these
#' Designed to be extensible to other edge probability models
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
    initialize = function(params,
                          filtration = NULL,
                          opts = list(),
                          decay_fun = NULL,
                          model = NULL,
                          bipartite = FALSE) {
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
    },
    
    #' Legacy wrapper to match old PMF_mark_* function signatures
    as_legacy_fun = function() {
      parent_kernel <- self  # capture
      function(time,
               params,
               mark_filtration,
               mark = NULL,
               generate_mark = FALSE,
               generate_density = TRUE,
               grad = FALSE,
               new_edge_hash = NULL,
               truncation = NULL,
               ...) {
        
        # clone the kernel so per-call opts/overrides don't mutate parent
        ker <- parent_kernel$clone(deep = FALSE)
        
        # allow call-time overrides, else fall back to what's on the kernel
        if (!missing(params) && !is.null(params)) ker$params <- params
        if (!missing(mark_filtration) && !is.null(mark_filtration)) ker$filtration <- mark_filtration
        
        # update opts for this call (keep other opts intact)
        ker$opts <- utils::modifyList(
          ker$opts,
          list(
            generate_mark    = generate_mark,
            generate_density = generate_density,
            truncation       = truncation,
            new_edge_hash    = new_edge_hash
            # note: grad, formula_RHS, mark_decay, max_node_time stay as set on ker$opts
          ),
          keep.null = TRUE
        )
        
        # compute
        out <- ker$compute(time = time, mark = mark)
        
        # return legacy field names/shape
        list(
          mark_density            = out$mark_density,
          log_mark_density        = out$log_mark_density,
          edge_probs              = out$edge_probs,
          mark_sample             = out$mark_sample,
          mark_sample_density     = out$mark_sample_density,
          log_mark_sample_density = out$log_mark_sample_density
        )
      }
    }
  ),
  
  private = list(
    #' -------- Shared setup (non-bipartite & bipartite) --------
    setup_common = function(time, mark) {
      if (self$bipartite) return(private$setup_bipartite(time, mark))
      # ---- non-bipartite setup (adapted from your mark_setup) ----
    
      if (is.null(mark)) {
        mark <- filtration_to_net(self$filtration, time, equals = TRUE)
      }else{
      }
      # note that if mark must be happening at time t
      last_net <- filtration_to_net(self$filtration, time, equals = FALSE)
      # add warning if mark supplied AND the current ent already has a mark at that time:
      if(max(mark %v% "time") > time){
        warning("You supplied a mark that has nodes with time >= current time; this is likely an error.")
      }
      new_net <- network::network.copy(last_net)
      
      if (!is.null(last_net) && (last_net %n% "n") != 0) {
        new_nodes <- mark %n% "n"
        old_nodes <- last_net %n% "n"
        to_add   <- max(0L, new_nodes - old_nodes)
        network::add.vertices(new_net, nv = to_add)
        set.vertex.attribute(new_net,
                             "time",
                             c(get.vertex.attribute(last_net, "time"),rep(time, to_add)))
      } else {
        last_net <- NULL
        new_net  <- network::network(matrix(0, 1, 1), directed = FALSE)
        set.vertex.attribute(new_net, "time", time)
        old_nodes <- 0
        new_nodes <- 1
      }
      cand <- private$candidate_edges(last_net, new_nodes)
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
    
    candidate_edges = function(last_net, new_nodes) {
      old_nodes <- last_net %n% "n"
      
      
      if(old_nodes == new_nodes){
        poss_tails <- 1:old_nodes
      }else{
        poss_tails <- seq.int(old_nodes + 1, new_nodes)
      }
      
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
      if(max(setup$tails) > (setup$new_net %n% 'n') | max(setup$heads) > (setup$new_net %n% 'n')){
        stop("Internal error: Edge indices exceed number of nodes in mark_sample.")
      }
      
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
      node_dens <- private$node_growth_density(mark,time, setup$last_net)
      
      if (length(edge_probs) == 1) {
        log_mark_density <- log(node_dens)
      } else {
        log_mark_density <- sum(log(edge_probs[in_mark])) + sum(log(1 - edge_probs[!in_mark])) + log(node_dens)
      }
      list(mark_density = exp(log_mark_density),
           log_mark_density = log_mark_density,
           edge_probs = edge_probs)
    },
    
    generate_mark = function(time, setup) {
      last_net <- setup$mark
      if (!is.null(last_net) && (last_net %n% "n") >= 1) {
        mark_sample <- network::network.copy(last_net)
        old_nodes <- last_net %n% "n"
        # growth policy
        new_nodes <- private$node_growth_sample(time, last_net)
        if (new_nodes > 0) {
          network::add.vertices(mark_sample, new_nodes)
          set.vertex.attribute(mark_sample, "time",
                               c((last_net %v% "time"), rep(time, new_nodes)))
        }
        new_size <- mark_sample %n% "n"
        poss <- private$candidate_edges(last_net, new_size)
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
        
        if(max(grown_setup$tails) > (mark_sample %n% 'n') | max(grown_setup$heads) > (mark_sample %n% 'n')){
          stop("Internal error: Edge indices exceed number of nodes in mark_sample.")
        }
        
        probs <- private$edge_probs_from_setup(grown_setup, time, generation = TRUE)
        add <- stats::runif(length(probs)) < probs
        
        if(any(heads[add] > mark_sample %n% 'n') | any(tails[add] > mark_sample %n% 'n')){
          stop("Internal error: Edge indices exceed number of nodes in mark_sample.")
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
        mark_sample <- if (is.null(last_net)) network::network(matrix(1), directed = FALSE) else network::network.copy(last_net)
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
    node_growth_density = function(mark,time, last_net){
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
      if (is.null(last_net) || (last_net %n% "n") <= 2) return(c(1))
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
          return(stats::rpois(1, self$params$node_lambda))
        } else {
          return(1)
        }
      }
    },
    
    node_growth_density = function(mark,time, last_net){
      new_nodes <- mark %n% "n" - (last_net %n% "n" %||% 0)
      if (!is.null(self$params$node_lambda)) {
        return(stats::dpois(new_nodes, self$params$node_lambda))
      } else {
        return(1)
      }
    }
  )
)

#' @details CS Kernel (Change-Statistic HawkesNet)
CSKernel <- R6::R6Class(
  "CSKernel",
  inherit = MarkKernel,
  private = list(
    # ---- Edge probabilities for density or generation ----
    # - For density: diffs = time - node_times[heads]
    # - For generation: diffs = node_times[tails] - node_times[heads]
    edge_probs_from_setup = function(setup, time, generation = FALSE) {
      # If we don't have enough history, fall back to neutral prob 1
      last_net <- setup$last_net %||% setup$new_net
      if (is.null(last_net) || (last_net %n% "n") <= 2) return(c(1))
      
      net <- setup$new_net
      # Ensure a minimum of 4 nodes  before computing change stat
      if ((net %n% "n") < 4) {
        net <- network::add.vertices(net, 4 - (net %n% "n"))
      }
      
      # Safety check: avoid invalid edge indices
      if (length(setup$tails) > 0 && max(setup$tails) > (net %n% "n")) {
        stop("accidentally adding an edge into a network that doesn't have that node yet")
      }
      
      # Clean attributes that could upset the C++ backend
      delete.vertex.attribute(net, "na")
      
      # Get or build ERNM model
      m <- self$model
      if (is.null(m)) {
        frm <- as.formula(paste("net ~ ", self$opts$formula_RHS))
        m <- ernm::createCppModel(frm)
      } else {
        m$setNetwork(ernm::as.BinaryNet(net))
      }
      m$calculate()
      
      cs <- m$computeChangeStats(setup$tails, setup$heads)
      if (length(cs) == 0) {
        stop("these parameters result in full networks; you probably don't want this")
      }
      
      # Logistic probabilities from change stats
      beta <- as.numeric(self$params$CS_params)
      p    <- plogis(as.numeric(cs %*% beta))
      
      # Decay factor
      if (identical(self$opts$mark_decay, "activity")) {
        node_times <- get_latest_times(setup$new_net)
      } else {
        node_times <- setup$new_net %v% "time"
      }
      
      if (!generation) {
        # Legacy density: diffs = time - node_times[heads]
        diffs  <- time - node_times[setup$heads]
      } else {
        # Legacy generation: diffs = node_times[tails] - node_times[heads]
        diffs <- vapply(seq_along(setup$tails), function(i) {
          node_times[setup$tails[i]] - node_times[setup$heads[i]]
        }, numeric(1))
      }
      
      factor <- exp(-self$params$beta_edges * diffs)
      p * factor
    },
    
    # ---- Growth hooks (Poisson growth like legacy; bounded by max_node_time) ----
    node_growth_sample = function(time, last_net) {
      # if size < 4: add padding to 4
      # else if time > max_node_time: 0
      # else: rpois(1, node_lambda)
      n_now <- last_net %n% "n"
      if (n_now < 4) {
        return(4 - n_now)
      }
      if (!is.null(self$opts$max_node_time) && time > self$opts$max_node_time) {
        return(0L)
      }
      as.integer(stats::rpois(1, self$params$node_lambda %||% 0))
    },
    
    node_growth_density = function(mark, time, last_net) {
      # if no previous net, neutral factor
      if (is.null(last_net)) return(1)
      
      new_nodes <- (mark %n% "n") - (last_net %n% "n")
      
      # cutoff: after max_node_time the factor is 1 (i.e., log contribution = 0)
      mt <- self$opts$max_node_time %||% Inf
      if (is.finite(mt) && time > mt) return(1)
      
      # if no lambda provided, neutral
      lambda <- self$params$node_lambda
      if (is.null(lambda)) return(1)
      stats::dpois(new_nodes, lambda)
    },
    
    # ---- Override generate_mark to mirror legacy generation exactly ----
    generate_mark = function(time, setup) {
      last_net <- setup$mark
      if (is.null(last_net) || (last_net %n% "n") < 1) {
        # Initialize
        mark_sample <- if (is.null(last_net)) network::network(matrix(1), directed = FALSE) else last_net
        set.vertex.attribute(mark_sample, "time", c(mark_sample %v% "time", time))
        return(list(mark_sample = mark_sample,
                    mark_sample_density = 1,
                    log_mark_sample_density = 0))
      }
      
      mark_sample <- last_net
      old_nodes   <- last_net %n% "n"
      
      # Legacy growth
      add_nodes <- if (old_nodes < 4) {
        4 - old_nodes
      } else if (!is.null(self$opts$max_node_time) && time > self$opts$max_node_time) {
        0L
      } else {
        as.integer(stats::rpois(1, self$params$node_lambda %||% 0))
      }
      
      if (add_nodes > 0) {
        network::add.vertices(mark_sample, add_nodes)
        set.vertex.attribute(mark_sample, "time",
                             c((last_net %v% "time"), rep(time, add_nodes)))
      }
      new_size <- mark_sample %n% "n"
      
      # Legacy truncation window during generation:
      trunc <- self$opts$truncation %||% 0L
      poss_tails <- seq.int(old_nodes - trunc, new_size)
      poss_heads <- seq.int(new_size - trunc - 1L, new_size)
      poss_tails <- poss_tails[poss_tails > 0]
      poss_heads <- poss_heads[poss_heads > 0]
      
      poss_edges <- expand.grid(poss_tails, poss_heads)
      poss_edges <- poss_edges[poss_edges[, 1] > poss_edges[, 2], , drop = FALSE]
      tails <- poss_edges[, 1]; heads <- poss_edges[, 2]
      
      # Remove edges already in old network
      if (!is.null(last_net) && length(tails) > 0) {
        in_old <- vapply(seq_along(heads), function(i)
          length(get.edgeIDs(last_net, heads[i], tails[i])) != 0, logical(1))
        tails <- tails[!in_old]; heads <- heads[!in_old]
      }
      
      # If no candidates, we're done
      if (length(tails) == 0) {
        return(list(mark_sample = mark_sample,
                    mark_sample_density = 1,
                    log_mark_sample_density = 0))
      }
      
      delete.vertex.attribute(mark_sample, "na")
      
      # Build a "grown" setup to compute generation probs with correct diffs (tails-heads)
      grown_setup <- list(
        mark = last_net,
        new_net = mark_sample,
        last_net = last_net,
        tails = tails,
        heads = heads
      )
      
      if(max(grown_setup$tails) > (mark_sample %n% 'n') | max(grown_setup$heads) > (mark_sample %n% 'n')){
        stop("Internal error: Edge indices exceed number of nodes in mark_sample.")
      }
      
      probs <- private$edge_probs_from_setup(grown_setup, time, generation = TRUE)
      add <- stats::runif(length(probs)) < probs
      
      network::add.edges(mark_sample, heads[add], tails[add])
      set.edge.attribute(mark_sample, "time",
                         c(mark_sample %e% "time", rep(time, sum(add))))
      
      # Sample log-density: product of Bernoullis × Poisson(new-old) (if before cutoff)
      logp <- sum(log(probs[add])) + sum(log(1 - probs[!add]))
      if (is.null(self$opts$max_node_time) || time <= self$opts$max_node_time) {
        if (!is.null(self$params$node_lambda)) {
          logp <- logp + log(stats::dpois((new_size - old_nodes), self$params$node_lambda))
        }
      }
      list(mark_sample = mark_sample,
           mark_sample_density = exp(logp),
           log_mark_sample_density = logp)
    }
  )
)



#' @details BA_bipartite Kernel



