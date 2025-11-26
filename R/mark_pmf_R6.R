#' MarkKernel R6 Class
#'
#' @description
#' Base class for mark-update kernels on dynamic mark spaces. This class updates a network structure by adding nodes and edges according to supplied probabilities and filtration. Subclasses provide specific edge probability logic.
#'
#' @docType class
#' @format An \code{R6Class} object.
#'
#' @field params List of numeric parameters controlling model behavior.
#' @field filtration Filtration object, representing network history.
#' @field opts List of options (e.g., generate_mark, generate_density, truncation, max_node_time, etc.).
#' @field bipartite Logical; TRUE if bipartite network.
#' @field model Optional compiled model or adapter, e.g. for ERNM.
#' @field decay_fun Function controlling temporal decay of probabilities.
#'
#' @section Methods:
#' \describe{
#'   \item{initialize}{Create a new MarkKernel object.}
#'   \item{pmf}{Evaluate only: returns log pmf etc. (no sampling).}
#'   \item{sample}{Sample only: returns sampled mark & its log density (no provided mark density).}
#'   \item{compute}{Evaluate and (optionally) sample in one go (back-compat API).}
#'   \item{as_legacy_fun}{Legacy wrapper to match old PMF_mark_* function signatures.}
#' }
#'
#' @seealso
#' \link{BAKernel}, \link{CSKernel}
MarkKernel <- R6::R6Class(
                      "MarkKernel",
                      public = list(
                          #' @description Initialize a MarkKernel object.
                          #' @param params List of numeric parameters.
                          #' @param filtration Filtration object.
                          #' @param opts List of options (e.g., generate_mark, generate_density, truncation, max_node_time, etc.).
                          #' @param decay_fun Function for temporal decay of probabilities.
                          #' @param model Optional compiled model or adapter.
                          #' @param bipartite Logical indicating if network is bipartite.
                          #' @return A new MarkKernel object.
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
                          #' @description Evaluate only: returns log pmf etc. (no sampling).
                          #' @param time Numeric time for evaluation.
                          #' @param mark Optional network mark object.
                          #' @return List containing log pmf and related quantities.
                          pmf = function(time, mark = NULL) {
                              setup <- private$setup_common(time, mark)
                              private$eval_density(time, setup)
                          },
                          #' @description Sample only: returns sampled mark & its log density (no provided mark density).
                          #' @param time Numeric time for sampling.
                          #' @param mark Optional network mark object.
                          #' @return List containing the sampled mark and its log density.
                          sample = function(time, mark = NULL) {
                              setup <- private$setup_common(time, mark)
                              private$generate_mark(time, setup)
                          },
                          #' @description Evaluate and (optionally) sample in one go (back-compat API).
                          #' @param time Numeric time for computation.
                          #' @param mark Optional network mark object.
                          #' @return Combined results from evaluation and (optionally) sampling.
                          compute = function(time, mark = NULL) {
                              setup <- private$setup_common(time, mark)
                              out_d <- private$eval_density(time, setup)
                              if(max(setup$new_net %v% 'time') >= time){
                                  out_s <- NULL
                              }else{
                                  out_s <- private$generate_mark(time, setup)
                              }
                              c(out_d, out_s)
                          },
                          #' @description Legacy wrapper to match old PMF_mark_* function signatures.
                          #' @return A function with legacy-compatible arguments and return shape.
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
                                  if (to_add > 0) network::add.vertices(new_net, nv = to_add)
                                  network::set.vertex.attribute(new_net,
                                                       "time",
                                                       c(network::get.vertex.attribute(last_net, "time"),
                                                         rep(time, to_add)))
                              } else {
                                  last_net <- NULL
                                  new_net  <- network::network(matrix(0, 1, 1), directed = FALSE)
                                  network::set.vertex.attribute(new_net, "time", time)
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
                          
                          ##============================
                          ##  Bipartite
                          ## ===========================
                          setup_bipartite = function(time, mark) {
                              if(is.null(mark)){
                                  mark <- filtration_to_net(self$filtration, time, equals = TRUE)
                              }
                              last_net <- filtration_to_net(self$filtration, time, equals = FALSE)
                              new_net <- network::network.copy(last_net)
                              if(is.null(last_net) || (last_net %n% "n") == 0){
                                  new_net <- network::network(matrix(0, 1, 1), directed = FALSE, bipartite = 0)
                                  network::set.vertex.attribute(new_net, "time", time)
                                  network::set.vertex.attribute(new_net, "type", "type_j") # was event 
                                  old_nodes <- 0
                                  new_nodes <- 1
                                  last_net <- network::network(matrix(0,0,0), directed = FALSE, bipartite = 0)
                              } else {
                                  new_nodes <- mark %n% "n"
                                  old_nodes <- last_net %n% "n"
                                  network::add.vertices(new_net, nv = new_nodes - old_nodes)
                                  network::set.vertex.attribute(new_net, "time",
                                                                c(network::get.vertex.attribute(last_net, "time"),
                                                                  rep(time, new_nodes - old_nodes)))
                                  types <- network::get.vertex.attribute(last_net, "type")
                                  if(is.null(types)) types <- rep("type_i", old_nodes)        # was "perp"
                                  new_types <- c(types, rep("type_j", new_nodes - old_nodes)) # event to type_j
                                  network::set.vertex.attribute(new_net, "type", new_types)
                              }
                              
                              cand <- private$candidate_edges(last_net, new_nodes, bipartite = TRUE)
                              list(mark = mark,
                                   new_net = new_net,
                                   last_net = last_net,
                                   poss_tails = cand$poss_tails,
                                   poss_heads = cand$poss_heads,
                                   poss_edges = cand$poss_edges,
                                   tails = cand$tails,
                                   heads = cand$heads,
                                   new_nodes = new_nodes,
                                   old_nodes = old_nodes)
                              
                          },
                          
                          candidate_edges = function(last_net, new_nodes, bipartite = FALSE) {
                              if(is.null(last_net)){
                                  tails <- integer(0)
                                  heads <- integer(0)
                                  return(list(poss_tails = tails,
                                              poss_heads = heads,
                                              poss_edges = matrix(numeric(0), ncol=2),
                                              tails = tails,
                                              heads = heads))
                              }
                              old_nodes <- last_net %n% "n"
                              
                              if(old_nodes == new_nodes){
                                  poss_tails <- 1:old_nodes
                              }else{
                                  poss_tails <- seq.int(old_nodes + 1, new_nodes)
                              }
                              if(!bipartite){
                                  poss_tails <- poss_tails[poss_tails > 0]
                                  poss_heads <- seq_len(old_nodes)
                                  poss_heads <- poss_heads[poss_heads > 0]
                                  poss_edges <- expand.grid(poss_tails, poss_heads)
                                  poss_edges <- poss_edges[poss_edges[,1] > poss_edges[,2], , drop = FALSE]
                                  tails <- poss_edges[,1]
                                  heads <- poss_edges[,2]
                              } else {
                                  type_i_nodes <- which(network::get.vertex.attribute(last_net, "type") == "type_i")
                                  type_j_nodes <- which(network::get.vertex.attribute(last_net, "type") == "type_j")
                                  poss_tails <- seq.int(from = old_nodes + 1, to = new_nodes)
                                  poss_tails <- poss_tails[poss_tails > 0]
                                  poss_heads <- type_i_nodes
                                  poss_heads <- poss_heads[poss_heads > 0]
                                  poss_edges <- expand.grid(poss_tails, poss_heads)
                                  colnames(poss_edges) <- c("tail", "head")
                                  tails <- poss_edges[, "tail"]
                                  heads <- poss_edges[, "head"]
                              }
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
                              n <- network.size(setup$new_net)
                              h <- setup$heads; t <- setup$tails
                              if (length(h) && length(t) &&
                                  (max(h, na.rm = TRUE) > n || max(t, na.rm = TRUE) > n ||
                                   min(h, na.rm = TRUE) < 1 || min(t, na.rm = TRUE) < 1)){
                                  stop(sprintf("Internal error: Edge indices outside [1,%d] in mark_sample.", n), call. = FALSE)
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
                                        # can only generate mark at a time if the last net is before that time:
                              if(max(last_net %v% 'time' >= time)){
                                  stop("You cannot generate a mark at time ", time, " because the supplied mark already has nodes at or after that time.")
                              }
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
                                  
                                  if((sum(is.na(heads[add])) + sum(is.na(tails[add]))) > 0) browser()
                                  if(any(heads[add] > mark_sample %n% 'n') | any(tails[add] > mark_sample %n% 'n')){
                                      browser()
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
                                  mark_sample <- if (is.null(last_net))
                                                     network::network(matrix(1), directed = FALSE) else network::network.copy(last_net)
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

#' BAKernel R6 Class
#'
#' @description
#' Implements a Barabási–Albert (BA) style kernel for dynamic networks, inheriting from MarkKernel. Edge probabilities are degree-based, with exponential aging on heads (nodes).
#'  Node growth is controlled by Poisson or fixed logic.
#'
#' @docType class
#' @section Methods:
#' \describe{
#'   \item{initialize}{Create a new BAKernel object.}
#'   \item{pmf}{Evaluate only: returns log pmf etc. (no sampling).}
#'   \item{sample}{Sample only: returns sampled mark & its log density (no provided mark density).}
#'   \item{compute}{Evaluate and (optionally) sample in one go (back-compat API).}
#'   \item{as_legacy_fun}{Legacy wrapper to match old PMF_mark_* function signatures.}
#' }
#' @seealso
#' \link{MarkKernel}, \link{CSKernel}
BAKernel <- R6::R6Class(
                    "BAKernel",
                    inherit = MarkKernel,
                    private = list(
                        edge_probs_from_setup = function(setup, time, generation = FALSE) {
                                        # degree-based with exponential aging on heads
                            last_net <- setup$last_net %||% setup$new_net
                            if (is.null(last_net) || (last_net %n% "n") <= 2) return(c(1))
                            times <- network::get.vertex.attribute(last_net, "time")
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
                            if(is.null(last_net)){
                                old_nodes <- 0
                            }else{
                                old_nodes <- last_net %n% "n"
                            }
                            new_nodes <- (mark %n% "n") - old_nodes
                            if (!is.null(self$params$node_lambda)) {
                                return(stats::dpois(new_nodes, self$params$node_lambda))
                            } else {
                                return(1)
                            }
                        }
                    )
                )

#' CSKernel R6 Class
#'
#' @description
#' Implements a Change-Statistic (CS) kernel for HawkesNet, inheriting from MarkKernel. Uses ERNM-style change statistics to compute edge probabilities. Node growth is controlled by a Poisson process or cutoff logic. Designed for models where edge probabilities depend on network statistics.
#'
#' @docType class
#' @format An \code{R6Class} object inheriting from \code{MarkKernel}.
#' @section Methods:
#' \describe{
#'   \item{initialize}{Create a new CSKernel object.}
#'   \item{pmf}{Evaluate only: returns log pmf etc. (no sampling).}
#'   \item{sample}{Sample only: returns sampled mark & its log density (no provided mark density).}
#'   \item{compute}{Evaluate and (optionally) sample in one go (back-compat API).}
#'   \item{as_legacy_fun}{Legacy wrapper to match old PMF_mark_* function signatures.}
#' }
#'
#' @seealso
#' \link{MarkKernel}, \link{BAKernel}
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
                                frm <- as.formula(paste0("net ~ ", self$opts$formula_RHS))
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
                            network::set.edge.attribute(mark_sample, "time",
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

#'  BA_BipartiteKernel R6 Class
#' @description
#' Implements a Barabási–Albert (BA) style kernel for bipartite dynamic networks, inheriting from MarkKernel.
#' @docType class
#' @seealso
#' \link{MarkKernel}, \link{BAKernel}

BABipartiteKernel <- R6::R6Class(
  "BABipartiteKernel",
  inherit = MarkKernel,
  private = list(
    edge_probs_from_setup = function(setup, time, generation = FALSE) {
      last_net <- setup$last_net %||% setup$new_net
      if(is.null(last_net) || (last_net %n% "n") <= 1) return(c(1))
      candidate_heads <- setup$heads
      if(length(candidate_heads) == 0) return(numeric(0))

      el <- network::as.edgelist(last_net, names = FALSE)
      edge_times <- if(nrow(el) > 0) last_net %e% "time" else numeric(0)
      times <- network::get.vertex.attribute(last_net, "time")
      degs  <- sna::degree(last_net)
      type_i_nodes <- which(network::get.vertex.attribute(last_net, "type") == "type_i")
      beta_e <- self$params$beta_edges %||% 0.1  # small default if missing
      beta_e <- max(0, beta_e)

      scores <- numeric(length(candidate_heads))
      for(i in seq_along(candidate_heads)) {
        h <- candidate_heads[i]
        if(!(h %in% type_i_nodes)) {
          scores[i] <- 1e-8
          next
        }

        incident_idx <- which(el[,1] == h | el[,2] == h)
        if(length(incident_idx) == 0) {
          scores[i] <- 1 + exp(-beta_e * (time - times[h]))
        } else {
          et <- edge_times[incident_idx]
          diffs <- pmax(time - et, 0)
          scores[i] <- sum(exp(-beta_e * diffs)) + 1  ## mimics BA offset
        }
      }

      probs <- scores / sum(scores)
      return(probs)
    },
    node_growth_sample = function(time, last_net) {
      K <- stats::rpois(1, lambda = self$params$lambda_new %||% 0)
      list(U = 0, V = K)
    },
    node_growth_density = function(mark, time, last_net) {
      old_nodes <- if(is.null(last_net)) 0 else last_net %n% "n"
      new_nodes <- (mark %n% "n") - old_nodes
      stats::dpois(new_nodes, self$params$lambda_new %||% 1)
    }

  )
)
