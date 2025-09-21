## set base class
setClass("mark_PMF",
  slots = c(
    name = "character",
    params = "list"
  ))
## subclasses
## BA
setClass("BA_PMF",
         contains = "mark_PMF",
         slots = c())
## CS
setClass("CS_PMF",
         contains = "mark_PMF",
         slots = c(formula_RHS,
                   mark_decay = 'node_entrance',
                   model = NULL,
                   max_node_time = NULL))
              
## generic function
setGeneric("PMF", function(object, time,
                           params,
                           mark_filtration,
                           mark = NULL,
                           generate_mark = FALSE,
                           generate_density = TRUE,
                           grad = FALSE,
                           new_edge_hash = NULL,
                           truncation = NULL,
                           ...) standardGeneric("PMF"))

## specific methods
## setup - internal function
setup <- function(time,
                  mark_filtration,
                  mark = NULL,
                  truncation = NULL) {
    ## get the current net and mark
   if (is.null(mark)) {
    mark <- filtration_to_net(mark_filtration, time, equals = TRUE)
  }
  last_net <- filtration_to_net(mark_filtration, time, equals = FALSE)
  new_net <- last_net
  if (last_net %n% "n" != 0) {
    new_nodes <- mark %n% "n"
    old_nodes <- last_net %n% "n"
    network::add.vertices(new_net, nv = new_nodes - old_nodes)
    set.vertex.attribute(new_net, "time",
      c(last_net %v% "time", rep(time, new_nodes))
    )
  } else {
    last_net <- NULL
    new_net <- network::network(matrix(0, 1, 1), directed = FALSE)
    set.vertex.attribute(new_net, "time", time)
    old_nodes <- 0
    new_nodes <- 1
  }
  if (is.null(truncation)) truncation <- 1
  poss_tails <- (old_nodes - truncation + 1):new_nodes
  poss_tails <- poss_tails[poss_tails > 0]
  poss_heads <- 1:old_nodes
  poss_heads <- poss_heads[poss_heads > 0]
  poss_edges <- expand.grid(poss_tails, poss_heads)
  poss_edges <- poss_edges[poss_edges[, 1] > poss_edges[, 2], ]
  tails <- poss_edges[, 1]
  heads <- poss_edges[, 2]
  if (!is.null(last_net)) {
    in_old_net <- sapply(seq_along(heads), function(i) {
      length(get.edgeIDs(last_net, heads[i], tails[i])) != 0
    })
    tails <- tails[!in_old_net]
    heads <- heads[!in_old_net]
  }
  list(mark = mark,
    last_net = last_net,
    new_net = new_net,
    old_nodes = old_nodes,
    new_nodes = new_nodes,
    tails = tails,
    heads = heads)
}
## mark density - internal function
mark_density <- function(probs,
                         mark,
                         last_net,
                         heads,
                         tails,
                         generate_mark = FALSE,
                         grad = FALSE,
                         time = NULL,
                         new_edge_hash = NULL,
                         params = NULL,
                         mark_decay = "node_entrance",
                         model = NULL,
                         max_node_time = NULL) {
  ## inits
  log_mark_density <- 0
  mark_density <- 1
  mark_grad <- 0
  decay_grad <- 0
  mark_sample <- last_net
  if (!is.null(probs) && length(probs) > 0) {
    if (is.null(new_edge_hash)) {
      in_mark <- sapply(seq_along(heads), function(i) {
        length(get.edgeIDs(mark, heads[i], tails[i])) != 0
      })
    } else {
      in_mark <- has_edge(heads, tails, new_edge_hash)
    }

    log_mark_density <- sum(log(probs[in_mark])) + sum(log(1 - probs[!in_mark]))
    mark_density <- exp(log_mark_density)
    if (grad) {
      mark_grad <- rep(0, length(probs))
      decay_grad <- rep(0, length(probs))
    }
  }

  # Optional: generate a new mark sample
  if (generate_mark) {
    mark_sample <- mark
  }
  list(
    mark_density = mark_density,
    log_mark_density = log_mark_density,
    edge_probs = probs,
    mark_grad = mark_grad,
    decay_grad = decay_grad,
    mark_sample = mark_sample,
    mark_sample_density = mark_density, # placeholder
    log_mark_sample_density = log_mark_density
  )
}
## PMFs
setMethod("PMF", "BA_PMF",
          function(object, x) {
              state <- setup(time, mark_filtration, mark, truncation)
              ark <- state$mark
              last_net <- state$last_net
              new_net <- state$new_net
              heads <- state$heads
              tails <- state$tails

                                        # --- BA-specific probability model ---
              if (!is.null(last_net) && (last_net %n% "n" > 2)) {
                  times <- last_net %v% "time"
                  degs <- degree(last_net) * exp(-params$beta_edges * (time - times))
                  total_deg <- sum(degs)
                  probs <- if (total_deg == 0) rep(1, length(heads)) else degs[heads] / total_deg
              } else {
                  probs <- rep(1, length(heads))
              }

              .compute_mark_density(probs,
                                    mark,
                                    last_net,
                                    heads,
                                    tails,
                                    generate_mark,
                                    grad,
                                    time,
                                    new_edge_hash,
                                    params)
          })


setMethod("PMF", "CS_PMF",
          function(object, x) {
              state <- setup(time, mark_filtration, mark, truncation)
              mark <- state$mark
              last_net <- state$last_net
              new_net <- state$new_net
              heads <- state$heads
              tails <- state$tails              
              probs <- NULL
              if (!is.null(last_net) && last_net %n% "n" > 2) {
                  delete.vertex.attribute(new_net, "na")
                  if (is.null(model)) {
                      model <- createCppModel(as.formula(paste("new_net ~ ", formula_RHS)))
                  } else {
                      model$setNetwork(as.BinaryNet(new_net))
                  }
                  model$calculate()
                  change_stats <- model$computeChangeStats(tails, heads)
                  probs <- apply(change_stats, 1, function(c) {
                      1 / (1 + exp(-sum(c * params$CS_params)))
                  })

                  node_times <- if (mark_decay == "activity")
                                    get_latest_times(new_net) else (new_net %v% "time")
                  diffs <- time - node_times[heads]
                  factor <- exp(-params$beta_edges * diffs)
                  probs <- probs * factor
              }

              mark_density(probs,
                           mark,
                           last_net,
                           heads,
                           tails,
                           generate_mark,
                           grad,
                           time,
                           new_edge_hash,
                           params,
                           mark_decay,
                           model,
                           max_node_time)
          })
