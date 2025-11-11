

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
  existing_keys <- hash::keys(edge_hash)
  keys_to_check <- paste(i, j, sep = "-")
  keys_to_check %in% existing_keys
}


#' Convert event list to network
#'
#' Constructs a network from a data frame of event interactions, where each event
#' specifies a source node, a target node, and a timestamp. Vertices and edges are
#' added as needed, and the time attribute is set for each.
#'
#' @param events Data frame or list. Must contain columns \code{i} (from node), \code{j} (to node),
#' and \code{t} (time of interaction).
#' @param net Network or NULL. Optionally, an existing network object to add events to.
#' If \code{NULL}, a new network is constructed. Default: \code{NULL}.
#' @param directed Logical. Whether the network should be directed. Default: \code{FALSE}.
#'
#' @return A \code{network} object with time-stamped edges and vertices. The network attribute
#' \code{n} is set to the largest node index appearing in the events.
#'
#' @details
#' Each row in \code{events_list} creates (if necessary) the specified nodes and adds an
#' edge, with time attributes set for both edge and involved vertices.
#' @examples
#' \dontrun{
#' if(interactive()){
#' events_list <- list(
#' i = c(1, 2, 3, 1),
#' j = c(2, 3, 1, 4),
#' t = c(5, 10, 15, 20))
#' net <- events_to_net(events_list)
#' plot(net)
#'  }
#' }
#' @seealso
#'  \code{\link[network]{network}}, \code{\link[network]{attribute.methods}}, \code{\link[network]{add.vertices}}
#' @rdname events_to_net
#' @export
events_to_net <- function(events,
                          net = NULL,
                          directed = FALSE){
    ## validate input
    if(!inherits(events, "data.frame") & !inherits(events, "list"))
        stop("'events' should be either a list or a data frame.")
    if(sum(c("i", "j", "t") %in% names(events)) != 3)
        stop("'events' should have elements 'i', 'j', & 't'. Please see details.")
    if (!is.numeric(events$i) || any(events$i %% 1 != 0))
        stop("'i' must be whole numbers (i.e., node indices)")
    if (!is.numeric(events$j) || any(events$j %% 1 != 0)) 
        stop("'j' must be whole numbers (i.e., node indices)")
    if (!is.numeric(events$t)) stop("'t' must be numeric (times)")
    if (length(events$i) != length(events$j) ||
        length(events$i) != length(events$t)) {
        stop("'i', 'j', and 't' must all have the same length")
    }
    ## rename to match initial
    events_list <- events
    ## df --> list
    if (inherits(events_list, "data.frame")) events_list <- as.list(events_list)
    stopifnot(is.list(events_list))
    ## node indecies do not have to be numbered consectutively
    for(k in seq_along(events_list$i)){
        i <- events_list$i[k]
        j <- events_list$j[k]
        t <- events_list$t[k]
        if(k == 1 & is.null(net)){
            net <- network::network(matrix(c(i, j), nrow = 1), directed = directed)
      if (!(i %in% get.vertex.attribute(net, "vertex.names"))) {
        network::set.vertex.attribute(net,"vertex.names", i, v = 1)
      }
      if (!(j %in% get.vertex.attribute(net, "vertex.names"))) {
        network::set.vertex.attribute(net,"vertex.names", j, v = 2)
      }
      network::set.vertex.attribute(net,"time", t, v = 1)
      network::set.vertex.attribute(net,"time", t, v = 2)
    } else {
      N <- net %n% 'n'
      over_i <- i - N
      over_j <- j - N
      if(over_i > 0){
        net <- network::add.vertices(net, over_i)
      }
      if(over_j > 0){
        net <- network::add.vertices(net, over_j)
      }
      existing_names <- get.vertex.attribute(net, "vertex.names")
      if (!(i %in% existing_names)) {
        v_index <- which(is.na(existing_names))[1]  
        network::set.vertex.attribute(net, "vertex.names", i, v = v_index)
        network::set.vertex.attribute(net, "time", t, v = v_index)
      } else {
        v_index <- which(existing_names == i)
        network::set.vertex.attribute(net, "time", t, v = v_index)
      }
      if (!(j %in% existing_names)) {
        v_index <- which(is.na(existing_names))[1]
        network::set.vertex.attribute(net, "vertex.names", j, v = v_index)
        network::set.vertex.attribute(net, "time", t, v = v_index)
      } else {
        v_index <- which(existing_names == j)
        network::set.vertex.attribute(net, "time", t, v = v_index)
      }
      add.edge(net, i, j)
    }
    e <- get.dyads.eids(net, i, j)
    if(!is.na(e[[1]])){
      network::set.edge.attribute(net, "time", t, e = e[[1]])
    }
  }
  network::set.network.attribute(net,'n',max(c(events_list$i,events_list$j)))
  return(net)
}


#' Convert event-participant data to bipartite network
#'
#' Creates a bipartite network from a data frame of participant-event
#' memberships, with time attributes for each edge and node.
#'
#' @return A bipartite \code{network} object, where vertices correspond to participants
#' and events, and edges represent participation, with time-stamped attributes.
#'
#' @details
#' The resulting network is undirected and bipartite (events and participants).
#' Edges are added only if not already present. Time attributes are set for each edge
#' @inheritParams events_to_net
#' @examples
#' \dontrun{
#' if(interactive()){
#' events_list <- list(
#' i = c("A", "B", "C", "D","D"),
#' j = c(2, 3, 1, 2, 1),
#' t = c(5, 10, 15, 20, 21))
#' bip <- events_to_bipartite_net(events_list)
#' plot(bip, vertex.col = ifelse(get.vertex.attribute(bip, "type") == "type_i", "#E41A1C", "#377EB8"))
#'  }
#' }
#' @rdname events_to_bipartite_net
#' @export
events_to_bipartite_net <- function(events){
    events_list <- events
    ## make sure node labels are numeric
    ## for bipartite setting
    type_i <- unique(events_list$i)
    type_j <- unique(events_list$j)
    type_i_map <- stats::setNames(seq_along(type_i), type_i)
    type_j_map <- stats::setNames(seq_along(type_j) + length(type_i), type_j)
    events_list <- list(
        i = as.numeric(type_i_map[events_list$i]),
        j = as.numeric(type_j_map[events_list$j]),
        t = events_list$t
    )
    net <- events_to_net(events_list)
    network::set.network.attribute(net, "bipartite", length(unique(events_list$i)))
    network::set.vertex.attribute(net, "name", c(unique(events_list$i), unique(events_list$j)))
    network::set.vertex.attribute(net, "type", c(rep("type_i", length(unique(events_list$i))),
                                                 rep("type_j", length(unique(events_list$j)))))
    return(net)
}

#' Internal function, takes a bipartite network as returned by \link{events_to_bipartite_net}
#' and makes it into a bipartite igraph; mainly useful for plotting.
#' @noRd
bn_ig <- function(net){
    edges <- network::as.matrix.network.edgelist(net)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% unique(edges[,2]), "type_j", "type_i")
    igraph::V(g)$time <- network::get.vertex.attribute(net, "time")
    igraph::E(g)$time <- network::get.edge.attribute(net, "time")
    return(g)
}
#' Plot example subcomponents of a bipartite igraph as returned by \link{bn_ig}
#' (internal function). if \code{grow} is `TRUE` then an animation is returned.
#' @noRd
plot_example_sub_component <- function(g, size, idx = 1,
                                       cols = c("#E41A1C", "#377EB8"), grow = FALSE,
                                       gif.name = "network_growth.gif", interval = 0.8, ...){
    x <- igraph::components(g)
    if(!size %in% x$csize){
        stop(paste("`size` must be one of", paste(names(table(x$csize)), collapse = ", ")))
    }
    if(idx > length(which(x$csize == size))){
        stop(paste("There are only", length(which(x$csize == size)), "components of size", size))
    }
    comp_id <- which(x$csize == size)[idx]
    nodes <- igraph::V(g)$name[x$membership == comp_id]
    sub_g <- igraph::induced_subgraph(g, vids = nodes)
    ## col; if type == "type_i" then cols[1]
    if(grow == FALSE){
        igraph::plot.igraph(sub_g, vertex.color = ifelse(igraph::V(sub_g)$type == "type_i", cols[1], cols[2]),
                            vertex.frame.color = ifelse(igraph::V(sub_g)$type == "type_i", cols[1], cols[2]), ...)
    }else{
        tidx <- igraph::V(sub_g)$time |> sort() |> unique()
        layout_static <- igraph::layout_with_fr(sub_g)
        animation::saveGIF({
            for (t in tidx) {
                eidx <- which(igraph::E(sub_g)$time <= t)
                g_t <- igraph::subgraph_from_edges(sub_g, eids = eidx, delete.vertices = TRUE)
                idx <- match(igraph::V(g_t)$name, igraph::V(sub_g)$name)
                igraph::plot.igraph(g_t,layout = matrix(layout_static[idx, ], ncol = 2),
                                    vertex.color = ifelse(igraph::V(g_t)$type == "type_i", cols[1], cols[2]),
                                    vertex.frame.color = ifelse(igraph::V(g_t)$type == "type_i", cols[1], cols[2]),
                                    xlim = range(layout_static[,1]),
                                    ylim = range(layout_static[,2]),
                                    rescale = FALSE, ...)
            }
        }, movie.name = gif.name, interval = interval)
    }
}

# Active edge IDs are the non-NULL slots of net$mel
.active_eids <- function(net) which(!vapply(net$mel, is.null, FALSE))

#' @title edge_ids_where
#' @description helper function to make sure we delete the right edges when a network has been
#' subject to deletion of vertices and edges
#' @param net network object
#' @param attr edge attribute to apply condition to 
#' @param condition function on values returns vector of Bools
#' @return return EIDs of the edges that meet tehe condition - for safe feeding into delete.edges


# Return TRUE internal edge IDs where attr satisfies `condition`
edge_ids_where <- function(net, attr, condition) {
  vals <- network::get.edge.attribute(net, attr)
  if (!length(vals)) return(integer(0))
  idx  <- which(condition(vals))           # indices in the *attribute vector*
  .active_eids(net)[idx]                   # map to real edge IDs in net$mel
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
filtration_to_net <- function(net, t, equals = FALSE){
  if (is.null(net)){
    net <- network::network(matrix(0, 0, 0), directed = FALSE, bipartite = 0)
    return(net)
  }
  bip <- network::get.network.attribute(net, "bipartite")
  e_delete <- edge_ids_where(net, "time", function(x) x > t)
  if (length(e_delete) > 0) delete.edges(net, e_delete)
  v_times <- tryCatch(network::get.vertex.attribute(net, "time"), error = function(e) NULL)
  if (!is.null(v_times)) delete.vertices(net, which(v_times > t))
  if (!equals){
    e_times <- tryCatch(network::get.edge.attribute(net,"time"), error = function(e) numeric(0))
    n_times <- tryCatch(network::get.vertex.attribute(net,"time"), error = function(e) numeric(0))
    times <- c(e_times, n_times)
    if (length(times) == 0 || all(is.na(times))) {
      t_to_delete <- t
    } else {
      t_to_delete <- max(times, na.rm = TRUE)
    }
    e_delete <- edge_ids_where(net, "time", function(x) x == t_to_delete)
    if (length(e_delete) > 0) network::delete.edges(net, e_delete)
    n_times <- tryCatch(network::get.vertex.attribute(net,"time"), error = function(e) numeric(0))
    if (length(n_times) > 0) network::delete.vertices(net, which(n_times == t_to_delete))
  }
  if ("vertex.names" %in% network::list.vertex.attributes(net)) network::delete.vertex.attribute(net,'vertex.names')
  if (!is.null(bip)) network::set.network.attribute(net, "bipartite", bip)
  return(net)
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
get_times <- function(net){
  node_times <- get.vertex.attribute(net,"time")
  edge_times <- get.edge.attribute(net,"time")
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

# helper function to get latest edge each node is involved with:
get_latest_times <- function(nw){
  el    <- as.matrix.network.edgelist(nw, names = FALSE)
  times <- get.edge.attribute(nw, "time")
  
  latest_edge <- sapply(seq_len(network.size(nw)), function(v) {
    inc <- which(el[,1] == v | el[,2] == v)
    if (length(inc) == 0) return(NA)       # no edges for this node
    inc[which.max(times[inc])]              # index of max-time edge
  })
  
  latest_times <- el    <- as.matrix.network.edgelist(nw, names = FALSE)
  times <- get.edge.attribute(nw, "time")
  
  latest_edge <- sapply(seq_len(network.size(nw)), function(v) {
    inc <- which(el[,1] == v | el[,2] == v)
    if (length(inc) == 0) return(NA)
    inc[which.max(times[inc])]
  })
  
  latest_times <- times[latest_edge]
  if(is.null(latest_times)) {
    latest_times <- rep(NA, network.size(nw))
  }
  node_times <- get_times(nw)
  latest_times[which(is.na(latest_times))] <- node_times$node_times[which(is.na(latest_times))]
  
  return(latest_times)
}


#' Null-coalescing infix operator
#'
#' Returns the left-hand side if it is not \code{NULL}, otherwise returns the right-hand side.
#' Useful for providing default values in expressions.
#'
#' @name %||%
#' @usage a %||% b
#' @param a An object, possibly \code{NULL}.
#' @param b A fallback value to return if \code{a} is \code{NULL}.
#' @return \code{a} if not \code{NULL}, otherwise \code{b}.
#' @examples
#' x <- NULL
#' y <- 10
#' x %||% y   # returns 10
#'
#' z <- 5
#' z %||% y   # returns 5
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
