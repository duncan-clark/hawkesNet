

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
#' @param events_list A data frame of interactions
#' between participants,  with elements \code{i} - from node,
#' \code{j} - to node, &
#' \code{t} - time of interaction.
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
events_to_net <- function(events_list,
                          net = NULL,
                          directed = FALSE){

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
#' @param x A data frame ...
#' @export
events_to_bipartite_net <- function(x, pid = "anon_person_id",
                                    eid = "event_id", time = "diff_date"){
    participants <- unique(x[[pid]])
    events <- unique(x[[eid]])
    all_nodes <- c(events, participants)
    net <- network::network.initialize(length(all_nodes),
                                       bipartite = length(events), directed = FALSE,
                                       multiple = FALSE)
    network::network.vertex.names(net) <- all_nodes
    for(i in seq_len(nrow(x))) {
        tail <- match(x[[eid]][i], all_nodes)     
        head <- match(x[[pid]][i], all_nodes)
        if(length(get.edgeIDs(net, tail, head)) == 0) {
            network::add.edge(net, tail = tail, head = head)
        }
    }
    return(net)
}

#' Internal function, takes a
#' bipartite network and makes it into a bipartite igraph
#' mainly useful for plotting
bn_ig <- function(net){
    edges <- network::as.matrix.network.edgelist(net)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% unique(edges[,2]), TRUE, FALSE)
    return(g)
}
#' Plot example subcomponents of a bipartite igraph
#' internal function
plot_example_sub_component <- function(g, size, idx = 1,
                                       cols = c("#E41A1C", "#377EB8"), ...){
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
    ## col; if type == TRUE then cols[1]
    igraph::plot.igraph(sub_g, vertex.color = ifelse(igraph::V(sub_g)$type, cols[1], cols[2]),
         vertex.frame.color = ifelse(igraph::V(sub_g)$type, cols[1], cols[2]), ...)
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
