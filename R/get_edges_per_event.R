#' Calculate the number of new edges added at each event time
#'
#' @param nw A network object with a "time" edge attribute
#' @return A data.frame with columns 'time' and 'n_edges'
#' @examples
#' net <- network::network(matrix(c(1,2, 1,3, 2,3), ncol = 2, byrow = TRUE),
#'                         directed = FALSE)
#' network::set.edge.attribute(net, "time", c(0.1, 0.1, 0.5))
#' get_edges_per_event(net)
#' @export
get_edges_per_event <- function(nw) {
  if (!is.network(nw)) stop("Input must be a network object.")
  
  # Extract edge times
  edge_times <- nw %e% "time"
  if (is.null(edge_times) || length(edge_times) == 0) {
    return(data.frame(time = numeric(0), n_edges = integer(0)))
  }
  
  # Count occurrences of each unique time
  counts <- table(edge_times)
  
  # Convert to data.frame
  df <- data.frame(
    time = as.numeric(names(counts)),
    n_edges = as.integer(counts)
  )
  
  # Sort by time
  df <- df[order(df$time), ]
  rownames(df) <- NULL
  
  return(df)
}

# Example usage (commented out):
# ep_event <- get_edges_per_event(net_day1)
# mean(ep_event$n_edges)
# hist(ep_event$n_edges, breaks = 0:max(ep_event$n_edges) + 0.5)
