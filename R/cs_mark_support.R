# A mark adds nodes and edges; it cannot replace attributes of existing nodes.
# Check the birth timestamps and modeled attributes independently of edge keys.
.cs_preserves_old_vertex_data <- function(old, proposed, params) {
  n <- network::network.size(old)
  if (!n) return(TRUE)
  for (attr in unique(c("time", names(params$vertex_categorical)))) {
    if (!(attr %in% network::list.vertex.attributes(old))) next
    if (!(attr %in% network::list.vertex.attributes(proposed))) return(FALSE)
    before <- network::get.vertex.attribute(old, attr)
    after <- network::get.vertex.attribute(proposed, attr)
    if (length(after) < n || !isTRUE(all.equal(before, after[seq_len(n)],
                                             tolerance = 0, check.attributes = FALSE))) {
      return(FALSE)
    }
  }
  TRUE
}
