library(hawkesNet)
library(dplyr)
library(network)

# 1. Load and jitter (exactly as in the script)
raw <- read.table(system.file("extdata", "ht09_contact_list.dat", package = "hawkesNet"))
df <- data.frame(
  time = raw$V1 / 20 / 3600,
  from = raw$V2,
  to   = raw$V3
)
set.seed(1)
df$time <- df$time + rnorm(nrow(df), 0, 0.01 / 3600)

# 2. Subset first session
subset_first_session <- function(df, gap_threshold = 1.0) {
  df <- df[order(df$time), , drop = FALSE]
  gaps <- diff(df$time)
  first_gap <- which(gaps > gap_threshold)[1]
  if (is.na(first_gap)) return(df)
  df[seq_len(first_gap), , drop = FALSE]
}
df_day1 <- subset_first_session(df, gap_threshold = 1.0)

# 3. Make net (exactly as in the script)
make_hypertext_net <- function(df, use_first_contact_only = TRUE, max_edges = 0L) {
  swap <- df$from > df$to
  df[swap, c("from", "to")] <- df[swap, c("to", "from")]
  df <- df %>% distinct()
  df$time <- as.numeric(df$time)
  df <- df[is.finite(df$time), ]
  df$time <- df$time - min(df$time)
  if (use_first_contact_only) {
    df <- df %>%
      dplyr::group_by(.data$from, .data$to) %>%
      dplyr::summarise(time = min(.data$time), .groups = "drop")
  }
  df <- df %>% arrange(time)
  nodes_raw <- sort(unique(c(df$from, df$to)))
  entry_time_by_node <- vapply(nodes_raw, function(v) {
    min(df$time[df$from == v | df$to == v], na.rm = TRUE)
  }, numeric(1))
  ord_nodes <- order(entry_time_by_node, nodes_raw)
  nodes <- nodes_raw[ord_nodes]
  id_map <- setNames(seq_along(nodes), nodes)
  df$tail <- unname(id_map[as.character(df$from)])
  df$head <- unname(id_map[as.character(df$to)])
  el <- as.matrix(df[, c("tail", "head")])
  net <- network::network(el, matrix.type = "edgelist", directed = FALSE)
  network::set.edge.attribute(net, "time", df$time)
  list(net = net, edges = df)
}

obj_day1 <- make_hypertext_net(df_day1, use_first_contact_only = TRUE)
net_day1 <- obj_day1$net
cat("Events in net_day1:", length(get_times(net_day1)$times), "\n")
cat("Edges in net_day1:", network.edgecount(net_day1), "\n")
cat("Time window:", range(get_times(net_day1)$times), "\n")
