# Speed test: waiting_times_between_formations
# Compare current approach vs ERNM-based approach
# Run from package root: Rscript inst/openalex_study/test_waiting_times_speed.R

devtools::load_all(".")
library(ernm)
library(network)
library(sna)

# Generate a test network for benchmarking
# Note: This only generates the network, does NOT run the full fit
set.seed(42)
TIME <- 10
params <- list(
  mu = 15, beta_overall = 2, K = 0.5,
  beta_edges = 1, node_lambda = 1.5,
  CS_params = c(-4.0, 1.0, 0.1, -0.03)  # Less negative to get denser network
)

cat("Generating test network (simulation only, no fitting)...\n")
sim <- sim_hawkesGrowthNet(
  params = params,
  time_window = c(0, TIME),
  PMF_mark = PMF_mark_CS,
  cond_intensity = cond_intensity,
  hashed_edges = TRUE,
  verbose = FALSE,
  mu_multiplier = 3,
  truncation = 100,
  formula_RHS = "edges + triangles() + star(c(2,3))"
)
test_net <- sim$net
n_events <- length(get_times(test_net)$times)
n_nodes <- network::network.size(test_net)
n_edges <- network::network.edgecount(test_net)
cat("Test network:", n_nodes, "nodes,", n_events, "events,", n_edges, "edges\n")
cat("Average degree:", round(2 * n_edges / n_nodes, 2), "\n\n")

# Current approach (manual incremental)
waiting_times_current <- function(net, time_attr = "time") {
  out <- list(triangle = numeric(), star2 = numeric(), star3 = numeric())
  el <- network::as.edgelist(net)
  if (nrow(el) == 0) return(out)
  
  edge_times <- NULL
  if (time_attr %in% network::list.edge.attributes(net)) {
    edge_times <- network::get.edge.attribute(net, time_attr)
  }
  if (is.null(edge_times) || length(edge_times) != nrow(el)) {
    vtimes <- network::get.vertex.attribute(net, time_attr)
    if (is.null(vtimes) || all(is.na(vtimes))) return(out)
    edge_times <- pmax(vtimes[el[, 1]], vtimes[el[, 2]], na.rm = TRUE)
  }
  
  ord <- order(edge_times)
  el <- el[ord, , drop = FALSE]
  edge_times <- edge_times[ord]
  event_times <- unique(edge_times)
  n <- network::network.size(net)
  g <- network::network.initialize(n, directed = network::is.directed(net))
  
  t_tri <- numeric()
  t_2s <- numeric()
  t_3s <- numeric()
  n_tri_prev <- 0
  n_2s_prev <- 0
  n_3s_prev <- 0
  
  for (t_cur in event_times) {
    idx <- which(edge_times == t_cur)
    
    # Incremental triangle counting
    new_edges <- el[idx, , drop = FALSE]
    new_triangles <- 0L
    for (j in seq_len(nrow(new_edges))) {
      u <- new_edges[j, 1]
      v <- new_edges[j, 2]
      nbrs_u <- network::get.neighborhood(g, u, type = "combined")
      nbrs_v <- network::get.neighborhood(g, v, type = "combined")
      common <- intersect(nbrs_u, nbrs_v)
      new_triangles <- new_triangles + length(common)
    }
    
    for (j in idx) {
      network::add.edges(g, tail = el[j, 1], head = el[j, 2])
    }
    
    n_tri <- n_tri_prev + new_triangles
    if (n_tri > n_tri_prev) t_tri <- c(t_tri, t_cur)
    n_tri_prev <- n_tri
    
    degs <- sna::degree(g, gmode = "graph")
    n_2s <- sum(choose(degs, 2))
    n_3s <- sum(choose(degs, 3))
    if (n_2s > n_2s_prev) t_2s <- c(t_2s, t_cur)
    if (n_3s > n_3s_prev) t_3s <- c(t_3s, t_cur)
    n_2s_prev <- n_2s
    n_3s_prev <- n_3s
  }
  
  out$triangle <- if (length(t_tri) >= 2) diff(t_tri) else numeric()
  out$star2 <- if (length(t_2s) >= 2) diff(t_2s) else numeric()
  out$star3 <- if (length(t_3s) >= 2) diff(t_3s) else numeric()
  out
}

# ERNM-based approach (optimized with incremental change stats)
waiting_times_ernm <- function(net, time_attr = "time") {
  out <- list(triangle = numeric(), star2 = numeric(), star3 = numeric())
  el <- network::as.edgelist(net)
  if (nrow(el) == 0) return(out)
  
  edge_times <- NULL
  if (time_attr %in% network::list.edge.attributes(net)) {
    edge_times <- network::get.edge.attribute(net, time_attr)
  }
  if (is.null(edge_times) || length(edge_times) != nrow(el)) {
    vtimes <- network::get.vertex.attribute(net, time_attr)
    if (is.null(vtimes) || all(is.na(vtimes))) return(out)
    edge_times <- pmax(vtimes[el[, 1]], vtimes[el[, 2]], na.rm = TRUE)
  }
  
  ord <- order(edge_times)
  el <- el[ord, , drop = FALSE]
  edge_times <- edge_times[ord]
  event_times <- unique(edge_times)
  n <- network::network.size(net)
  g <- network::network.initialize(n, directed = network::is.directed(net))
  
  # Create ERNM model once and initialize with empty network
  model <- ernm::createCppModel(g ~ triangles + star(c(2,3)))
  model$setNetwork(ernm::as.BinaryNet(g))
  model$calculate()
  
  t_tri <- numeric()
  t_2s <- numeric()
  t_3s <- numeric()
  n_tri_prev <- 0L
  n_2s_prev <- 0L
  n_3s_prev <- 0L
  
  for (t_cur in event_times) {
    idx <- which(edge_times == t_cur)
    new_edges <- el[idx, , drop = FALSE]
    
    if (nrow(new_edges) > 0) {
      # Optimized: use computeChangeStats on batch BEFORE adding edges
      # This computes the change for all edges at once (faster than one-by-one)
      # Note: computeChangeStats handles edge interactions correctly when called on a batch
      change_stats <- model$computeChangeStats(new_edges[, 1], new_edges[, 2])
      # change_stats is a matrix: rows = edges, columns = [triangles, star(2), star(3)]
      delta_tri <- sum(change_stats[, 1])
      delta_2s <- sum(change_stats[, 2])
      delta_3s <- sum(change_stats[, 3])
      
      # Add all edges to network
      for (j in seq_len(nrow(new_edges))) {
        network::add.edges(g, tail = new_edges[j, 1], head = new_edges[j, 2])
      }
      
      # Update model once after adding all edges (much faster than per-edge updates)
      model$setNetwork(ernm::as.BinaryNet(g))
      model$calculate()
      
      # Update statistics incrementally using computed changes
      n_tri <- n_tri_prev + as.integer(round(delta_tri))
      n_2s <- n_2s_prev + as.integer(round(delta_2s))
      n_3s <- n_3s_prev + as.integer(round(delta_3s))
    } else {
      n_tri <- n_tri_prev
      n_2s <- n_2s_prev
      n_3s <- n_3s_prev
    }
    
    if (n_tri > n_tri_prev) t_tri <- c(t_tri, t_cur)
    if (n_2s > n_2s_prev) t_2s <- c(t_2s, t_cur)
    if (n_3s > n_3s_prev) t_3s <- c(t_3s, t_cur)
    
    n_tri_prev <- n_tri
    n_2s_prev <- n_2s
    n_3s_prev <- n_3s
  }
  
  out$triangle <- if (length(t_tri) >= 2) diff(t_tri) else numeric()
  out$star2 <- if (length(t_2s) >= 2) diff(t_2s) else numeric()
  out$star3 <- if (length(t_3s) >= 2) diff(t_3s) else numeric()
  out
}

# Benchmark both approaches
cat("=== Benchmarking ===\n")
cat("Current approach (manual incremental)...\n")
t_current <- proc.time()
result_current <- waiting_times_current(test_net)
elapsed_current <- (proc.time() - t_current)[3]
cat("  Time:", round(elapsed_current, 3), "s\n")
cat("  Triangles:", length(result_current$triangle), "| 2-stars:", length(result_current$star2), 
    "| 3-stars:", length(result_current$star3), "\n\n")

cat("ERNM approach...\n")
t_ernm <- proc.time()
result_ernm <- waiting_times_ernm(test_net)
elapsed_ernm <- (proc.time() - t_ernm)[3]
cat("  Time:", round(elapsed_ernm, 3), "s\n")
cat("  Triangles:", length(result_ernm$triangle), "| 2-stars:", length(result_ernm$star2), 
    "| 3-stars:", length(result_ernm$star3), "\n\n")

cat("=== Results ===\n")
cat("  Speedup:", round(elapsed_current / elapsed_ernm, 2), "x\n")
cat("  Results match:", 
    all.equal(result_current$triangle, result_ernm$triangle, tolerance = 1e-6),
    "&",
    all.equal(result_current$star2, result_ernm$star2, tolerance = 1e-6),
    "&",
    all.equal(result_current$star3, result_ernm$star3, tolerance = 1e-6), "\n")
