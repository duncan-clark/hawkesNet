# Test script: Benchmark waiting_times_between_formations on OpenAlex data
# Run from package root: Rscript inst/openalex_study/test_waiting_times_openalex.R

devtools::load_all(".")
library(ernm)
library(network)
library(sna)

# Paths: run from package root
PKG_ROOT <- getwd()
source(file.path(PKG_ROOT, "inst", "openalex_study", "get_network_openalex.R"))

# =============================================================================
# Config
# =============================================================================
EMAIL <- Sys.getenv("OPENALEX_EMAIL", "duncan-clark@outlook.com")
SEARCH_STRING <- Sys.getenv("OPENALEX_STRING", "Hawkes")
PAGES <- as.integer(Sys.getenv("OPENALEX_PAGES", 100))
PER_PAGE <- 100L
MIN_DATE <- "1971-04-01"
MAX_DATE <- "2020-01-01"
FORMULA_RHS <- "edges + triangles + star(c(2,3))"

# =============================================================================
# Load OpenAlex network
# =============================================================================
cat("=== Loading OpenAlex Network ===\n")
cat("  Search:", SEARCH_STRING, "| Pages:", PAGES, "\n")
cat("  Date range:", MIN_DATE, "to", MAX_DATE, "\n\n")

t_load <- proc.time()
out <- get_network(
  email = EMAIL,
  string = SEARCH_STRING,
  pages = PAGES,
  per_page = PER_PAGE,
  min_date = MIN_DATE,
  max_date = MAX_DATE
)
net_raw <- out$net
edges <- out$edges

# Process network same as OpenAlex study script
network::set.vertex.attribute(net_raw, "time", net_raw %v% "time_scaled")
network::set.edge.attribute(net_raw, "time", net_raw %e% "time_scaled")
net_raw <- hawkesGrowthNet::normalize_times_01(net_raw, attr = "time", keep_na = TRUE)

elapsed_load <- (proc.time() - t_load)[3]

n_nodes <- network::network.size(net_raw)
n_edges <- network::network.edgecount(net_raw)
n_events <- length(hawkesGrowthNet::get_times(net_raw)$times)

cat("\nNetwork loaded:", round(elapsed_load, 1), "s\n")
cat("  Nodes:", n_nodes, "\n")
cat("  Edges:", n_edges, "\n")
cat("  Events:", n_events, "\n")
cat("  Average degree:", round(2 * n_edges / n_nodes, 2), "\n\n")

# =============================================================================
# Define waiting_times_between_formations function
# =============================================================================
waiting_times_between_formations <- function(net, time_attr = "time", 
                                             formula_RHS = "triangles + star(c(2,3))") {
  # Get edge list and edge times
  el <- network::as.edgelist(net)
  
  # Determine number of statistics by creating a temporary model
  n <- network::network.size(net)
  g_temp <- network::network.initialize(n, directed = network::is.directed(net))
  # Create formula with g_temp in the environment
  formula_str <- paste0("g_temp ~ ", formula_RHS)
  formula_obj <- as.formula(formula_str)
  model_temp <- ernm::createCppModel(formula_obj)
  model_temp$setNetwork(ernm::as.BinaryNet(g_temp))
  model_temp$calculate()
  n_stats <- length(model_temp$statistics())
  rm(model_temp, g_temp)
  
  if (nrow(el) == 0) {
    # Return empty structure matching formula
    return(setNames(rep(list(numeric()), n_stats), paste0("stat", seq_len(n_stats))))
  }

  # Try edge times first; fall back to vertex times for edge ordering
  edge_times <- NULL
  if (time_attr %in% network::list.edge.attributes(net)) {
    edge_times <- network::get.edge.attribute(net, time_attr)
  }
  if (is.null(edge_times) || length(edge_times) != nrow(el)) {
    # Use vertex times: assign each edge the max time of its endpoints
    vtimes <- network::get.vertex.attribute(net, time_attr)
    if (is.null(vtimes) || all(is.na(vtimes))) {
      return(setNames(rep(list(numeric()), n_stats), paste0("stat", seq_len(n_stats))))
    }
    edge_times <- pmax(vtimes[el[, 1]], vtimes[el[, 2]], na.rm = TRUE)
  }

  # Sort edges by time
  ord <- order(edge_times)
  el <- el[ord, , drop = FALSE]
  edge_times <- edge_times[ord]

  # Group edges by unique event times
  event_times <- unique(edge_times)
  g <- network::network.initialize(n, directed = network::is.directed(net))

  # Initialize tracking for all statistics
  stat_prev <- rep(0, n_stats)
  stat_times <- rep(list(numeric()), n_stats)

  # Create ERNM model once and initialize with empty network
  # Create formula with g now that it exists
  formula_str <- paste0("g ~ ", formula_RHS)
  formula_obj <- as.formula(formula_str)
  model <- ernm::createCppModel(formula_obj)
  model$setNetwork(ernm::as.BinaryNet(g))
  model$calculate()
  
  for (t_cur in event_times) {
    idx <- which(edge_times == t_cur)
    new_edges <- el[idx, , drop = FALSE]
    
    if (nrow(new_edges) > 0) {
      # Incremental approach: add edges one at a time, get change stats, update model
      # This correctly handles interactions (e.g., triangles formed by multiple edges)
      stat_curr <- stat_prev
      for (j in seq_len(nrow(new_edges))) {
        tail_j <- new_edges[j, 1]
        head_j <- new_edges[j, 2]
        
        # Get change stats BEFORE adding this edge (from current network state)
        change_stats_j <- model$computeChangeStats(tail_j, head_j)
        # change_stats_j is a matrix: rows = 1 (single edge), columns = statistics
        delta_stats <- change_stats_j[1, ]
        
        # Add edge to network
        network::add.edges(g, tail = tail_j, head = head_j)
        
        # Update model to reflect this new edge (for next iteration's change stats)
        model$setNetwork(ernm::as.BinaryNet(g))
        model$calculate()
        
        # Update statistics incrementally
        stat_curr <- stat_curr + as.integer(round(delta_stats))
      }
    } else {
      stat_curr <- stat_prev
    }
    
    # Track when each statistic increases
    for (i in seq_len(n_stats)) {
      if (stat_curr[i] > stat_prev[i]) {
        stat_times[[i]] <- c(stat_times[[i]], t_cur)
      }
    }
    
    stat_prev <- stat_curr
  }

  # Compute waiting times (differences between formation times)
  out <- lapply(stat_times, function(t_vec) {
    if (length(t_vec) >= 2) diff(t_vec) else numeric()
  })
  
  # Name outputs based on formula
  # For backward compatibility, use triangle, star2, star3 if formula matches
  if (formula_RHS == "triangles + star(c(2,3))") {
    names(out) <- c("triangle", "star2", "star3")
  } else if (formula_RHS == "edges + triangles + star(c(2,3))") {
    names(out) <- c("edges", "triangle", "star2", "star3")
  } else {
    # Generic naming: try to extract statistic names from formula
    # This is a simple heuristic - may not work for all formulas
    stat_names <- trimws(strsplit(formula_RHS, "\\+")[[1]])
    if (length(stat_names) == length(out)) {
      names(out) <- stat_names
    } else {
      names(out) <- paste0("stat", seq_len(length(out)))
    }
  }
  
  out
}

# =============================================================================
# Benchmark waiting_times_between_formations
# =============================================================================
cat("=== Benchmarking waiting_times_between_formations ===\n")
cat("  Formula:", FORMULA_RHS, "\n\n")

# Run benchmark
cat("  Computing waiting times...\n")
t_wait <- proc.time()
wait_results <- waiting_times_between_formations(net_raw, formula_RHS = FORMULA_RHS)
elapsed_wait <- (proc.time() - t_wait)[3]

cat("\nWaiting times computed:", round(elapsed_wait, 1), "s\n")
cat("  Results:\n")
for (stat_name in names(wait_results)) {
  n_formations <- length(wait_results[[stat_name]])
  if (n_formations > 0) {
    mean_wait <- mean(wait_results[[stat_name]])
    cat("    ", stat_name, ":", n_formations, "waiting times, mean =", 
        round(mean_wait, 4), "\n")
  } else {
    cat("    ", stat_name, ": no formations\n")
  }
}

cat("\n=== Benchmark Complete ===\n")
cat("  Total time:", round(elapsed_load + elapsed_wait, 1), "s\n")
cat("  Network load:", round(elapsed_load, 1), "s\n")
cat("  Waiting times:", round(elapsed_wait, 1), "s\n")
