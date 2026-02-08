# Speed test for GOF function using OpenAlex network
# Run from package root: Rscript inst/openalex_study/test_gof_speed.R

devtools::load_all(".")
library(ernm)
library(network)
library(sna)
library(parallel)

# Paths: run from package root
PKG_ROOT <- getwd()
source(file.path(PKG_ROOT, "inst", "openalex_study", "get_network_openalex.R"))

# =============================================================================
# Config
# =============================================================================
EMAIL <- Sys.getenv("OPENALEX_EMAIL", "duncan-clark@outlook.com")
SEARCH_STRING <- Sys.getenv("OPENALEX_STRING", "Hawkes")
PAGES <- as.integer(Sys.getenv("OPENALEX_PAGES", 100))  # Full OpenAlex network
PER_PAGE <- 100L
MIN_DATE <- "1971-04-01"
MAX_DATE <- "2020-01-01"
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 7))  # Use available cores
FORMULA_RHS <- "edges + triangles + star(c(2,3))"
N_SIM <- 2L  # Just 2 sims for speed test

# =============================================================================
# Load and process OpenAlex network
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
# Set up initial parameters and fit (mock fit for testing)
# =============================================================================
cat("=== Setting up parameters and mock fit ===\n")
time_window_01 <- c(0, 1)

# Prepare inhomogeneous background
cat("  Preparing inhomogeneous background (KDE)...\n")
t_kde <- proc.time()
inhom_bg <- tryCatch(
  prepare_inhomogeneous_background(net_raw, time_attr = "time", bw = NULL, grid_n = 2048),
  error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
)
elapsed_kde <- (proc.time() - t_kde)[3]
cat("  KDE background:", round(elapsed_kde, 1), "s\n")

# Set up initial parameters
exp_cs <- expected_params_PMF_mark_CS(net_raw, FORMULA_RHS)
n_cs <- if (!is.na(exp_cs$CS_params_length)) exp_cs$CS_params_length else 4L
params_init_inhom <- list(
  mu = 1,
  beta_overall = 1,
  K = 0.5,
  beta_edges = 1,
  node_lambda = 1,
  CS_params = c(-10, rep(0, n_cs - 1)),
  vertex_categorical = list(gender = c(female = 0.1, male = 0.5)),
  vertex_categorical_levels = list(gender = c("female", "male", "unknown"))
)

# Create a mock fit object (for testing - normally this would come from actual fitting)
# We'll use the initial parameters as "fitted" parameters for this speed test
mock_fit <- list(
  fit = list(
    par = unlist(params_init_inhom[names(params_init_inhom) != "vertex_categorical_levels"])
  )
)

# Get cond_intensity closure (needed for simulation)
cat("  Creating cond_intensity closure...\n")
t_closure <- proc.time()
# Use regular cond_intensity for speed test (inhomogeneous needs special handling)
cond_intensity_func <- cond_intensity
elapsed_closure <- (proc.time() - t_closure)[3]
cat("  Closure creation:", round(elapsed_closure, 1), "s\n\n")

# =============================================================================
# Test GOF function with timing breakdown
# =============================================================================
cat("=== Testing GOF Function ===\n")
cat("  Simulations:", N_SIM, "\n")
cat("  Cores:", N_CORES, "\n")
cat("  Formula:", FORMULA_RHS, "\n\n")

t_gof_total <- proc.time()

# Break down the GOF function to time individual parts
cat("--- GOF Step 1: Parameter reconstruction ---\n")
t_step <- proc.time()
skel <- params_init_inhom
skel$vertex_categorical_levels <- NULL
pfit <- relist(mock_fit$fit$par, skeleton = skel)
pfit$vertex_categorical_levels <- params_init_inhom$vertex_categorical_levels
pfit$K <- params_init_inhom$K
if (!is.null(inhom_bg)) {
  Tval <- time_window_01[2] - time_window_01[1]
  pfit$mu <- inhom_bg$integral_bg / Tval
} else {
  pfit$mu <- params_init_inhom$mu
}
pfit$K <- min(max(pfit$K, 0.001), 0.999)
pfit$mu <- max(pfit$mu, 0.001)
pfit$node_lambda <- max(pfit$node_lambda, 1)
pfit$beta_overall <- max(pfit$beta_overall, 0.001)
pfit$beta_edges <- max(pfit$beta_edges, 0.001)
if (!is.null(params_init_inhom$vertex_categorical)) {
  if (is.null(pfit$vertex_categorical)) {
    pfit$vertex_categorical <- params_init_inhom$vertex_categorical
  }
}
elapsed_param <- (proc.time() - t_step)[3]
cat("  Parameter reconstruction:", round(elapsed_param, 1), "s\n\n")

cat("--- GOF Step 2: Simulations ---\n")
t_sim <- proc.time()
sim_results <- parallel::mclapply(seq_len(N_SIM), function(i) {
  t_one <- proc.time()
  s <- tryCatch(
    sim_hawkesGrowthNet(
      params = pfit,
      time_window = c(0, 0.05),
      PMF_mark = PMF_mark_CS,
      cond_intensity = cond_intensity_func,
      formula_RHS = FORMULA_RHS,
      truncation = 100L,
      mark_decay = "activity",
      max_node_time = 1,
      hashed_edges = TRUE,
      verbose = FALSE,
      mu_multiplier = 5,
      stop_on_full_network = FALSE
    ),
    error = function(e) { return(list(net = NULL, error = e$message)) }
  )
  elapsed_one <- (proc.time() - t_one)[3]
  if (!is.null(s$net)) {
    cat("    Sim", i, ":", network::network.size(s$net), "nodes,",
        network::network.edgecount(s$net), "edges (", round(elapsed_one, 1), "s)\n")
  } else {
    cat("    Sim", i, "FAILED:", ifelse(is.null(s$error), "unknown", s$error), "\n")
  }
  if (is.null(s$net)) {
    return(list(net = NULL, error = ifelse(is.null(s$error), "unknown", s$error)))
  }
  return(list(net = s$net, error = NULL))
}, mc.cores = N_CORES)

sim_nets <- lapply(sim_results, function(x) x$net)
sim_nets <- sim_nets[!sapply(sim_nets, is.null)]
n_success <- length(sim_nets)
elapsed_sim <- (proc.time() - t_sim)[3]
cat("  Simulations:", round(elapsed_sim, 1), "s (", n_success, "succeeded,", N_SIM - n_success, "failed)\n\n")

if (length(sim_nets) > 0) {
  cat("--- GOF Step 3: Computing statistics ---\n")
  t_stats <- proc.time()
  max_deg <- 15
  k_esp <- 15
  
  cat("  Computing observed statistics...\n")
  t_obs <- proc.time()
  # Access helper functions from gof.R (they're internal but available via package)
  degree_obs <- hawkesGrowthNet:::degree_dist(net_raw, max_deg)
  esp_obs <- hawkesGrowthNet:::esp_dist(net_raw, k_esp)
  geodist_obs <- hawkesGrowthNet:::geodist_dist(net_raw)
  wait_obs <- hawkesGrowthNet:::waiting_times_between_formations(net_raw, formula_RHS = FORMULA_RHS)
  elapsed_obs <- (proc.time() - t_obs)[3]
  cat("    Observed stats:", round(elapsed_obs, 1), "s\n")
  
  cat("  Computing simulated statistics (parallelized)...\n")
  t_sim_stats <- proc.time()
  
  cat("    Degree distributions...\n")
  t_deg <- proc.time()
  degree_sim <- do.call(rbind, parallel::mclapply(sim_nets, function(n) {
    hawkesGrowthNet:::degree_dist(n, max_deg)
  }, mc.cores = N_CORES))
  elapsed_deg <- (proc.time() - t_deg)[3]
  cat("      Degree:", round(elapsed_deg, 1), "s\n")
  
  cat("    ESP distributions...\n")
  t_esp <- proc.time()
  esp_sim <- do.call(rbind, parallel::mclapply(sim_nets, function(n) {
    hawkesGrowthNet:::esp_dist(n, k_esp)
  }, mc.cores = N_CORES))
  elapsed_esp <- (proc.time() - t_esp)[3]
  cat("      ESP:", round(elapsed_esp, 1), "s\n")
  
  cat("    Geodesic distances...\n")
  t_geod <- proc.time()
  geodist_sim <- parallel::mclapply(sim_nets, function(n) {
    hawkesGrowthNet:::geodist_dist(n)
  }, mc.cores = N_CORES)
  elapsed_geod <- (proc.time() - t_geod)[3]
  cat("      Geodesic:", round(elapsed_geod, 1), "s\n")
  
  cat("    Waiting times...\n")
  t_wait <- proc.time()
  wait_sim <- parallel::mclapply(sim_nets, function(n) {
    hawkesGrowthNet:::waiting_times_between_formations(n, formula_RHS = FORMULA_RHS)
  }, mc.cores = N_CORES)
  elapsed_wait <- (proc.time() - t_wait)[3]
  cat("      Waiting times:", round(elapsed_wait, 1), "s\n")
  
  elapsed_sim_stats <- (proc.time() - t_sim_stats)[3]
  elapsed_stats <- (proc.time() - t_stats)[3]
  cat("  Simulated stats:", round(elapsed_sim_stats, 1), "s\n")
  cat("  Total statistics:", round(elapsed_stats, 1), "s\n\n")
}

elapsed_gof_total <- (proc.time() - t_gof_total)[3]

cat("=== GOF Speed Test Complete ===\n")
cat("  Total GOF time:", round(elapsed_gof_total, 1), "s\n")
cat("  Breakdown:\n")
cat("    Parameter reconstruction:", round(elapsed_param, 1), "s\n")
cat("    Simulations:", round(elapsed_sim, 1), "s\n")
if (length(sim_nets) > 0) {
  cat("    Statistics computation:", round(elapsed_stats, 1), "s\n")
  cat("      - Observed:", round(elapsed_obs, 1), "s\n")
  cat("      - Simulated:", round(elapsed_sim_stats, 1), "s\n")
  cat("        - Degree:", round(elapsed_deg, 1), "s\n")
  cat("        - ESP:", round(elapsed_esp, 1), "s\n")
  cat("        - Geodesic:", round(elapsed_geod, 1), "s\n")
  cat("        - Waiting times:", round(elapsed_wait, 1), "s\n")
}
