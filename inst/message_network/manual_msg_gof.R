# =============================================================================
# Message Network Hawkes Study (Modernized)
# =============================================================================
# Data: https://snap.stanford.edu/data/CollegeMsg.html
# =============================================================================

library(hawkesNet)
library(network)
library(sna)
library(ernm)
library(dplyr)
library(parallel)

# --- 1. Configuration ---
N_CORES <- 1 # Local run
TRUNCATION <- 100
GROWTH_ONLY <- FALSE # As requested
FORMULA_RHS <- "edges + triangles + star(2) + star(3)"
DATA_PATH <- "hawkesNet/data/CollegeMsg.txt"

# --- 2. Load and Preprocess Data ---
if (!file.exists(DATA_PATH)) {
  # Try to find it in the workspace
  DATA_PATH <- "data/CollegeMsg.txt"
}

if (!file.exists(DATA_PATH)) {
  stop("CollegeMsg.txt not found. Please ensure it is in 'data/' or 'hawkesNet/data/'.")
}

cat("--- Loading Message Network Data ---\n")
dat <- read.table(DATA_PATH)
names(dat) <- c('from', 'to', 'time')

# Convert to days and normalize
dat$time <- (dat$time - min(dat$time)) / (24 * 60 * 60)

# Make undirected (canonical ordering) and distinct
swap <- dat$from > dat$to
dat[swap, c("from", "to")] <- dat[swap, c("to", "from")]
dat <- dat %>% distinct()

# Get first interaction time for each pair
edges <- dat %>%
  group_by(from, to) %>%
  summarize(time = min(time), .groups = 'drop') %>%
  arrange(time)

# Get node arrival times
node_times <- bind_rows(
  edges %>% select(id = from, time),
  edges %>% select(id = to, time)
) %>%
  group_by(id) %>%
  summarize(time = min(time), .groups = 'drop') %>%
  arrange(time)

# Build network object
net <- as.network(edges %>% select(from, to), matrix.type = "edgelist", directed = FALSE)
set.edge.attribute(net, "time", edges$time)
set.vertex.attribute(net, "time", node_times$time)

# Use a subset for a quicker "manual" run if needed, or full
# Let's take the first 14 days as in the old script
T_LIMIT <- 14
net_sub <- filtration_to_net(net, T_LIMIT)
delete.vertices(net_sub, isolates(net_sub))

cat(sprintf("Sub-network (up to day %d): %d events, %d nodes\n", 
            T_LIMIT, length(get_times(net_sub)$times), network.size(net_sub)))

# --- 3. Fit Model ---
cat("\n--- Preparing Background ---\n")
inhom_bg <- prepare_inhomogeneous_background(net_sub, time_attr = "time")
mu_init <- inhom_bg$integral_bg

# Initial parameters: edges, triangles, star2, star3 (4 CS params)
params_init <- list(
  mu = mu_init,
  beta_overall = 1.0,
  K = 0.5,
  beta_edges = 1.0,
  node_lambda = 1.0,
  CS_params = c(-5, 0, 0, 0)
)

cat("\n--- Fitting Message Network Model ---\n")
fit_msg <- fit_hawkesNet(
  params_init = params_init,
  time_window = c(0, T_LIMIT),
  mark_filtration = net_sub,
  PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec,
  integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_RHS,
  truncation = TRUNCATION,
  growth_only = GROWTH_ONLY,
  fixed_params = c("K", "mu"),
  cores = N_CORES,
  combine_intensity = TRUE,
  maxit = 1000
)

print(fit_msg$fit_table)

# --- 4. GOF ---
cat("\n--- Running GOF ---\n")
# Note: n_sim=2 for a quick check
gof_msg <- gof(
  fit = fit_msg,
  net_obs = net_sub,
  params_init = params_init,
  PMF_mark = PMF_mark_CS,
  cond_intensity = cond_intensity,
  formula_RHS = FORMULA_RHS,
  time_window = c(0, T_LIMIT),
  inhom_bg = inhom_bg,
  n_sim = 2,
  cores = N_CORES,
  seed_events = 20,
  growth_only = GROWTH_ONLY
)

cat("\nDone. Results in 'fit_msg' and 'gof_msg'.\n")
