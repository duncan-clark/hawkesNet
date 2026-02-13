# =============================================================================
# Minimal OpenAlex Hawkes Study & Explorer
# =============================================================================
# Run interactively: source("inst/openalex_study/manual_gof.R")
# =============================================================================

library(hawkesNet)
library(parallel)
# Load source files directly for verification to ensure we use the latest changes
source("R/hawkesNet.R")
source("R/mark_PMF.R")
source("R/gof.R")
source("R/utils.R")
source("R/temporal_hawkes.R")
source("R/kde_background.R")

library(network)
library(sna)
library(ernm)
library(dplyr)

# --- 1. Fetch Data ---
EMAIL <- "duncan-clark@outlook.com"
SEARCH_STRING <- "Hawkes Process"
TOPIC <- "Point processes and geometric inequalities"
PAGES <- 2
N_CORES <- 1 # Use 1 core locally to avoid complex PSOCK exports in minimal script

cat("--- Fetching OpenAlex Data ---\n")
source("inst/openalex_study/get_network_openalex.R")
out <- get_network(email = EMAIL, pages = PAGES, per_page = 100L,
                   string = SEARCH_STRING, min_date = "1971-04-01", 
                   max_date = "2020-01-01", topics_include = c(TOPIC))

net_raw <- out$net
set.vertex.attribute(net_raw, "time", net_raw %v% "time_scaled")
set.edge.attribute(net_raw, "time", net_raw %e% "time_scaled")
net_raw <- normalize_times_01(net_raw, attr = "time", keep_na = TRUE)

cat(sprintf("Network: %d events, %d nodes\n", 
            length(get_times(net_raw)$times), network.size(net_raw)))

# --- 2. Fit Models ---
TRUNCATION <- 300
GROWTH_ONLY <- TRUE
FORMULA_STRUCT <- "edges + triangles + gwdegree(0.5)"
FORMULA_MATCH  <- "edges + triangles + gwdegree(0.5) + nodeMatch('gender')"

cat("\n--- Preparing Background ---\n")
inhom_bg <- prepare_inhomogeneous_background(net_raw, time_attr = "time")
mu_init <- inhom_bg$integral_bg

# Helper for params
make_p <- function(n_cs, mu, gender = FALSE) {
  p <- list(mu = mu, beta_overall = 1, K = 0.5, beta_edges = 1, 
            node_lambda = 1, CS_params = c(-10, rep(0, n_cs - 1)))
  if (gender) {
    p$vertex_categorical <- list(gender = c(female = 0.1, male = 0.5))
    p$vertex_categorical_levels <- list(gender = c("female", "male", "unknown"))
  }
  p
}

cat("\n--- Fitting Structural Model ---\n")
fit_struct <- fit_hawkesNet(
  params_init = make_p(3, mu_init),
  time_window = c(0, 1), mark_filtration = net_raw, PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec, integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_STRUCT, truncation = TRUNCATION, growth_only = GROWTH_ONLY,
  fixed_params = c("K", "mu"), cores = N_CORES, combine_intensity = TRUE
)
print(fit_struct$fit_table)

cat("\n--- Fitting NodeMatch Model ---\n")
fit_match <- fit_hawkesNet(
  params_init = make_p(4, mu_init, gender = TRUE),
  time_window = c(0, 1), mark_filtration = net_raw, PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec, integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_MATCH, truncation = TRUNCATION, growth_only = GROWTH_ONLY,
  fixed_params = c("K", "mu"), cores = N_CORES, combine_intensity = TRUE
)
print(fit_match$fit_table)

# --- 3. GOF ---
cat("\n--- Running GOF (Structural) ---\n")
gof_struct <- gof(fit = fit_struct, net_obs = net_raw, params_init = make_p(3, mu_init),
                  PMF_mark = PMF_mark_CS, cond_intensity = cond_intensity,
                  formula_RHS = FORMULA_STRUCT, time_window = c(0, 1), 
                  inhom_bg = inhom_bg, n_sim = 2, cores = N_CORES, 
                  seed_events = 20, growth_only = GROWTH_ONLY)

cat("\n--- Running GOF (NodeMatch) ---\n")
gof_match <- gof(fit = fit_match, net_obs = net_raw, params_init = make_p(4, mu_init, gender=T),
                 PMF_mark = PMF_mark_CS, cond_intensity = cond_intensity,
                 formula_RHS = FORMULA_MATCH, time_window = c(0, 1), 
                 inhom_bg = inhom_bg, n_sim = 2, cores = N_CORES, 
                 seed_events = 20, growth_only = GROWTH_ONLY)

cat("\nDone. Objects 'fit_struct', 'fit_match', 'gof_struct', 'gof_match', 'net_raw' available.\n")
