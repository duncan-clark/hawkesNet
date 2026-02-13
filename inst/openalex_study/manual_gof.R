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
FORMULA_STRUCT <- "edges + gwdegree(0.1)"
FORMULA_STRUCT_GW001 <- "edges + gwdegree(0.01)"
FORMULA_STRUCT_DEGSTARS <- "degree(0) + star(c(2,3,4,5))"
FORMULA_MATCH  <- "edges + gwdegree(0.1) + nodeMatch('gender')"

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
exp_struct <- expected_params_PMF_mark_CS(net_raw, FORMULA_STRUCT)
n_cs_struct <- if (!is.na(exp_struct$CS_params_length)) exp_struct$CS_params_length else 2L
fit_struct <- fit_hawkesNet(
  params_init = make_p(n_cs_struct, mu_init),
  time_window = c(0, 1), mark_filtration = net_raw, PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec, integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_STRUCT, truncation = TRUNCATION, growth_only = GROWTH_ONLY,
  fixed_params = c("K", "mu"), cores = N_CORES, combine_intensity = TRUE,
  maxit = 1000
)
print(fit_struct$fit_table)

cat("\n--- Fitting Structural Model (gwdegree(0.01)) ---\n")
exp_struct_gw001 <- expected_params_PMF_mark_CS(net_raw, FORMULA_STRUCT_GW001)
n_cs_struct_gw001 <- if (!is.na(exp_struct_gw001$CS_params_length)) exp_struct_gw001$CS_params_length else 2L
fit_struct_gw001 <- fit_hawkesNet(
  params_init = make_p(n_cs_struct_gw001, mu_init),
  time_window = c(0, 1), mark_filtration = net_raw, PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec, integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_STRUCT_GW001, truncation = TRUNCATION, growth_only = GROWTH_ONLY,
  fixed_params = c("K", "mu"), cores = N_CORES, combine_intensity = TRUE,
  maxit = 1000
)
print(fit_struct_gw001$fit_table)

cat("\n--- Fitting Structural Model (degree(0) + star(2:5)) ---\n")
exp_struct_degstars <- expected_params_PMF_mark_CS(net_raw, FORMULA_STRUCT_DEGSTARS)
n_cs_struct_degstars <- if (!is.na(exp_struct_degstars$CS_params_length)) exp_struct_degstars$CS_params_length else 5L
fit_struct_degstars <- fit_hawkesNet(
  params_init = make_p(n_cs_struct_degstars, mu_init),
  time_window = c(0, 1), mark_filtration = net_raw, PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec, integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_STRUCT_DEGSTARS, truncation = TRUNCATION, growth_only = GROWTH_ONLY,
  fixed_params = c("K", "mu"), cores = N_CORES, combine_intensity = TRUE,
  maxit = 1000
)
print(fit_struct_degstars$fit_table)

cat("\n--- Fitting NodeMatch Model ---\n")
exp_match <- expected_params_PMF_mark_CS(net_raw, FORMULA_MATCH)
n_cs_match <- if (!is.na(exp_match$CS_params_length)) exp_match$CS_params_length else 3L
fit_match <- fit_hawkesNet(
  params_init = make_p(n_cs_match, mu_init, gender = TRUE),
  time_window = c(0, 1), mark_filtration = net_raw, PMF_mark = PMF_mark_CS,
  mu_vec = inhom_bg$mu_vec, integral_bg = inhom_bg$integral_bg,
  formula_RHS = FORMULA_MATCH, truncation = TRUNCATION, growth_only = GROWTH_ONLY,
  fixed_params = c("K", "mu"), cores = N_CORES, combine_intensity = TRUE,
  maxit = 1000
)
print(fit_match$fit_table)

# --- 3. GOF ---
cat("\n--- Running GOF (Structural) ---\n")
gof_struct <- gof(fit = fit_struct, net_obs = net_raw, params_init = make_p(n_cs_struct, mu_init),
                  PMF_mark = PMF_mark_CS, cond_intensity = cond_intensity,
                  formula_RHS = FORMULA_STRUCT, time_window = c(0, 1), 
                  inhom_bg = inhom_bg, n_sim = 2, cores = N_CORES, 
                  seed_events = 20, growth_only = GROWTH_ONLY)

cat("\n--- Running GOF (Structural gwdegree(0.01)) ---\n")
gof_struct_gw001 <- gof(fit = fit_struct_gw001, net_obs = net_raw, params_init = make_p(n_cs_struct_gw001, mu_init),
                 PMF_mark = PMF_mark_CS, cond_intensity = cond_intensity,
                 formula_RHS = FORMULA_STRUCT_GW001, time_window = c(0, 1),
                 inhom_bg = inhom_bg, n_sim = 2, cores = N_CORES,
                 seed_events = 20, growth_only = GROWTH_ONLY)

cat("\n--- Running GOF (Structural degree+stars) ---\n")
gof_struct_degstars <- gof(fit = fit_struct_degstars, net_obs = net_raw, params_init = make_p(n_cs_struct_degstars, mu_init),
                    PMF_mark = PMF_mark_CS, cond_intensity = cond_intensity,
                    formula_RHS = FORMULA_STRUCT_DEGSTARS, time_window = c(0, 1),
                    inhom_bg = inhom_bg, n_sim = 2, cores = N_CORES,
                    seed_events = 20, growth_only = GROWTH_ONLY)

cat("\n--- Running GOF (NodeMatch) ---\n")
gof_match <- gof(fit = fit_match, net_obs = net_raw, params_init = make_p(n_cs_match, mu_init, gender=T),
                 PMF_mark = PMF_mark_CS, cond_intensity = cond_intensity,
                 formula_RHS = FORMULA_MATCH, time_window = c(0, 1), 
                 inhom_bg = inhom_bg, n_sim = 2, cores = N_CORES, 
                 seed_events = 20, growth_only = GROWTH_ONLY)

cat("\n--- Fitting BA Model ---\n")
fit_ba <- fit_hawkesNet(
  params_init = list(mu = mu_init, beta_overall = 1, K = 0.5, beta_edges = 1, m = 1),
  time_window = c(0, 1), mark_filtration = net_raw, PMF_mark = PMF_mark_BA,
  mu_vec = inhom_bg$mu_vec, integral_bg = inhom_bg$integral_bg,
  truncation = TRUNCATION, fixed_params = c("mu"), cores = N_CORES,
  maxit = 1000
)
print(fit_ba$fit_table)

cat("\n--- Running GOF (BA) ---\n")
gof_ba <- gof(fit = fit_ba, net_obs = net_raw, 
              params_init = list(mu = mu_init, beta_overall = 1, K = 0.5, beta_edges = 1, m = 1),
              PMF_mark = PMF_mark_BA, cond_intensity = cond_intensity,
              formula_RHS = "",
              time_window = c(0, 1), inhom_bg = inhom_bg, n_sim = 2, 
              cores = N_CORES, seed_events = 20)

# Display plots if ggplot2 is available
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  if (!is.null(gof_struct$plots$degree_plot)) print(gof_struct$plots$degree_plot + labs(subtitle = "Structural"))
  if (!is.null(gof_struct_gw001$plots$degree_plot)) print(gof_struct_gw001$plots$degree_plot + labs(subtitle = "Structural gwdegree(0.01)"))
  if (!is.null(gof_struct_degstars$plots$degree_plot)) print(gof_struct_degstars$plots$degree_plot + labs(subtitle = "Structural degree+stars"))
  if (!is.null(gof_match$plots$degree_plot))  print(gof_match$plots$degree_plot + labs(subtitle = "NodeMatch"))
  if (!is.null(gof_ba$plots$degree_plot))     print(gof_ba$plots$degree_plot + labs(subtitle = "BA"))
}

cat("\nDone. Objects 'fit_struct', 'fit_struct_gw001', 'fit_struct_degstars', 'fit_match', 'fit_ba', 'gof_struct', 'gof_struct_gw001', 'gof_struct_degstars', 'gof_match', 'gof_ba', 'net_raw' available.\n")
