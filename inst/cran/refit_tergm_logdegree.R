# Run from the project root: Rscript hawkesNet/inst/cran/refit_tergm_logdegree.R OUTPUT train|full
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply an output directory and a fitting window")
out <- args[[1L]]
window <- match.arg(args[[2L]], c("train", "full"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (dir.exists("hawkesNet")) pkgload::load_all("hawkesNet", quiet = TRUE, debug = FALSE) else library(hawkesNet)

root <- "paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
archive <- "paper_1/6_full_rewrite/diagnostics/exponential_rates_20261105/fits"
data <- cran_update_data(root, "TERGM-logdegree")
data$release_fits <- readRDS(file.path(archive, "release_fits.rds"))
data$release_values <- lapply(data$release_fits, "[[", "values")
old <- readRDS(file.path(archive, paste0("TERGM_", window, ".rds")))

# With the four new covariates set to zero, recover the archived likelihood.
embedded <- c(old$parameters, rep(0, 4))
reference <- loglik_cran_updates(data, embedded, TRUE)
stopifnot(max(abs(reference$loglik - old$contributions$edge)) < 1e-8)
start <- embedded
start[data$lower == data$upper] <- data$lower[data$lower == data$upper]
alternative <- start
alternative[20:23] <- c(0.4, 0.4, -0.4, -0.2)
direction <- cos(seq_along(start)) * (data$lower != data$upper)
h <- 1e-5
score <- sum(loglik_cran_updates(data, alternative)$gradient * direction)
plus <- loglik_cran_updates(data, alternative + h * direction, TRUE)$loglik
minus <- loglik_cran_updates(data, alternative - h * direction, TRUE)$loglik
difference <- sum(plus - minus) / (2 * h)
stopifnot(abs(score - difference) / max(1, abs(score), abs(difference)) < 1e-5)
message("Exact archive reproduction and directional score check passed. Starting ", window, " fit.")

name <- paste0("TERGM-logdegree_", window)
checkpoint <- file.path(out, paste0(name, "_checkpoint.rds"))
if (file.exists(checkpoint)) start <- readRDS(checkpoint)
fit <- fit_cran_updates(data, window, parameters = list(start, alternative), checkpoint = checkpoint)
stopifnot(fit$parameter_count == 34L, all(fit$parameters[data$lower == data$upper] == 0))
saveRDS(fit, file.path(out, paste0(name, ".rds")))
jsonlite::write_json(list(window = window, parameter_count = fit$parameter_count,
  loglik = fit$update_loglik, predictive = fit$predictive,
  parameters = as.list(fit$parameters), bound_parameters = fit$bound_parameters,
  fixed_parameters = fit$fixed_parameters, projected_gradient = fit$projected_gradient,
  multistart_loglik = fit$multistart_loglik),
  file.path(out, paste0(name, ".json")), auto_unbox = TRUE, pretty = TRUE, digits = 16)
