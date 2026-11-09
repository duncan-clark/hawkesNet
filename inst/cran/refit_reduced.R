# From the workspace root:
# Rscript hawkesNet/inst/cran/refit_reduced.R OUTPUT LAW [train|full|both|check]
# Retain archived fits; OUTPUT must be a new directory for this specification.
args <- commandArgs(TRUE)
if(length(args) < 2L) stop("Supply output directory and law")
out <- args[1]
law <- match.arg(args[2], c("CS-1", "CS-pair", "CS-pair-migration-add", "TERGM-logdegree"))
mode <- if(length(args) < 3L) "both" else match.arg(args[3], c("train", "full", "both", "check"))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
pkgload::load_all("hawkesNet", quiet=TRUE, debug=FALSE)
root <- "paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
diagnostics <- "paper_1/6_full_rewrite/diagnostics"
original <- file.path(diagnostics, "exponential_rates_20261105/fits")
archive <- switch(law,
  `CS-pair-migration-add`=file.path(diagnostics, "migration_addition_comparison_20261106"),
  `TERGM-logdegree`=file.path(diagnostics, "tergm_logdegree_20261107"), original)
m <- cran_update_data(root, law, term_specification="reduced")
m$release_fits <- readRDS(file.path(original, "release_fits.rds"))
m$release_values <- lapply(m$release_fits, "[[", "values")
fixed <- m$lower == m$upper
expected_count <- c(`CS-1`=33L, `CS-pair`=36L, `CS-pair-migration-add`=42L, `TERGM-logdegree`=33L)[[law]]
stopifnot(sum(!fixed) + 15L == expected_count,
          all(m$lower[grepl("cyclic$", m$names)] == 0),
          all(m$upper[grepl("cyclic$", m$names)] == 0))
if(law %in% c("CS-pair", "CS-pair-migration-add"))
  stopifnot(m$lower[match("migration", m$names)] == 0,
            m$upper[match("migration", m$names)] == 0)
windows <- if(mode %in% c("both", "check")) c("train", "full") else mode
for(window in windows) {
  name <- paste0(law, "_", window)
  old <- readRDS(file.path(archive, paste0(name, ".rds")))
  stopifnot(identical(names(old$parameters), m$names))
  # The unchanged native evaluator still reproduces the archived law exactly.
  old_eval <- loglik_cran_updates(m, old$parameters, TRUE)
  reproduction_error <- max(abs(old_eval$loglik - old$contributions$edge))
  stopifnot(reproduction_error < 1e-7)
  start <- pmin(m$upper, pmax(m$lower, old$parameters))
  direction <- sin(seq_along(start)) * (!fixed) / pmax(1, abs(start))
  h <- 1e-7
  evaluation <- loglik_cran_updates(m, start, TRUE)
  score <- sum(colSums(evaluation$gradient) * direction)
  difference <- sum(loglik_cran_updates(m, start+h*direction, TRUE)$loglik -
                    loglik_cran_updates(m, start-h*direction, TRUE)$loglik)/(2*h)
  gradient_error <- abs(score-difference)/max(1, abs(score), abs(difference))
  message(name, ": ", expected_count, " parameters; archive error ", signif(reproduction_error,3),
          "; reduced gradient error ", signif(gradient_error,3))
  stopifnot(gradient_error < 1e-5, all(is.finite(evaluation$loglik)))
  record <- list(law=law, window=window, term_specification="reduced",
    fixed_parameters=m$names[fixed], parameter_count=expected_count,
    archive=file.path(archive, paste0(name, ".rds")),
    archive_reproduction_error=reproduction_error, directional_gradient_error=gradient_error,
    in_stars="Both retained; no change to in-star terms", status="validated")
  jsonlite::write_json(record, file.path(out, paste0(name, "_validation.json")),
                       pretty=TRUE, auto_unbox=TRUE, digits=16)
  if(mode == "check") next
  dest <- file.path(out, paste0(name, ".rds"))
  if(file.exists(dest)) {
    fit <- readRDS(dest)
    stopifnot(identical(fit$term_specification,"reduced"), fit$parameter_count == expected_count)
    next
  }
  checkpoint <- file.path(out, paste0(name, "_checkpoint.rds"))
  alternative <- start + 0.05 * direction
  if(file.exists(checkpoint)) start <- readRDS(checkpoint)
  fit <- fit_cran_updates(m, window, parameters=list(start,alternative),
                         checkpoint=checkpoint, control=list(polish_tol=2e-6))
  stopifnot(fit$parameter_count == expected_count, all(fit$parameters[fixed] == m$lower[fixed]))
  saveRDS(fit, dest)
  record$status <- "fitted"
  jsonlite::write_json(c(record, list(loglik=fit$update_loglik,
    predictive=fit$predictive, parameters=as.list(fit$parameters),
    bounds=fit$bound_parameters, projected_gradient=fit$projected_gradient,
    multistart_loglik=fit$multistart_loglik)),
    file.path(out, paste0(name, ".json")), pretty=TRUE, auto_unbox=TRUE, digits=16)
}
