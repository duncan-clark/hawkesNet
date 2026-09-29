#!/usr/bin/env Rscript
# Small, reproducible simulation/refit checks for the revised mark laws.
# This is not a consistency study. No historical outputs are read or replaced.
# Example, from any working directory:
# Rscript /path/to/hawkesNet/inst/simulation_study/smoke_collapsed_poisson.R \
#   --output-dir=/tmp/hawkesnet-collapsed-smoke --models=BA

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(models = "BA", horizons = "2,4", replicates = "2",
                 seed = "20260921", maxit = "200", fit = "mark_scale")
opts <- defaults
for (arg in args) {
  if (!grepl("^--[^=]+=.+$", arg)) stop("Use --name=value arguments: ", arg)
  key <- sub("^--([^=]+)=.*$", "\\1", arg)
  if (!key %in% c(names(defaults), "output-dir")) stop("Unknown argument: ", key)
  opts[[key]] <- sub("^--[^=]+=", "", arg)
}
if (is.null(opts[["output-dir"]])) {
  stop("Supply --output-dir=/tmp/a-new-directory; historical study outputs must be preserved.")
}
models <- strsplit(opts$models, ",", fixed = TRUE)[[1L]]
horizons <- as.numeric(strsplit(opts$horizons, ",", fixed = TRUE)[[1L]])
n_reps <- as.integer(opts$replicates)
seed <- as.integer(opts$seed)
maxit <- as.integer(opts$maxit)
stopifnot(length(models) > 0L, all(models %in% c("BA", "CS")),
          !anyDuplicated(models), length(horizons) > 0L, !anyDuplicated(horizons),
          all(is.finite(horizons)), all(horizons > 0),
          length(n_reps) == 1L, !is.na(n_reps), n_reps > 0L,
          length(seed) == 1L, !is.na(seed), seed >= 0L,
          length(maxit) == 1L, !is.na(maxit), maxit > 0L,
          opts$fit %in% c("mark_scale", "joint"))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file using Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
pkg_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
for (dep in c("pkgload", "ernm", "network", "data.table", "hash", "sna")) {
  if (!requireNamespace(dep, quietly = TRUE)) stop("Missing dependency: ", dep)
}
suppressPackageStartupMessages(pkgload::load_all(pkg_root, quiet = TRUE, helpers = FALSE))

output_dir <- path.expand(opts[["output-dir"]])
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty; choose a new directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
cat("Source:", pkg_root, "\nOutput:", output_dir, "\n")

params_by_model <- list(
  BA = list(mu = 8, beta_overall = 1.25, K = 0.35, beta_edges = 0.6, m = 1.4),
  CS = list(mu = 8, beta_overall = 1.25, K = 0.35, beta_edges = 0.6,
            node_lambda = 1, m = 1.4, CS_params = c(0, 0.4, 0.02, -0.01))
)
mark_options <- list(
  BA = list(truncation = 20L, mark_decay = "node_entrance"),
  # CS normalizes an interacting law over complete unordered edge subsets.
  # Four candidate vertices keep its exact state enumeration deliberately small.
  CS = list(truncation = 4L, max_candidates = 12L, mark_decay = "node_entrance",
            growth_only = FALSE, formula_RHS = "edges + triangles + star(c(2,3))")
)
source_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
manifest <- list(
  purpose = "Bounded smoke/conditional recovery check; not evidence of consistency",
  options = opts, models = models, horizons = horizons, replicates = n_reps,
  params = params_by_model[models], mark_options = mark_options[models],
  ground_intensity = "mu + K * sum(exp(-beta_overall * (t - t_j)))",
  branching_ratios = vapply(params_by_model[models], function(p) p$K / p$beta_overall, numeric(1)),
  seed_rule = "base seed + model index in c(BA,CS) * 100000 + horizon index * 1000 + replicate",
  fit_scope = if (opts$fit == "mark_scale") "Only m estimated; all other parameters fixed at truth" else
    "Joint exploratory fit; CS edges coefficient fixed at its generating value",
  package_source = pkg_root, source_md5 = tools::md5sum(source_files),
  session = sessionInfo(), started_at = Sys.time()
)
saveRDS(manifest, file.path(output_dir, "manifest.rds"))
capture.output(str(manifest, max.level = 3L), file = file.path(output_dir, "manifest.txt"))

# Exhaustive normalization on a tiny BA update space; also exercise the cached
# density against direct evaluation. This would reject the former Poisson-count
# times unconditional-Bernoulli expression before fitting any simulated path.
checks <- data.frame(check = character(), value = double(), passed = logical())
if ("BA" %in% models) {
  base_net <- network::network.initialize(3L, directed = FALSE)
  network::add.edges(base_net, c(1L, 2L), c(2L, 3L))
  network::set.vertex.attribute(base_net, "time", c(0.1, 0.2, 0.3))
  network::set.edge.attribute(base_net, "time", c(0.2, 0.3))
  masks <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 3L)))
  pmfs <- numeric(nrow(masks))
  cache_errors <- numeric(nrow(masks))
  for (i in seq_len(nrow(masks))) {
    observed <- network::network.copy(base_net)
    network::add.vertices(observed, 1L)
    network::set.vertex.attribute(observed, "time", c(0.1, 0.2, 0.3, 1))
    heads <- which(masks[i, ])
    if (length(heads)) {
      network::add.edges(observed, heads, rep(4L, length(heads)))
      network::set.edge.attribute(observed, "time", c(0.2, 0.3, rep(1, length(heads))))
    }
    density <- PMF_mark_BA(time = 1, params = params_by_model$BA,
                          mark_filtration = observed, mark = observed)
    pmfs[i] <- exp(density$log_mark_density)
    cache_errors[i] <- abs(density$log_density_func(params_by_model$BA) - density$log_mark_density)
  }
  checks <- rbind(checks,
    data.frame(check = "BA tiny-space PMF sum", value = sum(pmfs),
               passed = isTRUE(abs(sum(pmfs) - 1) < 1e-10)),
    data.frame(check = "BA tiny-space cached/direct max error", value = max(cache_errors),
               passed = isTRUE(max(cache_errors) < 1e-10)))
}
write.csv(checks, file.path(output_dir, "checks.csv"), row.names = FALSE)
if (nrow(checks) && !all(checks$passed)) {
  print(checks)
  stop("Mark-law preflight failed; no recovery fits were run.")
}

attempts <- list()
estimates <- list()
for (model in models) for (hi in seq_along(horizons)) for (rep in seq_len(n_reps)) {
  horizon <- horizons[hi]
  task_seed <- seed + match(model, c("BA", "CS")) * 100000L + hi * 1000L + rep
  set.seed(task_seed)
  params_true <- params_by_model[[model]]
  params_init <- params_true
  # Fixed parameters are copied exactly from truth, never randomly perturbed.
  params_init$m <- 0.9
  if (opts$fit == "joint") {
    params_init$mu <- 6
    params_init$beta_overall <- 1
    params_init$K <- 0.25
    params_init$beta_edges <- 0.4
    if (model == "CS") params_init$CS_params[-1L] <- 0
  }
  fixed <- if (opts$fit == "mark_scale") setdiff(names(params_true), "m") else
    if (model == "CS") "CS_params1" else NULL
  pmf <- get(paste0("PMF_mark_", model), envir = asNamespace("hawkesNet"))
  extra <- mark_options[[model]]
  warnings_seen <- character()
  rec <- list(model = model, horizon = horizon, replicate = rep, seed = task_seed,
              fit_scope = opts$fit, events = NA_integer_, recovered_events = NA_integer_,
              nodes = NA_integer_, edges = NA_integer_, loglik_true = NA_real_,
              loglik_init = NA_real_, loglik_fitted = NA_real_, cache_error = NA_real_,
              optimizer_error = NA_real_, convergence = NA_integer_,
              boundary_parameters = "",
              simulation_seconds = NA_real_, fit_seconds = NA_real_,
              status = "error", error = "", warnings = "")
  result <- tryCatch(withCallingHandlers({
    cat(sprintf("%s T=%g replicate=%d seed=%d\n", model, horizon, rep, task_seed))
    tick <- proc.time()[["elapsed"]]
    sim <- do.call(sim_hawkesNet, c(list(params = params_true, time_window = c(0, horizon),
      PMF_mark = pmf, cond_intensity = cond_intensity, hashed_edges = TRUE,
      verbose = FALSE, mu_multiplier = 3, joint_accept = FALSE,
      stop_on_full_network = FALSE), extra))
    rec$simulation_seconds <- proc.time()[["elapsed"]] - tick
    rec$events <- length(sim$events$t)
    rec$recovered_events <- length(get_times(sim$net)$times)
    rec$nodes <- network::network.size(sim$net)
    rec$edges <- network::network.edgecount(sim$net)
    if (rec$events != rec$recovered_events) {
      stop("Simulated event times cannot all be recovered from the stored filtration: ",
           rec$events, " simulated versus ", rec$recovered_events, " recoverable.")
    }
    if (rec$events < 2L) stop("Too few events for a meaningful refit smoke check.")
    likelihood <- function(params, combined = FALSE, cache = NULL) {
      do.call(loglik_hawkesNet, c(list(params = params, time_window = c(0, horizon),
        mark_filtration = sim$net, PMF_mark = pmf, verbose = FALSE, cores = 1L,
        combine_intensity = combined, intens_funcs = cache), extra))
    }
    ll_true <- likelihood(params_true)
    ll_combined <- likelihood(params_true, combined = TRUE)
    ll_init <- likelihood(params_init, cache = ll_combined$intens_funcs)
    rec$loglik_true <- ll_true$loglik
    rec$loglik_init <- ll_init$loglik
    rec$cache_error <- abs(ll_true$loglik - ll_combined$loglik)
    if (!all(is.finite(c(rec$loglik_true, rec$loglik_init, rec$cache_error)))) {
      stop("Non-finite likelihood or cache comparison.")
    }
    if (rec$cache_error > 1e-7 * (1 + abs(rec$loglik_true))) stop("Direct/combined likelihood mismatch.")
    tick <- proc.time()[["elapsed"]]
    fit <- do.call(fit_hawkesNet, c(list(params_init = params_init,
      time_window = c(0, horizon), mark_filtration = sim$net, PMF_mark = pmf,
      fixed_params = fixed, method = "L-BFGS-B", maxit = maxit, trace = 0,
      cores = 1L, cache_intensity = TRUE, combine_intensity = TRUE,
      verbose = FALSE, get_hessian = FALSE, run_sim = FALSE), extra))
    rec$fit_seconds <- proc.time()[["elapsed"]] - tick
    rec$convergence <- fit$fit$convergence
    bounds <- getFromNamespace("build_optim_bounds", "hawkesNet")(names(fit$fit$par))
    at_bound <- fit$fit$par <= bounds$lower + 2e-6 | fit$fit$par >= bounds$upper - 2e-6
    rec$boundary_parameters <- paste(names(fit$fit$par)[at_bound], collapse = ",")
    rec$loglik_fitted <- likelihood(fit$params)$loglik
    # fit$value is the maximized log likelihood (optim uses fnscale = -1).
    rec$optimizer_error <- abs(rec$loglik_fitted - fit$fit$value)
    if (!all(is.finite(c(fit$fit$par, rec$loglik_fitted, rec$optimizer_error)))) {
      stop("Non-finite fitted parameters or likelihood.")
    }
    if (rec$optimizer_error > 1e-7 * (1 + abs(rec$loglik_fitted))) {
      stop("Fitted cached objective differs from fresh direct likelihood.")
    }
    if (rec$loglik_fitted < rec$loglik_init - 1e-6) stop("Fit reduced the initial log likelihood.")
    free_names <- names(fit$fit$par)
    truth <- unlist(params_true)[free_names]
    stopifnot(all(is.finite(truth)))
    estimates[[length(estimates) + 1L]] <- data.frame(
      model = model, horizon = horizon, replicate = rep, parameter = free_names,
      truth = unname(truth), estimate = unname(fit$fit$par),
      error = unname(fit$fit$par - truth), convergence = rec$convergence)
    rec$status <- if (rec$convergence == 0L) "ok" else "optimizer_not_converged"
    fit$intens_funcs <- NULL
    list(simulation = sim, fit = fit, params_true = params_true,
         params_init = params_init, fixed = fixed)
  }, warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    rec$error <<- conditionMessage(e)
    NULL
  })
  rec$warnings <- paste(unique(warnings_seen), collapse = " | ")
  attempts[[length(attempts) + 1L]] <- as.data.frame(rec, stringsAsFactors = FALSE)
  saveRDS(list(record = rec, result = result), file.path(output_dir,
    sprintf("%s_T%s_rep%d.rds", model, format(horizon, scientific = FALSE), rep)))
  write.csv(do.call(rbind, attempts), file.path(output_dir, "attempts.csv"), row.names = FALSE)
  if (length(estimates)) write.csv(do.call(rbind, estimates),
    file.path(output_dir, "estimates.csv"), row.names = FALSE)
}
all_attempts <- do.call(rbind, attempts)
print(all_attempts[, c("model", "horizon", "replicate", "events", "status", "error")], row.names = FALSE)
if (length(estimates)) print(do.call(rbind, estimates), row.names = FALSE)
cat("Every attempted replicate was retained. This run does not test consistency.\n")
if (any(all_attempts$status != "ok")) quit(save = "no", status = 1L)
