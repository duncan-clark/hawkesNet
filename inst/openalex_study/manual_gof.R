# =============================================================================
# Interactive Model Fitting & GOF for OpenAlex results
# =============================================================================
# Run interactively from package root:
#   source("inst/openalex_study/manual_gof.R")
#
# This script:
#   1. Loads saved results (or fetches fresh data)
#   2. Lets you re-fit models interactively (fast — ~10-30s per fit)
#   3. Runs GOF simulations and compares observed vs simulated
#   4. Leaves all objects in the global env for exploration
#
# After sourcing, useful objects:
#   net_obs, inhom_bg, pfit, fit_result, sim_results, dat
#   Helper functions: refit(), quick_sim(), compare_nets()
# =============================================================================

library(hawkesNet)
library(network)
library(sna)
library(ernm)

# %||% available since R 4.0; define fallback for older R
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (!is.null(a)) a else b

cat("=== Interactive Hawkes-Net Model Explorer ===\n\n")

# =============================================================================
# 1. Load data
# =============================================================================
PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) Sys.getenv("SLURM_SUBMIT_DIR") else getwd()
rds_path <- file.path(PKG_ROOT, "cluster_output", "results_openalex_full.RDS")
if (!file.exists(rds_path)) rds_path <- "cluster_output/results_openalex_full.RDS"
if (!file.exists(rds_path)) {
  # Try parent
  rds_path <- file.path(dirname(PKG_ROOT), "cluster_output", "results_openalex_full.RDS")
}
stopifnot("Cannot find results_openalex_full.RDS" = file.exists(rds_path))
cat("Loading results from:", rds_path, "\n")
dat <- readRDS(rds_path)

# Extract key objects
net_obs    <- dat$net_raw
inhom_bg   <- dat$inhom_bg
TRUNCATION <- 300L
GROWTH_ONLY <- TRUE
time_window <- c(0, 1)

# Basic info
cat("\n=== Observed Network ===\n")
cat("  Nodes:", network.size(net_obs), "\n")
cat("  Edges:", network.edgecount(net_obs), "\n")
cat("  Events:", length(get_times(net_obs)$times), "\n")
cat("  Mean degree:", round(mean(degree(net_obs, gmode = "graph")), 2), "\n")
if ("gender" %in% list.vertex.attributes(net_obs)) {
  cat("  Gender:", paste(names(table(net_obs %v% "gender")), "=",
                         table(net_obs %v% "gender"), collapse = ", "), "\n")
}

# Available models from cluster run
cat("\n=== Available Fits from Cluster Run ===\n")
for (m in c("structural", "nodematch", "nodemix")) {
  fit_name <- paste0("fit_inhom_", m)
  f <- dat[[fit_name]]
  if (!is.null(f)) {
    cat(sprintf("  %-12s: converged=%d, iters=%d, loglik=%.2f\n",
                m, f$fit$convergence, f$fit$counts[1], f$fit$value))
  } else {
    cat(sprintf("  %-12s: (not available)\n", m))
  }
}

# =============================================================================
# 2. Config — edit these before re-sourcing or call refit() directly
# =============================================================================
N_CORES <- min(parallel::detectCores() - 1L, 8L)
MAX_ITER <- 10000L
N_GOF <- 25L
SEED_EVENTS <- 20L

# =============================================================================
# 3. Helper: reconstruct fitted params from a fit object
# =============================================================================
reconstruct_params <- function(fit, params_init, inhom_bg = NULL,
                               time_window = c(0, 1), verbose = TRUE) {
  par_vec <- fit$fit$par
  par_names <- names(par_vec)
  fixed <- fit$fixed_params

  pfit <- params_init
  pfit$vertex_categorical_levels <- params_init$vertex_categorical_levels

  # Scalars
  for (nm in c("mu", "beta_overall", "K", "beta_edges", "node_lambda", "m")) {
    if (nm %in% par_names) pfit[[nm]] <- par_vec[nm]
  }

  # CS_params
  cs_idx <- grep("^CS_params[0-9]+$", par_names)
  if (length(cs_idx) > 0) {
    cs_nums <- as.integer(sub("^CS_params", "", par_names[cs_idx]))
    for (j in seq_along(cs_idx)) {
      k <- cs_nums[j]
      if (k >= 1L && k <= length(pfit$CS_params)) pfit$CS_params[k] <- par_vec[cs_idx[j]]
    }
  }

  # vertex_categorical
  vc_idx <- grep("^vertex_categorical\\.", par_names)
  if (length(vc_idx) > 0 && !is.null(pfit$vertex_categorical)) {
    for (j in vc_idx) {
      parts <- strsplit(par_names[j], "\\.")[[1]]
      if (length(parts) >= 3) {
        attr_name <- parts[2]
        level_name <- paste(parts[3:length(parts)], collapse = ".")
        if (!is.null(pfit$vertex_categorical[[attr_name]])) {
          pfit$vertex_categorical[[attr_name]][level_name] <- par_vec[j]
        }
      }
    }
  }

  # Inhomogeneous background
  if (!is.null(inhom_bg) && !is.null(inhom_bg$mu_fit) && !is.null(inhom_bg$mu_fit$mu_fun)) {
    Tval <- time_window[2] - time_window[1]
    pfit$mu <- inhom_bg$integral_bg / Tval
  }

  # Clamp
  pfit$K <- min(max(pfit$K, 0.001), 0.999)
  pfit$mu <- max(pfit$mu, 0.001)
  pfit$node_lambda <- max(pfit$node_lambda, 0.1)
  pfit$beta_overall <- max(pfit$beta_overall, 0.001)
  pfit$beta_edges <- max(pfit$beta_edges, 0.001)

  # Repair vertex_categorical
  if (!is.null(pfit$vertex_categorical)) {
    pfit <- tryCatch(repair_vertex_categorical_params(pfit, eps = 1e-6),
                     error = function(e) pfit)
  }

  if (verbose) {
    cat("\n  Reconstructed parameters:\n")
    for (p in c("mu", "beta_overall", "K", "beta_edges", "node_lambda")) {
      if (!is.null(pfit[[p]])) {
        src <- if (!is.null(fixed) && p %in% fixed) "(FIXED)" else "(fitted)"
        cat(sprintf("    %-15s = %12.6f  %s\n", p, pfit[[p]], src))
      }
    }
    if (!is.null(pfit$CS_params)) {
      cat("    CS_params      =", paste(round(pfit$CS_params, 6), collapse = ", "), "\n")
    }
    if (!is.null(pfit$vertex_categorical)) {
      for (a in names(pfit$vertex_categorical)) {
        cat("    vertex_cat$", a, " =",
            paste(names(pfit$vertex_categorical[[a]]), "=",
                  round(pfit$vertex_categorical[[a]], 4), collapse = ", "), "\n")
      }
    }
  }
  pfit
}

# =============================================================================
# 4. Helper: print a fit summary
# =============================================================================
print_fit_summary <- function(fit_result, label = "") {
  if (is.null(fit_result)) { cat("  Fit is NULL\n"); return(invisible()) }
  if (nzchar(label)) cat("\n=== ", label, " ===\n") else cat("\n=== Fit Summary ===\n")
  cat("  Convergence:", fit_result$fit$convergence,
      if (fit_result$fit$convergence == 0) "(OK)" else "(WARNING: did not converge)", "\n")
  cat("  Function evaluations:", fit_result$fit$counts[1], "\n")
  cat("  Log-likelihood:", round(fit_result$fit$value, 4), "\n")
  if (!is.null(fit_result$fit_table)) {
    cat("\n")
    print(fit_result$fit_table, row.names = FALSE, right = FALSE)
    n_se <- sum(!is.na(fit_result$fit_table$std.error))
    if (n_se == 0) {
      cat("\n  *** All standard errors are NA — check Hessian diagnostics ***\n")
    } else if (n_se < nrow(fit_result$fit_table)) {
      cat(sprintf("\n  *** %d of %d SEs are NA ***\n",
                  nrow(fit_result$fit_table) - n_se, nrow(fit_result$fit_table)))
    }
  }
  if (!is.null(fit_result$hessian)) {
    eig <- eigen(-fit_result$hessian, symmetric = TRUE, only.values = TRUE)$values
    n_neg <- sum(eig < 0)
    cat(sprintf("\n  Hessian: %dx%d, eigenvalues of -H: %s\n",
                nrow(fit_result$hessian), ncol(fit_result$hessian),
                paste(sprintf("%.3g", eig), collapse = ", ")))
    if (n_neg > 0) cat("  *** WARNING:", n_neg, "negative eigenvalue(s) — MLE may not be a true max ***\n")
  }
  invisible(fit_result)
}

# =============================================================================
# 5. refit() — re-run a model fit interactively
# =============================================================================
#' Re-fit a model. Call as:
#'   fit_result <- refit("structural")
#'   fit_result <- refit("nodematch")
#'   fit_result <- refit("nodematch", maxit = 5000, trace = 1)
#'   fit_result <- refit("structural", params_init = my_custom_params)
#'
#' By default, initializes from the cluster-run results (warm start).
refit <- function(model = c("structural", "nodematch", "nodemix"),
                  params_init = NULL,
                  formula_RHS = NULL,
                  maxit = MAX_ITER,
                  trace = 0,
                  method = "Nelder-Mead",
                  reltol = 1e-8,
                  cores = N_CORES,
                  fixed_params = c("K", "mu"),
                  verbose = TRUE) {
  model <- match.arg(model)

  # --- Defaults by model type ---
  if (is.null(formula_RHS)) {
    formula_RHS <- switch(model,
      structural = dat$FORMULA_RHS_STRUCTURAL %||% "edges + triangles + gwdegree(0.5)",
      nodematch  = dat$FORMULA_RHS_NODEMATCH  %||% "edges + triangles + gwdegree(0.5) + nodeMatch('gender')",
      nodemix    = dat$FORMULA_RHS_NODEMIX    %||% "edges + triangles + gwdegree(0.5) + nodeMix('gender')"
    )
  }
  cat("\n=== Refitting:", model, "model ===\n")
  cat("  Formula:", formula_RHS, "\n")
  cat("  Method:", method, "| maxit:", maxit, "| cores:", cores, "\n")

  # --- Build params_init ---
  if (is.null(params_init)) {
    # Try warm-start from saved fit
    saved_fit <- dat[[paste0("fit_inhom_", model)]]
    saved_init <- dat[[paste0("params_init_", model)]]
    if (!is.null(saved_fit) && !is.null(saved_init)) {
      cat("  Warm-starting from saved cluster fit\n")
      params_init <- reconstruct_params(saved_fit, saved_init, inhom_bg,
                                         time_window, verbose = FALSE)
      # For optim, only pass the free params (exclude fixed)
      # But fit_hawkesNet handles this internally — just pass full params_init
      params_init <- saved_init
      # Overwrite with fitted values as starting point
      par_vec <- saved_fit$fit$par
      par_names <- names(par_vec)
      for (nm in c("beta_overall", "beta_edges", "node_lambda")) {
        if (nm %in% par_names) params_init[[nm]] <- par_vec[nm]
      }
      cs_idx <- grep("^CS_params[0-9]+$", par_names)
      if (length(cs_idx) > 0) {
        cs_nums <- as.integer(sub("^CS_params", "", par_names[cs_idx]))
        for (j in seq_along(cs_idx)) {
          k <- cs_nums[j]
          if (k >= 1L && k <= length(params_init$CS_params)) {
            params_init$CS_params[k] <- par_vec[cs_idx[j]]
          }
        }
      }
      vc_idx <- grep("^vertex_categorical\\.", par_names)
      if (length(vc_idx) > 0 && !is.null(params_init$vertex_categorical)) {
        for (j in vc_idx) {
          parts <- strsplit(par_names[j], "\\.")[[1]]
          if (length(parts) >= 3) {
            attr_name <- parts[2]
            level_name <- paste(parts[3:length(parts)], collapse = ".")
            if (!is.null(params_init$vertex_categorical[[attr_name]])) {
              params_init$vertex_categorical[[attr_name]][level_name] <- par_vec[j]
            }
          }
        }
      }
    } else {
      # Cold start
      cat("  Cold-starting (no saved fit for this model)\n")
      mu_init <- if (!is.null(inhom_bg)) inhom_bg$integral_bg else 100
      exp_cs <- expected_params_PMF_mark_CS(net_obs, formula_RHS)
      n_cs <- if (!is.na(exp_cs$CS_params_length)) exp_cs$CS_params_length else 3L
      params_init <- list(
        mu = mu_init,
        beta_overall = 1,
        K = 0.5,
        beta_edges = 1,
        node_lambda = 1,
        CS_params = c(-10, rep(0, n_cs - 1))
      )
      if (grepl("nodeMatch|nodeMix", formula_RHS)) {
        params_init$vertex_categorical <- list(gender = c(female = 0.1, male = 0.5))
        params_init$vertex_categorical_levels <- list(gender = c("female", "male", "unknown"))
      }
    }
  }

  cat("  Starting params:\n")
  for (nm in c("beta_overall", "beta_edges", "node_lambda")) {
    if (!is.null(params_init[[nm]])) cat(sprintf("    %-15s = %.6f\n", nm, params_init[[nm]]))
  }
  if (!is.null(params_init$CS_params)) {
    cat("    CS_params      =", paste(round(params_init$CS_params, 4), collapse = ", "), "\n")
  }

  # --- Build parscale ---
  n_cs <- length(params_init$CS_params)
  p_scale <- c(
    beta_overall = 0.1, beta_edges = 0.1, node_lambda = 1,
    setNames(rep(0.1, n_cs), paste0("CS_params", seq_len(n_cs)))
  )
  if (!is.null(params_init$vertex_categorical)) {
    for (a in names(params_init$vertex_categorical)) {
      for (lev in names(params_init$vertex_categorical[[a]])) {
        p_scale[paste0("vertex_categorical.", a, ".", lev)] <- 0.1
      }
    }
  }

  # --- Fit ---
  cat("\n  Fitting...\n")
  t0 <- proc.time()
  fit_result <- tryCatch(
    fit_hawkesNet(
      params_init = params_init,
      time_window = time_window,
      mark_filtration = net_obs,
      PMF_mark = PMF_mark_CS,
      mu_vec = if (!is.null(inhom_bg)) inhom_bg$mu_vec else NULL,
      integral_bg = if (!is.null(inhom_bg)) inhom_bg$integral_bg else NULL,
      formula_RHS = formula_RHS,
      truncation = TRUNCATION,
      mark_decay = "activity",
      growth_only = GROWTH_ONLY,
      max_node_time = 1,
      method = method,
      maxit = maxit,
      trace = trace,
      reltol = reltol,
      verbose = verbose,
      fixed_params = fixed_params,
      parscale = p_scale,
      cache_intensity = TRUE,
      combine_intensity = TRUE,
      cores = cores
    ),
    error = function(e) {
      cat("  *** FIT FAILED:", e$message, "***\n")
      NULL
    }
  )
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("\n  Fit completed in %.1f s\n", elapsed))

  if (!is.null(fit_result)) {
    print_fit_summary(fit_result, paste(model, "fit"))
    # Store globally for quick access
    assign("fit_result", fit_result, envir = .GlobalEnv)
    assign("pfit", reconstruct_params(fit_result, params_init, inhom_bg,
                                       time_window, verbose = TRUE),
           envir = .GlobalEnv)
    assign("formula_RHS_current", formula_RHS, envir = .GlobalEnv)
    assign("params_init_current", params_init, envir = .GlobalEnv)
  }
  invisible(fit_result)
}

# =============================================================================
# 6. quick_sim() — run a single simulation with current params
# =============================================================================
quick_sim <- function(pfit = NULL, growth_only = GROWTH_ONLY,
                      seed_events = SEED_EVENTS, verbose = TRUE) {
  if (is.null(pfit)) {
    pfit <- get0("pfit", envir = .GlobalEnv)
    if (is.null(pfit)) stop("No pfit available. Run refit() first or pass pfit.")
  }
  formula_RHS <- get0("formula_RHS_current", envir = .GlobalEnv)
  if (is.null(formula_RHS)) formula_RHS <- "edges + triangles + gwdegree(0.5)"

  # Seed
  seed_net <- NULL; seed_times <- NULL
  if (seed_events > 0) {
    all_times <- get_times(net_obs)$times
    if (length(all_times) >= seed_events) {
      t_seed <- all_times[seed_events]
      seed_net <- filtration_to_net(net_obs, t_seed, equals = TRUE)
      seed_times <- all_times[1:seed_events]
      if (verbose) cat(sprintf("  Seeding with first %d events (t <= %.4f), %d nodes, %d edges\n",
                               seed_events, t_seed, network.size(seed_net), network.edgecount(seed_net)))
    }
  }

  cat("  Simulating (growth_only =", growth_only, ")...\n")
  t0 <- proc.time()
  sim <- tryCatch(
    sim_hawkesNet(
      params = pfit,
      time_window = time_window,
      PMF_mark = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS = formula_RHS,
      truncation = TRUNCATION,
      mark_decay = "activity",
      growth_only = growth_only,
      max_node_time = 1,
      hashed_edges = TRUE,
      verbose = verbose,
      mu_multiplier = 5,
      stop_on_full_network = FALSE,
      inhom_bg = inhom_bg,
      seed_net = seed_net,
      seed_times = seed_times
    ),
    error = function(e) {
      cat("  *** SIMULATION FAILED:", e$message, "***\n")
      NULL
    }
  )
  elapsed <- (proc.time() - t0)[3]
  if (!is.null(sim) && !is.null(sim$net)) {
    cat(sprintf("  Done in %.1f s: %d nodes, %d edges, %d events, mean degree %.2f\n",
                elapsed, network.size(sim$net), network.edgecount(sim$net),
                length(sim$events$t), mean(degree(sim$net, gmode = "graph"))))
  }
  invisible(sim)
}

# =============================================================================
# 7. compare_nets() — compare observed vs one or more simulated nets
# =============================================================================
compare_nets <- function(..., labels = NULL) {
  sims <- list(...)
  if (is.null(labels)) labels <- paste0("Sim_", seq_along(sims))

  cat(sprintf("\n  %-20s %10s", "", "Observed"))
  for (l in labels) cat(sprintf(" %15s", l))
  cat("\n")

  row <- function(stat, obs_val, sim_vals) {
    cat(sprintf("  %-20s %10s", stat, if (is.na(obs_val)) "NA" else sprintf("%d", obs_val)))
    for (v in sim_vals) cat(sprintf(" %15s", if (is.na(v)) "NA" else sprintf("%d", v)))
    cat("\n")
  }
  rowf <- function(stat, obs_val, sim_vals) {
    cat(sprintf("  %-20s %10.2f", stat, obs_val))
    for (v in sim_vals) cat(sprintf(" %15.2f", v))
    cat("\n")
  }

  n_obs <- network.size(net_obs)
  e_obs <- network.edgecount(net_obs)
  d_obs <- mean(degree(net_obs, gmode = "graph"))
  ev_obs <- length(get_times(net_obs)$times)

  ns <- sapply(sims, function(s) if (!is.null(s$net)) network.size(s$net) else NA)
  es <- sapply(sims, function(s) if (!is.null(s$net)) network.edgecount(s$net) else NA)
  ds <- sapply(sims, function(s) if (!is.null(s$net)) mean(degree(s$net, gmode = "graph")) else NA)
  evs <- sapply(sims, function(s) if (!is.null(s$events)) length(s$events$t) else NA)

  row("Nodes", n_obs, ns)
  row("Edges", e_obs, es)
  rowf("Mean degree", d_obs, ds)
  row("Events", ev_obs, evs)
}

# =============================================================================
# 8. run_gof() — full GOF with parallel simulations
# =============================================================================
run_gof <- function(model = "structural", n_sim = N_GOF,
                    seed_events = SEED_EVENTS, cores = N_CORES) {
  fit <- dat[[paste0("fit_inhom_", model)]]
  params_init <- dat[[paste0("params_init_", model)]]
  formula_RHS <- switch(model,
    structural = dat$FORMULA_RHS_STRUCTURAL %||% "edges + triangles + gwdegree(0.5)",
    nodematch  = dat$FORMULA_RHS_NODEMATCH  %||% "edges + triangles + gwdegree(0.5) + nodeMatch('gender')",
    nodemix    = dat$FORMULA_RHS_NODEMIX    %||% "edges + triangles + gwdegree(0.5) + nodeMix('gender')"
  )

  # If we just re-fit, use the fresh results
  if (exists("fit_result", envir = .GlobalEnv) && exists("formula_RHS_current", envir = .GlobalEnv)) {
    if (get("formula_RHS_current", envir = .GlobalEnv) == formula_RHS) {
      cat("  Using fresh fit from refit() instead of saved cluster fit\n")
      fit <- get("fit_result", envir = .GlobalEnv)
      params_init <- get("params_init_current", envir = .GlobalEnv)
    }
  }

  if (is.null(fit)) { cat("  No fit available for model:", model, "\n"); return(invisible(NULL)) }

  cat(sprintf("\n=== GOF for %s model (%d sims, %d cores) ===\n", model, n_sim, cores))
  t0 <- proc.time()
  gof_result <- tryCatch(
    gof(
      fit = fit,
      net_obs = net_obs,
      params_init = params_init,
      PMF_mark = PMF_mark_CS,
      cond_intensity = cond_intensity,
      formula_RHS = formula_RHS,
      time_window = time_window,
      truncation = TRUNCATION,
      mark_decay = "activity",
      growth_only = GROWTH_ONLY,
      max_node_time = 1,
      inhom_bg = inhom_bg,
      n_sim = n_sim,
      cores = cores,
      seed_events = seed_events,
      verbose = TRUE
    ),
    error = function(e) {
      cat("  *** GOF FAILED:", e$message, "***\n")
      NULL
    }
  )
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  GOF completed in %.1f s\n", elapsed))
  if (!is.null(gof_result)) {
    assign("gof_result", gof_result, envir = .GlobalEnv)
    # Summary stats
    if (!is.null(gof_result$degree_obs) && !is.null(gof_result$degree_sim)) {
      cat("  Degree obs mean:", round(mean(gof_result$degree_obs), 2), "\n")
      cat("  Degree sim mean:", round(mean(rowMeans(gof_result$degree_sim, na.rm = TRUE), na.rm = TRUE), 2), "\n")
    }
    if (!is.null(gof_result$nets_sim)) {
      n_ok <- sum(!sapply(gof_result$nets_sim, is.null))
      cat(sprintf("  Successful simulations: %d / %d\n", n_ok, n_sim))
    }
  }
  invisible(gof_result)
}

# =============================================================================
# 9. Load the default model and print summary
# =============================================================================
MODEL <- "structural"   # Change to "nodematch" or "nodemix" as desired
fit_default <- dat[[paste0("fit_inhom_", MODEL)]]
params_init_default <- dat[[paste0("params_init_", MODEL)]]

if (!is.null(fit_default) && !is.null(params_init_default)) {
  formula_RHS_current <- switch(MODEL,
    structural = dat$FORMULA_RHS_STRUCTURAL %||% "edges + triangles + gwdegree(0.5)",
    nodematch  = dat$FORMULA_RHS_NODEMATCH  %||% "edges + triangles + gwdegree(0.5) + nodeMatch('gender')",
    nodemix    = dat$FORMULA_RHS_NODEMIX    %||% "edges + triangles + gwdegree(0.5) + nodeMix('gender')"
  )
  params_init_current <- params_init_default
  print_fit_summary(fit_default, paste("Saved", MODEL, "fit"))
  pfit <- reconstruct_params(fit_default, params_init_default, inhom_bg, time_window)
} else {
  cat("\n  No saved fit for", MODEL, "model.\n")
  cat("  Run: fit_result <- refit(\"", MODEL, "\")\n", sep = "")
}

# =============================================================================
# 10. Usage guide
# =============================================================================
cat("\n")
cat("=============================================================\n")
cat("  Interactive Commands\n")
cat("=============================================================\n")
cat("\n")
cat("  Re-fit models (fast: ~10-30s):\n")
cat("    fit_result <- refit('structural')\n")
cat("    fit_result <- refit('nodematch')\n")
cat("    fit_result <- refit('nodematch', trace = 1)        # see NM iterations\n")
cat("    fit_result <- refit('structural', maxit = 20000)   # more iterations\n")
cat("\n")
cat("  Re-fit with custom starting params:\n")
cat("    my_params <- list(mu = 300, beta_overall = 5, K = 0.5,\n")
cat("                      beta_edges = 0.5, node_lambda = 1.2,\n")
cat("                      CS_params = c(-2, 0, -2))\n")
cat("    fit_result <- refit('structural', params_init = my_params)\n")
cat("\n")
cat("  Quick simulation with fitted params:\n")
cat("    sim1 <- quick_sim()                     # uses pfit from refit()\n")
cat("    sim2 <- quick_sim(growth_only = FALSE)  # all-edges mode\n")
cat("    compare_nets(sim1, sim2, labels = c('Growth', 'AllEdges'))\n")
cat("\n")
cat("  Full GOF:\n")
cat("    gof_result <- run_gof('structural', n_sim = 25)\n")
cat("    gof_result <- run_gof('nodematch', n_sim = 50, cores = 4)\n")
cat("\n")
cat("  Inspect objects:\n")
cat("    pfit                    # current fitted parameters\n")
cat("    fit_result$fit_table    # estimates and SEs\n")
cat("    fit_result$hessian      # numerical Hessian\n")
cat("    net_obs                 # observed network\n")
cat("    dat                     # full saved results from cluster\n")
cat("=============================================================\n")
