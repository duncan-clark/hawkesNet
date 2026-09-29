# Three explicit mark laws: BA-1, CS-1 and CS-2.
# Run from the package root after pkgload::load_all(".") or library(hawkesNet).
# All three use the same normalized-mark simulator and cached likelihood fitter.

run_mark_construction_example <- function(model = c("BA-1", "CS-1", "CS-2"),
                                          horizon = 5, seed = 20260922L) {
  model <- match.arg(model)
  params <- list(mu=4, K=.25, beta_overall=1.2, beta_edges=.5, m=1.4)
  if (model == "BA-1") {
    pmf <- PMF_mark_BA
    options <- list(truncation=12L)
  } else {
    params$node_lambda <- .8
    params$CS_params <- c(0,.6)
    pmf <- PMF_mark_CS
    options <- list(formula_RHS="edges + triangles", truncation=4L,
      cs_mode=if (model == "CS-1") "independent" else "size_conditional")
  }
  set.seed(seed)
  sim <- do.call(sim_hawkesNet,c(list(params=params,time_window=c(0,horizon),
    PMF_mark=pmf,cond_intensity=cond_intensity,verbose=FALSE),options))
  initial <- params; initial$m <- .9
  # A one-parameter smoke fit: this demonstrates the fit interface, not
  # consistency. CS-2's edge-count coefficient cancels and must be fixed.
  fit <- do.call(fit_hawkesNet,c(list(params_init=initial,time_window=c(0,horizon),
    mark_filtration=sim$net,PMF_mark=pmf,fixed_params=setdiff(names(params),"m"),
    method="L-BFGS-B",maxit=100L,get_hessian=FALSE,verbose=FALSE),options))
  list(model=model,truth=params,options=options,simulation=sim,fit=fit)
}

# Example:
# result <- run_mark_construction_example("CS-2")
# result$fit$fit_table
# CS-1 uses logistic single-edge change-statistic weights followed by Poisson
# collapsing; m is the attempt mean. CS-2 draws K from Poisson(m) conditioned
# on K <= D, then weights complete K-edge updates by their whole change.
# Exact CS-2 work grows as choose(D,K); keep its candidate window modest.
