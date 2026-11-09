# Small complete-update workflow; no external data or output files.
# Run with source(system.file("examples", "example_network_updates.R",
#                              package = "hawkesNet"))
library(hawkesNet)

model <- network_update_model(mark = "PA", timing = "M2",
                              summary = "anchor_threshold")
set.seed(12)
path <- simulate_network_updates(model, horizon = 5)
print(loglik_network_updates(path, model)[c("loglik", "timing", "marks")])

# Only mu is estimated to keep this an interface example, not a recovery study.
fit <- fit_network_updates(path, model,
                           fixed = setdiff(names(model$parameters), "mu"))
print(fit$parameters)
print(c(events = length(path$times), convergence = fit$convergence))

# Single edits and a coordinated migration can produce the same final update.
# The probability marginalizes those latent proposal combinations.
migration <- poisson_proposals(
  groups = list("delete_old", "add_new", c("delete_old", "add_new")),
  rates = c(0.2, 0.05, 0.15)
)
print(exp(logprob_update(migration, c("delete_old", "add_new"))))
print(sample_update(migration))
