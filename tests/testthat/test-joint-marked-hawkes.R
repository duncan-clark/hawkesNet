test_that("joint marked state decays and accumulates valued history", {
    state <- JointMarkedState(actors = 1:3, beta_edges = log(2))
    state <- .jmw_decay_state(state, 0)
    state <- .jmw_add_event(state, 1, 2, weight = 1, amount = expm1(1))
    expect_equal(state$edge_w[1, 2], 1)
    expect_equal(state$edge_v[1, 2], 1)
    state <- .jmw_decay_state(state, 1)
    expect_equal(state$edge_w[1, 2], 0.5, tolerance = 1e-10)
})

test_that("joint mark probabilities are non-uniform under nonzero RHEM params", {
    params <- joint_marked_hawkes_default_params()
    params$RHEM_params[] <- 0
    params$RHEM_params["repetition"] <- 2
    params$RHEM_params["edge_value"] <- 1
    hits <- data.frame(
        t = c(0.1, 0.2, 0.3),
        i = c(1L, 1L, 2L),
        j = c(2L, 2L, 1L),
        weight = c(1, 1, 0.5),
        amount = c(10, 10, 2),
        stringsAsFactors = FALSE
    )
    actors <- 1:4
    path <- .jmw_state_path(hits, actors, params$beta_edges)
    state <- path[[3]]
    all <- .jmw_all_dyad_stats(state)
    eta <- as.vector(all$stats %*% as.numeric(params$RHEM_params[colnames(all$stats)]))
    probs <- exp(eta - max(eta))
    probs <- probs / sum(probs)
    expect_gt(max(probs) - min(probs), 0.05)
    idx_rep <- which(all$risk$i == 1 & all$risk$j == 2)
    expect_gt(probs[idx_rep], 1 / nrow(all$risk))
})

test_that("mark-dependent excitation changes ground intensity", {
    params <- joint_marked_hawkes_default_params(n_events = 20, duration = 1)
    params$gamma_value <- 1
    params$K <- 1
    params$mu <- 0.1
    small <- data.frame(
        t = c(0.1, 0.2),
        i = c(1L, 1L),
        j = c(2L, 2L),
        weight = c(0.2, 0.2),
        amount = c(1, 100),
        stringsAsFactors = FALSE
    )
    actors <- 1:3
    # Score second event under high previous amount vs low
    low <- small
    low$amount[1] <- 1
    high <- small
    high$amount[1] <- 1000
    s_low <- score_joint_marked_hawkes(low, params = params, actors = actors, score_idx = 2)
    s_high <- score_joint_marked_hawkes(high, params = params, actors = actors, score_idx = 2)
    expect_gt(s_high$log_ground[1], s_low$log_ground[1])
})

test_that("joint loglik is finite on simulated data", {
    set.seed(42)
    params <- joint_marked_hawkes_default_params()
    params$mu <- 25
    params$K <- 0.3
    sim <- sim_joint_marked_hawkes(params, time_window = c(0, 1), actors = 1:5, max_events = 60)
    expect_gte(nrow(sim$hits), 5)
    ll <- loglik_joint_marked_hawkes(params, sim$hits, sim$actors, c(0, 1))
    expect_true(is.finite(ll$loglik))
    expect_true(is.finite(ll$mark_loglik))
    expect_true(is.finite(ll$amount_loglik))
    expect_true(is.finite(ll$ground_loglik))
})

test_that("joint fit recovers positive value coupling on synthetic signal", {
    actors <- 1:4
    # Pure repeated dyad: mark MLE for repetition must be strongly positive
    hits <- data.frame(
        t = seq(0.02, 0.98, length.out = 25),
        i = 1L,
        j = 2L,
        amount = 50,
        weight = log1p(50),
        stringsAsFactors = FALSE
    )
    params <- joint_marked_hawkes_default_params(n_events = nrow(hits), duration = 1)
    params$RHEM_params[] <- 0
    params$gamma_value <- 0
    params$gamma_edge <- 0
    params$K <- 0
    params$beta_edges <- 0
    params$amount_edge_value <- 0
    params$amount_repetition <- 0

    obj <- function(rep) {
        p <- params
        p$RHEM_params["repetition"] <- rep
        loglik_joint_marked_hawkes(p, hits, actors, c(0, 1))$mark_loglik
    }
    fit1 <- stats::optimize(obj, interval = c(-1, 5), maximum = TRUE)
    expect_gt(fit1$maximum, 1)

    # Larger past amounts raise future ground intensity when gamma_value > 0
    params$gamma_value <- 1
    params$K <- 1
    low <- hits
    low$amount[1:5] <- 1
    low$weight[1:5] <- log1p(1)
    high <- hits
    high$amount[1:5] <- 1000
    high$weight[1:5] <- log1p(1000)
    s_low <- score_joint_marked_hawkes(low, params = params, actors = actors, score_idx = 6)
    s_high <- score_joint_marked_hawkes(high, params = params, actors = actors, score_idx = 6)
    expect_gt(s_high$log_ground[1], s_low$log_ground[1])

    # Amount model tracks edge_value history when amounts increase with repeats
    hits_amt <- hits
    hits_amt$amount <- 2 * seq_len(nrow(hits))
    hits_amt$weight <- log1p(hits_amt$amount)
    p_amt <- params
    p_amt$gamma_value <- 0
    p_amt$K <- 0
    p_amt$amount_intercept <- log1p(2)
    p_amt$amount_sd <- 0.5
    obj_a <- function(b) {
        p <- p_amt
        p$amount_edge_value <- b
        loglik_joint_marked_hawkes(p, hits_amt, actors, c(0, 1))$amount_loglik
    }
    fit_a <- stats::optimize(obj_a, interval = c(-0.5, 2), maximum = TRUE)
    expect_gt(fit_a$maximum, 0)
})

test_that("benchmark_joint_marked_hawkes returns scaling table", {
    tab <- benchmark_joint_marked_hawkes(
        actor_grid = c(5, 8),
        event_grid = c(10),
        reps = 1
    )
    expect_true(all(c("n_actors", "n_events", "risk_dyads", "seconds") %in% names(tab)))
    expect_equal(nrow(tab), 2)
    expect_true(all(is.finite(tab$seconds)))
})
