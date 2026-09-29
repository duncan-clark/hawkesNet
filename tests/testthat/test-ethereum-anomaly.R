study_ethereum_root <- function() {
    # Studies live outside the package: <project>/studies/ethereum
    candidates <- c(
        file.path(testthat::test_path("..", "..", ".."), "studies", "ethereum"),
        file.path(testthat::test_path("..", ".."), "..", "studies", "ethereum")
    )
    hit <- candidates[dir.exists(candidates)]
    if (length(hit) == 0) {
        return(NA_character_)
    }
    normalizePath(hit[1], mustWork = FALSE)
}

source_study <- function() {
    root <- study_ethereum_root()
    skip_if(is.na(root), "Sibling studies/ethereum not found")
    env <- new.env(parent = globalenv())
    sys.source(file.path(root, "ethereum_utils.R"), envir = env)
    sys.source(file.path(root, "ethereum_anomaly.R"), envir = env)
    env
}

make_toy_transfers <- function(n = 200, n_actors = 8, seed = 1) {
    set.seed(seed)
    addrs <- sprintf("0x%040x", seq_len(n_actors))
    t0 <- as.numeric(as.POSIXct("2024-03-01", tz = "UTC"))
    from <- sample(addrs, n, replace = TRUE)
    to <- sample(addrs, n, replace = TRUE)
    same <- from == to
    to[same] <- addrs[((match(from[same], addrs)) %% n_actors) + 1]
    data.frame(
        block_timestamp = as.POSIXct(t0 + seq_len(n) * 30, origin = "1970-01-01", tz = "UTC"),
        transaction_hash = paste0("0x", sprintf("%04x", seq_len(n))),
        log_index = seq_len(n),
        block_number = 1000L + seq_len(n),
        token_address = "0xa0b86991c6218b36c1d19d4a2e9eb0ce3606eb48",
        from_address = from,
        to_address = to,
        raw_value = as.character(as.integer(10^runif(n, 4, 8))),
        stringsAsFactors = FALSE
    )
}

test_that("leakage-safe prep uses train-only actors and calendar split", {
    env <- source_study()
    raw <- make_toy_transfers()
    transfers <- env$eth_standardize_transfers(raw)
    prepared <- env$eth_prepare_anomaly_data(
        transfers,
        labels = NULL,
        max_entities = 5,
        max_events = 120,
        train_frac = 0.6,
        valid_frac = 0.1
    )
    expect_true(prepared$metadata$leakage_safe)
    expect_lte(prepared$metadata$n_actors, 5)
    expect_true(all(c("train", "test") %in% unique(prepared$hits$split)))
    train_actors <- unique(c(prepared$hits$i[prepared$hits$split == "train"],
                             prepared$hits$j[prepared$hits$split == "train"]))
    all_actors <- unique(c(prepared$hits$i, prepared$hits$j))
    expect_true(all(all_actors %in% train_actors))
})

test_that("injection evaluation ranks amount anomalies", {
    env <- source_study()
    set.seed(1)
    hits <- data.frame(
        t = seq(0, 1, length.out = 40),
        i = rep(1:4, length.out = 40),
        j = rep(c(2, 3, 4, 1), length.out = 40),
        weight = runif(40),
        amount = rexp(40, 0.01),
        token_symbol = "USDC",
        event_id = paste0("e", 1:40),
        from_entity = "a",
        to_entity = "b",
        split = c(rep("train", 28), rep("test", 12)),
        stringsAsFactors = FALSE
    )
    inj <- env$eth_inject_anomalies(hits, n_inject = 5, kinds = "amount", seed = 2)
    expect_true(any(inj$hits$is_injected))
    expect_true(all(inj$hits$is_injected[inj$injected$event_index]))
    expect_equal(
        inj$hits$inject_kind[inj$injected$event_index],
        inj$injected$inject_kind
    )
    sc <- env$eth_score_amount_baseline(inj$hits)
    eval <- env$eth_evaluate_injections(sc, inj$hits)
    expect_equal(nrow(eval), 1)
    expect_true(is.finite(eval$auc_total) || is.finite(eval$auc_amount))
    auc <- if (is.finite(eval$auc_total[1])) eval$auc_total[1] else eval$auc_amount[1]
    expect_gt(auc, 0.7)
})

test_that("case study table and weak-label helpers run", {
    env <- source_study()
    hits <- data.frame(
        t = 1:10 / 10,
        i = c(1, 1, 2, 2, 3, 1, 2, 3, 1, 2),
        j = c(2, 3, 1, 3, 1, 3, 1, 2, 2, 3),
        amount = 10 * (1:10),
        token_symbol = "USDC",
        split = c(rep("train", 7), rep("test", 3)),
        stringsAsFactors = FALSE
    )
    scores <- data.frame(
        event_index = 8:10,
        model = "simple_rhem",
        score_total = c(3, 1, 2),
        score_mark = c(3, 1, 2),
        score_ground = c(0.1, 0.1, 0.1),
        score_amount = c(0.2, 0.2, 0.5),
        stringsAsFactors = FALSE
    )
    actor_map <- data.frame(
        actor_id = 1:3,
        entity = c("addr:0x1", "addr:0x2", "addr:0x3"),
        stringsAsFactors = FALSE
    )
    actor_map <- env$eth_merge_weak_labels(
        actor_map,
        labels_df = data.frame(
            address = "0x1",
            category = "phish-hack",
            stringsAsFactors = FALSE
        )
    )
    expect_true(any(actor_map$is_risk))
    cases <- env$eth_case_study_table(scores, hits, actor_map, top_n = 2)
    expect_equal(nrow(cases), 2)
    expect_equal(cases$rank[1], 1)
})
