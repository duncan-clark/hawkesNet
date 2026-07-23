source(test_path("..", "..", "inst", "ethereum_study", "ethereum_utils.R"))

ethereum_test_prepared <- function() {
    raw <- data.frame(
        block_timestamp = c(
            "2024-03-01 00:00:00",
            "2024-03-01 00:00:00",
            "2024-03-01 00:00:03",
            "2024-03-01 00:00:10",
            "2024-03-01 00:00:11",
            "2024-03-01 00:00:12"
        ),
        transaction_hash = paste0("0xtx", 1:6),
        log_index = 0:5,
        token_address = c(
            "0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48",
            "0xdAC17F958D2ee523a2206206994597C13D831ec7",
            "0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48",
            "0x6B175474E89094C44Da98b954EedeAC495271d0F",
            "0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48",
            "0xdAC17F958D2ee523a2206206994597C13D831ec7"
        ),
        from_address = c("0xA", "0xB", "0xA", "0xC", "0xB", "0xC"),
        to_address = c("0xB", "0xA", "0xC", "0xA", "0xC", "0xB"),
        raw_value = c(1000000, 2500000, 4000000, 1e18, 3000000, 1500000),
        stringsAsFactors = FALSE
    )
    transfers <- eth_standardize_transfers(raw) # nolint: object_usage_linter
    labels <- data.frame(
        address = eth_normalize_address(c("0xA", "0xB", "0xC")), # nolint: object_usage_linter
        entity = c("Exchange A", "Treasury B", "Protocol C"),
        category = c("exchange", "treasury", "protocol"),
        stringsAsFactors = FALSE
    )
    eth_prepare_hawkesnet_data( # nolint: object_usage_linter
        transfers,
        labels = labels,
        max_entities = 3,
        max_events = NULL
    )
}

test_that("Ethereum stablecoin transfers standardize to strict repeated hits", {
    raw <- data.frame(
        block_timestamp = c(
            "2024-03-01 00:00:00",
            "2024-03-01 00:00:00",
            "2024-03-01 00:00:03",
            "2024-03-01 00:00:10"
        ),
        transaction_hash = paste0("0xtx", 1:4),
        log_index = 0:3,
        token_address = c(
            "0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48",
            "0xdAC17F958D2ee523a2206206994597C13D831ec7",
            "0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48",
            "0x6B175474E89094C44Da98b954EedeAC495271d0F"
        ),
        from_address = c("0xA", "0xB", "0xA", "0xC"),
        to_address = c("0xB", "0xA", "0xC", "0xA"),
        raw_value = c(1000000, 2500000, 4000000, 1e18),
        stringsAsFactors = FALSE
    )

    transfers <- eth_standardize_transfers(raw)
    expect_equal(transfers$token_symbol, c("USDC", "USDT", "USDC", "DAI"))
    expect_equal(transfers$amount, c(1, 2.5, 4, 1))
    expect_true(all(diff(transfers$t) > 0))

    labels <- data.frame(
        address = eth_normalize_address(c("0xA", "0xB", "0xC")),
        entity = c("Exchange A", "Treasury B", "Protocol C"),
        category = c("exchange", "treasury", "protocol"),
        stringsAsFactors = FALSE
    )
    prepared <- eth_prepare_hawkesnet_data(
        transfers,
        labels = labels,
        max_entities = 3,
        max_events = NULL
    )

    expect_equal(nrow(prepared$hits), 4)
    expect_true(all(c("t", "i", "j", "weight") %in% names(prepared$hits)))
    expect_equal(nrow(prepared$actor_map), 3)
    expect_true(all(c("out_count", "in_count", "net_amount", "category") %in% names(prepared$node_covariates)))
    expect_equal(sum(prepared$hits$split == "train"), 2)
})

test_that("Ethereum comparison baselines return ranked held-out scores", {
    prepared <- ethereum_test_prepared()
    cmp <- eth_run_baseline_comparisons(prepared$hits, history_maxit = 2)

    expect_true(all(c("dyad_frequency", "dyad_recency", "sender_receiver_activity", "history_softmax") %in% cmp$metrics$model))
    expect_true(all(is.finite(cmp$scores$log_mark_density)))
    expect_true(all(c("mean_negative_log_mark", "top10", "mrr") %in% names(cmp$metrics)))
})

test_that("No-excitation RHEM is fit as HawkesNet with K fixed at zero", {
    skip_if_not_installed("ernm")
    prepared <- ethereum_test_prepared()
    fit <- eth_fit_no_excitation_rhem(
        prepared$hits,
        formula_RHS = NULL,
        maxit = 2,
        fixed_params = c("K", "beta_overall", "beta_edges")
    )
    scores <- eth_score_hawkes_marks(
        prepared$hits,
        fit,
        formula_RHS = NULL,
        model = "rhem_no_excitation_hawkesnet"
    )

    expect_equal(fit$params_init$K, 0)
    expect_true("K" %in% fit$fixed_params)
    expect_equal(unique(scores$model), "rhem_no_excitation_hawkesnet")
    expect_true(all(is.finite(scores$log_mark_density)))
})

test_that("Ethereum timeNet temporal formulas initialize and score", {
    skip_if_not_installed("timeNet")
    prepared <- ethereum_test_prepared()
    formula_RHS <- eth_timenet_formula(beta = 2, terms = c("edge", "recip", "sender"))
    params <- eth_initial_params_for_formula(prepared$hits, formula_RHS)

    expect_equal(names(params$RHEM_params),
                 c("decayedEdgeValue", "decayedRecipValue", "decayedSenderActivity"))

    fit <- eth_fit_no_excitation_rhem(
        prepared$hits,
        formula_RHS = formula_RHS,
        maxit = 1,
        fixed_params = c("K", "beta_overall", "beta_edges")
    )
    scores <- eth_score_hawkes_marks(
        prepared$hits,
        fit,
        formula_RHS = formula_RHS,
        model = "rhem_timenet_smoke"
    )

    expect_equal(unique(scores$model), "rhem_timenet_smoke")
    expect_true(all(is.finite(scores$log_mark_density)))
})

test_that("Package-level RHEM likelihood scores repeated-hit marks", {
    prepared <- ethereum_test_prepared()
    hits <- prepared$hits[, c("t", "i", "j", "weight"), drop = FALSE]
    params <- list(beta_edges = 0, RHEM_params = eth_simple_rhem_params())

    ll <- hawkesNet::loglik_RHEM(
        params = params,
        mark_filtration = hits,
        formula_RHS = NULL
    )
    scores <- hawkesNet::score_RHEM(
        mark_filtration = hits,
        params = params,
        formula_RHS = NULL,
        score_idx = seq_len(nrow(hits))
    )

    expect_true(is.finite(ll$loglik))
    expect_equal(nrow(scores), nrow(hits))
    expect_true(all(is.finite(scores$negative_log_mark)))
})

test_that("Ethereum BigQuery SQL targets stablecoin token transfers", {
    sql <- eth_build_stablecoin_sql("2024-03-01", "2024-03-02", max_rows = 10)
    expect_match(sql, "bigquery-public-data.crypto_ethereum.token_transfers", fixed = TRUE)
    expect_match(sql, "0xa0b86991c6218b36c1d19d4a2e9eb0ce3606eb48", fixed = TRUE)
    expect_match(sql, "0xdac17f958d2ee523a2206206994597c13d831ec7", fixed = TRUE)
    expect_match(sql, "LIMIT 10", fixed = TRUE)
})

test_that("Ethereum CSV reader preserves hexadecimal addresses", {
    path <- tempfile(fileext = ".csv")
    raw <- data.frame(
        block_timestamp = "2024-03-01 00:00:00",
        transaction_hash = "0xabc123",
        log_index = 0,
        token_address = "0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48",
        from_address = "0x00000000000000000000000000000000000000aa",
        to_address = "0x00000000000000000000000000000000000000bb",
        raw_value = "1000000",
        stringsAsFactors = FALSE
    )
    utils::write.csv(raw, path, row.names = FALSE)

    transfers <- eth_read_transfers(path)
    expect_equal(transfers$from_address, "0x00000000000000000000000000000000000000aa")
    expect_equal(transfers$to_address, "0x00000000000000000000000000000000000000bb")
    expect_equal(transfers$token_symbol, "USDC")
})

test_that("Ethereum label reader supports Etherscan JSON snapshots", {
    skip_if_not_installed("jsonlite")
    path <- tempfile(fileext = ".json")
    writeLines(
        '{"0xabc":{"name":"Exchange Hot Wallet","labels":["exchange","hot-wallet"]}}',
        path
    )

    labels <- eth_read_labels(path)
    expect_equal(labels$address, "0xabc")
    expect_equal(labels$entity, "Exchange Hot Wallet")
    expect_equal(labels$category, "exchange;hot-wallet")
})
