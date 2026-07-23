source(test_path("..", "..", "inst", "amlsim_study", "amlsim_utils.R"))

test_that("AMLSim sample-style transactions become repeated hits", {
    tx <- data.frame(
        TXN_ID = 1:5,
        ACCOUNT_ID = c("A", "A", "B", "C", "D"),
        COUNTER_PARTY_ACCOUNT_NUM = c("B", "C", "A", "A", "D"),
        TXN_SOURCE_TYPE_CODE = "WIRE",
        TXN_AMOUNT_ORIG = c(100, 200, 150, 75, 10),
        start = c(1, 1, 2, 3, 4),
        stringsAsFactors = FALSE
    )

    standardized <- amlsim_standardize_transactions(tx)
    expect_equal(nrow(standardized), 4)
    expect_true(all(diff(standardized$t) > 0))
    expect_false(any(standardized$from_id == standardized$to_id))

    hits <- amlsim_make_repeated_hits(standardized)
    expect_true(all(c("t", "i", "j", "weight", "is_sar") %in% names(hits)))
    expect_equal(nrow(attr(hits, "actor_map")), 3)
    expect_true(all(hits$weight > 0))
})

test_that("AMLSim simulator-log transactions preserve SAR labels", {
    tx <- data.frame(
        step = c(0, 0, 1),
        type = "TRANSFER",
        amount = c(100, 200, 300),
        nameOrig = c("1", "2", "3"),
        nameDest = c("2", "3", "1"),
        isSAR = c("0", "1", "0"),
        alertID = c("-1", "7", "-1"),
        stringsAsFactors = FALSE
    )

    standardized <- amlsim_standardize_transactions(tx)
    expect_equal(standardized$is_sar, c(FALSE, TRUE, FALSE))

    hits <- amlsim_make_repeated_hits(standardized, weight = "count")
    hits$split <- amlsim_temporal_split(hits, train_frac = 2 / 3)
    features <- amlsim_history_features(hits)

    expect_equal(nrow(features), 3)
    expect_true(all(c("repetition", "reciprocity", "cycle_closure") %in% names(features)))
    expect_equal(features$split, c("train", "train", "test"))
})

test_that("AMLSim timestamp tie-breaking preserves real temporal gaps", {
    t <- amlsim_strict_unit_time(c(0, 0, 10, 10, 10, 40))

    expect_true(all(diff(t) > 0))
    expect_gt(t[3] - t[2], t[2] - t[1])
    expect_gt(t[6] - t[5], t[4] - t[3])
})

test_that("AMLSim reader merges separate alert account labels", {
    out_dir <- tempdir()
    tx_path <- file.path(out_dir, "tx.csv")
    alert_path <- file.path(out_dir, "alerts.csv")

    utils::write.csv(data.frame(
        TXN_ID = 1:3,
        ACCOUNT_ID = c("A", "B", "C"),
        COUNTER_PARTY_ACCOUNT_NUM = c("B", "C", "A"),
        TXN_SOURCE_TYPE_CODE = "WIRE",
        TXN_AMOUNT_ORIG = c(100, 200, 300),
        start = 1:3,
        stringsAsFactors = FALSE
    ), tx_path, row.names = FALSE)
    utils::write.csv(data.frame(
        ALERT_KEY = 1,
        ACCOUNT_ID = "C",
        Escalated_To_Case_Investigation = "YES",
        stringsAsFactors = FALSE
    ), alert_path, row.names = FALSE)

    standardized <- amlsim_read_transactions(out_dir)
    expect_equal(standardized$is_sar, c(FALSE, TRUE, TRUE))
})

test_that("AMLSim graph-generator transactions are supported without Java output", {
    out_dir <- tempdir()
    tx_path <- file.path(out_dir, "transactions.csv")
    alert_path <- file.path(out_dir, "alert_members.csv")

    utils::write.csv(data.frame(
        id = 1:4,
        src = c(1, 2, 3, 4),
        dst = c(2, 3, 4, 1),
        ttype = "TRANSFER"
    ), tx_path, row.names = FALSE)
    utils::write.csv(data.frame(
        alertID = 1,
        accountID = 3,
        isSAR = "true"
    ), alert_path, row.names = FALSE)

    standardized <- amlsim_read_transactions(out_dir)
    expect_equal(standardized$amount, rep(1, 4))
    expect_true(all(diff(standardized$t) > 0))
    expect_equal(standardized$is_sar, c(FALSE, TRUE, TRUE, FALSE))
})

test_that("AMLSim SAR metrics handle ranked scores", {
    predictions <- data.frame(
        event_id = 1:4,
        score = c(0.1, 0.9, 0.8, 0.2),
        is_sar = c(FALSE, TRUE, TRUE, FALSE),
        model = "toy",
        stringsAsFactors = FALSE
    )

    metrics <- amlsim_evaluate_sar(predictions)
    expect_equal(metrics$model, "toy")
    expect_equal(metrics$roc_auc, 1)
    expect_equal(metrics$average_precision, 1)
})
