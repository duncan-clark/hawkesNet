amlsim_first_existing <- function(paths) {
    paths <- paths[file.exists(paths)]
    if (length(paths) == 0) {
        return(NULL)
    }
    paths[1]
}

amlsim_require_hawkes_deps <- function() {
    if (!requireNamespace("ernm", quietly = TRUE)) {
        stop("The `ernm` package is required for exact valued HawkesNet fitting.")
    }
    if (!requireNamespace("network", quietly = TRUE)) {
        stop("The `network` package is required for exact valued HawkesNet fitting.")
    }
    suppressPackageStartupMessages(require("ernm", character.only = TRUE))
    suppressPackageStartupMessages(require("network", character.only = TRUE))
    invisible(TRUE)
}

amlsim_pick_col <- function(dat, candidates, required = TRUE) {
    nms <- names(dat)
    idx <- match(tolower(candidates), tolower(nms), nomatch = 0)
    idx <- idx[idx > 0]
    if (length(idx) > 0) {
        return(nms[idx[1]])
    }
    if (required) {
        stop("Could not find any required column: ", paste(candidates, collapse = ", "))
    }
    NULL
}

amlsim_parse_flag <- function(x) {
    if (is.null(x)) {
        return(logical(0))
    }
    if (is.logical(x)) {
        return(x)
    }
    tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

amlsim_time_numeric <- function(x) {
    if (inherits(x, c("POSIXct", "POSIXt", "Date"))) {
        return(as.numeric(x))
    }
    if (is.numeric(x) || is.integer(x)) {
        return(as.numeric(x))
    }

    numeric_x <- suppressWarnings(as.numeric(x))
    if (sum(!is.na(numeric_x)) >= length(x) / 2) {
        return(numeric_x)
    }

    parsed <- suppressWarnings(as.POSIXct(
        x,
        tz = "UTC",
        tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d", "%Y%m%d", "%Y-%m-%dT%H:%M:%S")
    ))
    as.numeric(parsed)
}

amlsim_strict_unit_time <- function(time_num) {
    if (length(time_num) == 0) {
        return(numeric(0))
    }
    if (anyNA(time_num)) {
        time_num[is.na(time_num)] <- seq_len(sum(is.na(time_num))) + max(time_num, na.rm = TRUE)
    }
    if (length(time_num) == 1) {
        return(0.5)
    }

    # AMLSim can emit many transactions in the same simulation step. Preserve
    # real gaps between steps/dates, and only spread ties within the smallest
    # observed gap so the Hawkes likelihood sees strict event times.
    unique_times <- sort(unique(time_num))
    gaps <- diff(unique_times)
    positive_gaps <- gaps[gaps > 0]
    min_gap <- if (length(positive_gaps) > 0) min(positive_gaps) else 1
    tie_counts <- table(time_num)
    max_ties <- max(as.integer(tie_counts))
    tie_width <- min_gap / (max_ties + 1)
    tie_index <- ave(seq_along(time_num), time_num, FUN = seq_along) - 1
    t <- time_num + tie_index * tie_width
    t <- t - min(t)
    if (max(t) == 0) {
        return(seq(0, 1, length.out = length(t)))
    }
    t / max(t)
}

amlsim_resolve_transactions_path <- function(path) {
    if (missing(path) || is.null(path)) {
        stop("`path` must be an AMLSim transaction CSV file or output directory.")
    }
    if (!dir.exists(path)) {
        return(path)
    }
    candidates <- file.path(path, c(
        "transactions.csv",
        "tx.csv",
        "tx_log.csv",
        "alert_transactions.csv"
    ))
    resolved <- amlsim_first_existing(candidates)
    if (is.null(resolved)) {
        stop("No AMLSim transaction CSV found in directory: ", path)
    }
    resolved
}

amlsim_read_transactions <- function(path) {
    tx_path <- amlsim_resolve_transactions_path(path)
    dat <- utils::read.csv(tx_path, stringsAsFactors = FALSE, check.names = FALSE)
    tx <- amlsim_standardize_transactions(dat, source_path = tx_path)
    output_dir <- if (dir.exists(path)) path else dirname(tx_path)
    amlsim_apply_external_alert_labels(tx, output_dir)
}

amlsim_standardize_transactions <- function(dat, source_path = NA_character_) {
    nms_lower <- tolower(names(dat))

    if (all(c("step", "nameorig", "namedest", "amount") %in% nms_lower)) {
        id_col <- amlsim_pick_col(dat, c("tran_id", "tx_id", "TXN_ID"), required = FALSE)
        from_col <- amlsim_pick_col(dat, c("nameOrig"))
        to_col <- amlsim_pick_col(dat, c("nameDest"))
        amount_col <- amlsim_pick_col(dat, c("amount"))
        time_col <- amlsim_pick_col(dat, c("step"))
        type_col <- amlsim_pick_col(dat, c("type"), required = FALSE)
        sar_col <- amlsim_pick_col(dat, c("isSAR"), required = FALSE)
        alert_col <- amlsim_pick_col(dat, c("alertID"), required = FALSE)
    } else if (all(c("id", "src", "dst", "ttype") %in% nms_lower)) {
        id_col <- amlsim_pick_col(dat, c("id"))
        from_col <- amlsim_pick_col(dat, c("src"))
        to_col <- amlsim_pick_col(dat, c("dst"))
        amount_col <- NULL
        time_col <- NULL
        type_col <- amlsim_pick_col(dat, c("ttype"), required = FALSE)
        sar_col <- NULL
        alert_col <- NULL
    } else if (all(c("account_id", "counter_party_account_num", "txn_amount_orig") %in% nms_lower)) {
        id_col <- amlsim_pick_col(dat, c("TXN_ID"), required = FALSE)
        from_col <- amlsim_pick_col(dat, c("ACCOUNT_ID"))
        to_col <- amlsim_pick_col(dat, c("COUNTER_PARTY_ACCOUNT_NUM"))
        amount_col <- amlsim_pick_col(dat, c("TXN_AMOUNT_ORIG"))
        time_col <- amlsim_pick_col(dat, c("start", "end"))
        type_col <- amlsim_pick_col(dat, c("TXN_SOURCE_TYPE_CODE"), required = FALSE)
        sar_col <- NULL
        alert_col <- NULL
    } else {
        id_col <- amlsim_pick_col(dat, c("tran_id", "tx_id", "TXN_ID"), required = FALSE)
        from_col <- amlsim_pick_col(dat, c("orig_acct", "orig_id", "nameOrig", "ACCOUNT_ID"))
        to_col <- amlsim_pick_col(dat, c("bene_acct", "dest_acct", "dest_id", "nameDest", "COUNTER_PARTY_ACCOUNT_NUM"))
        amount_col <- amlsim_pick_col(dat, c("base_amt", "amount", "TXN_AMOUNT_ORIG"), required = FALSE)
        time_col <- amlsim_pick_col(dat, c("tran_timestamp", "timestamp", "step", "start", "end"), required = FALSE)
        type_col <- amlsim_pick_col(dat, c("tx_type", "transaction_type", "type", "TXN_SOURCE_TYPE_CODE"), required = FALSE)
        sar_col <- amlsim_pick_col(dat, c("is_sar", "isSAR", "sar_flag"), required = FALSE)
        alert_col <- amlsim_pick_col(dat, c("alert_id", "alertID"), required = FALSE)
    }

    event_id <- if (!is.null(id_col)) dat[[id_col]] else seq_len(nrow(dat))
    time_raw <- if (!is.null(time_col)) dat[[time_col]] else seq_len(nrow(dat))
    time_num <- amlsim_time_numeric(time_raw)
    amount <- if (!is.null(amount_col)) suppressWarnings(as.numeric(dat[[amount_col]])) else rep(1, nrow(dat))
    amount[is.na(amount)] <- 0

    is_sar <- if (!is.null(sar_col)) {
        amlsim_parse_flag(dat[[sar_col]])
    } else {
        rep(FALSE, nrow(dat))
    }
    alert_id <- if (!is.null(alert_col)) as.character(dat[[alert_col]]) else rep(NA_character_, nrow(dat))
    alert_present <- !is.na(alert_id) & !(alert_id %in% c("", "-1", "NA", "NaN"))
    is_sar <- is_sar | alert_present

    out <- data.frame(
        event_id = as.character(event_id),
        t_raw = as.character(time_raw),
        time_num = time_num,
        from_id = as.character(dat[[from_col]]),
        to_id = as.character(dat[[to_col]]),
        amount = amount,
        tx_type = if (!is.null(type_col)) as.character(dat[[type_col]]) else NA_character_,
        is_sar = is_sar,
        alert_id = alert_id,
        source_path = source_path,
        stringsAsFactors = FALSE
    )
    out <- out[out$from_id != out$to_id, , drop = FALSE]
    out <- out[order(out$time_num, seq_len(nrow(out))), , drop = FALSE]
    out$t <- amlsim_strict_unit_time(out$time_num)
    rownames(out) <- NULL
    out
}

amlsim_apply_external_alert_labels <- function(transactions, output_dir) {
    if (!dir.exists(output_dir)) {
        return(transactions)
    }

    alert_tx_path <- amlsim_first_existing(file.path(output_dir, c("alert_transactions.csv", "alert_tx.csv")))
    if (!is.null(alert_tx_path)) {
        alert_dat <- utils::read.csv(alert_tx_path, stringsAsFactors = FALSE, check.names = FALSE)
        alert_tx <- tryCatch(
            amlsim_standardize_transactions(alert_dat, source_path = alert_tx_path),
            error = function(e) NULL
        )
        if (!is.null(alert_tx)) {
            alert_ids <- unique(alert_tx$event_id)
            alert_keys <- unique(paste(alert_tx$from_id, alert_tx$to_id, alert_tx$t_raw, sep = "\r"))
            tx_keys <- paste(transactions$from_id, transactions$to_id, transactions$t_raw, sep = "\r")
            transactions$is_sar <- transactions$is_sar |
                transactions$event_id %in% alert_ids |
                tx_keys %in% alert_keys
        }
    }

    alert_accounts_path <- amlsim_first_existing(file.path(output_dir, c("sar_accounts.csv", "alert_accounts.csv", "alert_members.csv", "alerts.csv")))
    if (!is.null(alert_accounts_path)) {
        acct_dat <- utils::read.csv(alert_accounts_path, stringsAsFactors = FALSE, check.names = FALSE)
        acct_col <- amlsim_pick_col(acct_dat, c("acct_id", "account_id", "accountID", "ACCOUNT_ID"), required = FALSE)
        escalated_col <- amlsim_pick_col(acct_dat, c("is_sar", "isSAR", "Escalated_To_Case_Investigation"), required = FALSE)
        if (!is.null(acct_col)) {
            keep <- rep(TRUE, nrow(acct_dat))
            if (!is.null(escalated_col)) {
                keep <- amlsim_parse_flag(acct_dat[[escalated_col]])
            }
            alert_accounts <- unique(as.character(acct_dat[[acct_col]][keep]))
            transactions$is_sar <- transactions$is_sar |
                transactions$from_id %in% alert_accounts |
                transactions$to_id %in% alert_accounts
        }
    }

    transactions
}

amlsim_make_repeated_hits <- function(transactions,
                                      weight = c("log_amount", "amount", "count"),
                                      rescale_weights = TRUE) {
    weight <- match.arg(weight)
    actors <- sort(unique(c(transactions$from_id, transactions$to_id)))

    hit_weight <- switch(weight,
        log_amount = log1p(transactions$amount),
        amount = transactions$amount,
        count = rep(1, nrow(transactions))
    )
    hit_weight[is.na(hit_weight) | hit_weight < 0] <- 0
    if (rescale_weights && any(hit_weight > 0)) {
        hit_weight <- hit_weight / mean(hit_weight[hit_weight > 0])
    }

    hits <- data.frame(
        t = transactions$t,
        i = match(transactions$from_id, actors),
        j = match(transactions$to_id, actors),
        weight = hit_weight,
        is_sar = transactions$is_sar,
        amount = transactions$amount,
        tx_type = transactions$tx_type,
        from_id = transactions$from_id,
        to_id = transactions$to_id,
        event_id = transactions$event_id,
        stringsAsFactors = FALSE
    )
    attr(hits, "actor_map") <- data.frame(
        actor = seq_along(actors),
        account_id = actors,
        stringsAsFactors = FALSE
    )
    hits
}

amlsim_filter_active_actors <- function(hits, max_actors = NULL, max_events = NULL) {
    if (!is.null(max_actors)) {
        actor_ids <- sort(unique(c(hits$i, hits$j)))
        activity <- stats::setNames(rep(0, length(actor_ids)), actor_ids)
        tab <- table(c(hits$i, hits$j))
        activity[names(tab)] <- as.numeric(tab)

        sar_actors <- unique(c(hits$i[hits$is_sar], hits$j[hits$is_sar]))
        activity[as.character(sar_actors)] <- activity[as.character(sar_actors)] + max(activity) + 1
        keep_actors <- as.integer(names(sort(activity, decreasing = TRUE))[seq_len(min(max_actors, length(activity)))])
        hits <- hits[hits$i %in% keep_actors & hits$j %in% keep_actors, , drop = FALSE]
    }
    hits <- hits[order(hits$t), , drop = FALSE]
    if (!is.null(max_events) && nrow(hits) > max_events) {
        hits <- hits[seq_len(max_events), , drop = FALSE]
    }
    rownames(hits) <- NULL
    hits$t <- amlsim_strict_unit_time(hits$t)
    hits
}

amlsim_temporal_split <- function(hits, train_frac = 0.7) {
    stopifnot(train_frac > 0, train_frac < 1)
    cut <- max(1, floor(nrow(hits) * train_frac))
    split <- rep("test", nrow(hits))
    split[seq_len(cut)] <- "train"
    split
}

amlsim_default_rhem_params <- function() {
    c(edgeValue = 0,
      recipValue = 0,
      senderValueActivity = 0,
      receiverValueActivity = 0,
      transitiveValue = 0,
      cycleValue = 0,
      commonSourceValue = 0,
      commonTargetValue = 0)
}

amlsim_default_formula <- function() {
    paste(names(amlsim_default_rhem_params()), collapse = " + ")
}

amlsim_initial_params <- function(hits) {
    n_events <- max(1, nrow(hits))
    duration <- max(diff(range(hits$t)), 1e-6)
    list(mu = n_events / duration,
         beta_overall = 1,
         K = 0.1,
         beta_edges = 0,
         RHEM_params = amlsim_default_rhem_params())
}

amlsim_fit_hawkes <- function(hits,
                              actors = sort(unique(c(hits$i, hits$j))),
                              formula_RHS = amlsim_default_formula(),
                              train_normals_only = TRUE,
                              maxit = 50,
                              trace = 0,
                              cache_max_risk_dyads = 10000,
                              streaming_risk_threshold = 100000,
                              fixed_params = "beta_edges") {
    amlsim_require_hawkes_deps()
    train_hits <- hits[hits$split == "train", , drop = FALSE]
    if (train_normals_only && "is_sar" %in% names(train_hits)) {
        train_hits <- train_hits[!train_hits$is_sar, , drop = FALSE]
    }
    if (nrow(train_hits) < 2) {
        stop("Need at least two training hits to fit valued HawkesNet.")
    }

    params_init <- amlsim_initial_params(train_hits)
    fixed_params <- intersect(fixed_params, names(params_init))
    optimized_skeleton <- params_init
    optimized_skeleton[fixed_params] <- NULL
    risk_dyads <- length(actors) * (length(actors) - 1)
    rhem_stats_cache <- if (!is.null(cache_max_risk_dyads) && risk_dyads <= cache_max_risk_dyads) {
        new.env(parent = emptyenv())
    } else {
        NULL
    }
    fit <- hawkesNet::fit_hawkesNet(
        params_init = params_init,
        time_window = c(0, 1),
        mark_filtration = train_hits[, c("t", "i", "j", "weight"), drop = FALSE],
        PMF_mark = hawkesNet::PMF_mark,
        type = "RHEM",
        actors = actors,
        formula_RHS = formula_RHS,
        rhem_stats_cache = rhem_stats_cache,
        rhem_streaming_threshold = streaming_risk_threshold,
        fixed_params = if (length(fixed_params) == 0) NULL else fixed_params,
        trace = trace,
        maxit = maxit,
        get_hessian = FALSE
    )
    fit$params_init <- params_init
    fit$optimized_skeleton <- optimized_skeleton
    fit$formula_RHS <- formula_RHS
    fit$rhem_stats_cache_size <- if (is.null(rhem_stats_cache)) NA_integer_ else length(ls(rhem_stats_cache))
    fit$fixed_params <- fixed_params
    fit
}

amlsim_score_hawkes_sar <- function(hits,
                                    fit,
                                    actors = sort(unique(c(hits$i, hits$j))),
                                    formula_RHS = amlsim_default_formula(),
                                    streaming_risk_threshold = 100000,
                                    score_type = c("full_intensity", "mark")) {
    amlsim_require_hawkes_deps()
    score_type <- match.arg(score_type)
    estimated <- utils::relist(fit$fit$par, skeleton = fit$optimized_skeleton)
    params <- fit$params_init
    params[names(estimated)] <- estimated
    params$RHEM_params <- stats::setNames(as.numeric(params$RHEM_params), names(fit$params_init$RHEM_params))

    test_idx <- which(hits$split == "test")
    scores <- rep(NA_real_, length(test_idx))
    for (idx_pos in seq_along(test_idx)) {
        idx <- test_idx[idx_pos]
        history <- hits[seq_len(idx), c("t", "i", "j", "weight"), drop = FALSE]
        pmf <- hawkesNet::PMF_mark_RHEM(
            time = hits$t[idx],
            params = params,
            mark_filtration = history,
            actors = actors,
            formula_RHS = formula_RHS,
            rhem_streaming_threshold = streaming_risk_threshold
        )
        if (score_type == "full_intensity") {
            prior_times <- hits$t[seq_len(idx - 1)]
            decays <- exp(-params$beta_overall * (hits$t[idx] - prior_times))
            ground_intensity <- params$mu + params$K * sum(decays)
            scores[idx_pos] <- -(pmf$log_mark_density + log(ground_intensity))
        } else {
            scores[idx_pos] <- -pmf$log_mark_density
        }
    }

    data.frame(
        event_id = hits$event_id[test_idx],
        score = scores,
        is_sar = hits$is_sar[test_idx],
        model = "valued_hawkesnet_exact",
        stringsAsFactors = FALSE
    )
}

amlsim_history_features <- function(hits) {
    actors <- sort(unique(c(hits$i, hits$j)))
    n <- length(actors)
    actor_index <- stats::setNames(seq_along(actors), actors)
    state <- matrix(0, nrow = n, ncol = n)
    row_sum <- numeric(n)
    col_sum <- numeric(n)
    n_hits <- nrow(hits)

    repetition <- reciprocity <- sender_activity <- receiver_activity <- numeric(n_hits)
    transitive_closure <- cycle_closure <- common_source <- common_target <- numeric(n_hits)

    for (row in seq_len(n_hits)) {
        i <- actor_index[as.character(hits$i[row])]
        j <- actor_index[as.character(hits$j[row])]

        repetition[row] <- state[i, j]
        reciprocity[row] <- state[j, i]
        sender_activity[row] <- row_sum[i]
        receiver_activity[row] <- col_sum[j]
        transitive_closure[row] <- sum(state[i, ] * state[, j])
        cycle_closure[row] <- sum(state[j, ] * state[, i])
        common_source[row] <- row_sum[i] - state[i, j]
        common_target[row] <- col_sum[j] - state[i, j]

        weight <- hits$weight[row]
        state[i, j] <- state[i, j] + weight
        row_sum[i] <- row_sum[i] + weight
        col_sum[j] <- col_sum[j] + weight
    }

    data.frame(
        event_id = hits$event_id,
        is_sar = hits$is_sar,
        split = hits$split,
        log_amount = log1p(hits$amount),
        repetition = repetition,
        reciprocity = reciprocity,
        sender_activity = sender_activity,
        receiver_activity = receiver_activity,
        transitive_closure = transitive_closure,
        cycle_closure = cycle_closure,
        common_source = common_source,
        common_target = common_target,
        stringsAsFactors = FALSE
    )
}

amlsim_fit_glm_baseline <- function(features) {
    train <- features[features$split == "train", , drop = FALSE]
    x_cols <- setdiff(names(features), c("event_id", "is_sar", "split"))
    if (length(unique(train$is_sar)) < 2) {
        return(list(type = "constant", rate = mean(train$is_sar), x_cols = x_cols))
    }
    stats::glm(
        stats::as.formula(paste("is_sar ~", paste(x_cols, collapse = " + "))),
        data = train,
        family = stats::binomial()
    )
}

amlsim_predict_glm_baseline <- function(model, features) {
    test <- features[features$split == "test", , drop = FALSE]
    if (is.list(model) && identical(model$type, "constant")) {
        scores <- rep(model$rate, nrow(test))
    } else {
        scores <- stats::predict(model, newdata = test, type = "response")
    }
    data.frame(
        event_id = test$event_id,
        score = as.numeric(scores),
        is_sar = test$is_sar,
        model = "glm_history_features",
        stringsAsFactors = FALSE
    )
}

amlsim_fit_rpart_baseline <- function(features) {
    if (!requireNamespace("rpart", quietly = TRUE)) {
        return(NULL)
    }
    train <- features[features$split == "train", , drop = FALSE]
    if (length(unique(train$is_sar)) < 2) {
        return(NULL)
    }
    x_cols <- setdiff(names(features), c("event_id", "is_sar", "split"))
    train$is_sar <- factor(train$is_sar)
    rpart::rpart(
        stats::as.formula(paste("is_sar ~", paste(x_cols, collapse = " + "))),
        data = train,
        method = "class"
    )
}

amlsim_predict_rpart_baseline <- function(model, features) {
    if (is.null(model)) {
        return(NULL)
    }
    test <- features[features$split == "test", , drop = FALSE]
    pred <- stats::predict(model, newdata = test, type = "prob")
    yes_col <- intersect(colnames(pred), "TRUE")
    scores <- if (length(yes_col) == 1) pred[, yes_col] else rep(0, nrow(test))
    data.frame(
        event_id = test$event_id,
        score = as.numeric(scores),
        is_sar = test$is_sar,
        model = "rpart_history_features",
        stringsAsFactors = FALSE
    )
}

amlsim_auc_roc <- function(labels, scores) {
    labels <- as.logical(labels)
    keep <- !is.na(labels) & !is.na(scores)
    labels <- labels[keep]
    scores <- scores[keep]
    n_pos <- sum(labels)
    n_neg <- sum(!labels)
    if (n_pos == 0 || n_neg == 0) {
        return(NA_real_)
    }
    ranks <- rank(scores, ties.method = "average")
    (sum(ranks[labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

amlsim_average_precision <- function(labels, scores) {
    labels <- as.logical(labels)
    keep <- !is.na(labels) & !is.na(scores)
    labels <- labels[keep]
    scores <- scores[keep]
    n_pos <- sum(labels)
    if (n_pos == 0) {
        return(NA_real_)
    }
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    precision <- cumsum(labels) / seq_along(labels)
    sum(precision[labels]) / n_pos
}

amlsim_evaluate_sar <- function(predictions) {
    models <- unique(predictions$model)
    do.call(rbind, lapply(models, function(model_name) {
        pred <- predictions[predictions$model == model_name, , drop = FALSE]
        data.frame(
            model = model_name,
            n = nrow(pred),
            n_sar = sum(pred$is_sar, na.rm = TRUE),
            roc_auc = amlsim_auc_roc(pred$is_sar, pred$score),
            average_precision = amlsim_average_precision(pred$is_sar, pred$score),
            stringsAsFactors = FALSE
        )
    }))
}

amlsim_run_sar_study <- function(path,
                                 max_actors = 100,
                                 max_events = 1000,
                                 train_frac = 0.7,
                                 run_hawkes = TRUE,
                                 run_ml = TRUE,
                                 hawkes_maxit = 50,
                                 weight = "log_amount",
                                 cache_max_risk_dyads = 10000,
                                 streaming_risk_threshold = 100000,
                                 fixed_hawkes_params = "beta_edges",
                                 hawkes_score_type = c("full_intensity", "mark")) {
    hawkes_score_type <- match.arg(hawkes_score_type)
    tx <- amlsim_read_transactions(path)
    hits <- amlsim_make_repeated_hits(tx, weight = weight)
    hits <- amlsim_filter_active_actors(hits, max_actors = max_actors, max_events = max_events)
    hits$split <- amlsim_temporal_split(hits, train_frac = train_frac)
    actors <- sort(unique(c(hits$i, hits$j)))

    predictions <- list()
    timings <- list()
    errors <- list()

    if (run_hawkes) {
        timings$hawkes_fit <- system.time({
            hawkes_fit <- tryCatch(
                amlsim_fit_hawkes(hits, actors = actors, maxit = hawkes_maxit,
                                  cache_max_risk_dyads = cache_max_risk_dyads,
                                  streaming_risk_threshold = streaming_risk_threshold,
                                  fixed_params = fixed_hawkes_params),
                error = function(e) e
            )
        })
        if (inherits(hawkes_fit, "error")) {
            errors$hawkes <- conditionMessage(hawkes_fit)
        } else {
            timings$hawkes_score <- system.time({
                predictions$hawkes <- amlsim_score_hawkes_sar(
                    hits, hawkes_fit, actors = actors,
                    streaming_risk_threshold = streaming_risk_threshold,
                    score_type = hawkes_score_type
                )
            })
        }
    }

    if (run_ml) {
        features <- amlsim_history_features(hits)
        glm_fit <- amlsim_fit_glm_baseline(features)
        predictions$glm <- amlsim_predict_glm_baseline(glm_fit, features)

        rpart_fit <- amlsim_fit_rpart_baseline(features)
        predictions$rpart <- amlsim_predict_rpart_baseline(rpart_fit, features)
    }

    predictions <- predictions[!vapply(predictions, is.null, logical(1))]
    predictions <- if (length(predictions) > 0) do.call(rbind, predictions) else data.frame()
    list(
        summary = data.frame(
            n_events = nrow(hits),
            n_actors = length(actors),
            n_train = sum(hits$split == "train"),
            n_test = sum(hits$split == "test"),
            n_sar = sum(hits$is_sar),
            exact_risk_dyads = length(actors) * (length(actors) - 1),
            stringsAsFactors = FALSE
        ),
        metrics = if (nrow(predictions) > 0) amlsim_evaluate_sar(predictions) else data.frame(),
        predictions = predictions,
        timings = timings,
        errors = errors,
        hits = hits
    )
}

amlsim_run_exact_scalability_grid <- function(path,
                                              actor_grid = c(25, 50, 100),
                                              event_grid = c(250, 500, 1000),
                                              hawkes_maxit = 5) {
    grid <- expand.grid(max_actors = actor_grid, max_events = event_grid)
    results <- vector("list", nrow(grid))
    for (idx in seq_len(nrow(grid))) {
        elapsed <- system.time({
            result <- tryCatch(
                amlsim_run_sar_study(
                    path = path,
                    max_actors = grid$max_actors[idx],
                    max_events = grid$max_events[idx],
                    run_hawkes = TRUE,
                    run_ml = FALSE,
                    hawkes_maxit = hawkes_maxit
                ),
                error = function(e) e
            )
        })
        if (inherits(result, "error")) {
            results[[idx]] <- data.frame(
                grid[idx, ],
                elapsed = elapsed[["elapsed"]],
                n_events = NA_integer_,
                n_actors = NA_integer_,
                exact_risk_dyads = NA_integer_,
                error = conditionMessage(result)
            )
        } else {
            results[[idx]] <- data.frame(
                grid[idx, ],
                elapsed = elapsed[["elapsed"]],
                result$summary[, c("n_events", "n_actors", "exact_risk_dyads")],
                error = if (length(result$errors)) paste(unlist(result$errors), collapse = "; ") else NA_character_
            )
        }
    }
    do.call(rbind, results)
}
