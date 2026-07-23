eth_stablecoin_registry <- function() {
    data.frame(
        symbol = c("USDC", "USDT", "DAI"),
        token_address = tolower(c(
            "0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48",
            "0xdAC17F958D2ee523a2206206994597C13D831ec7",
            "0x6B175474E89094C44Da98b954EedeAC495271d0F"
        )),
        decimals = c(6L, 6L, 18L),
        stringsAsFactors = FALSE
    )
}

eth_normalize_address <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x[x %in% c("", "na", "nan", "null", "none")] <- NA_character_
    x
}

eth_pick_col <- function(dat, candidates, required = TRUE) {
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

eth_strict_unit_time <- function(time_num) {
    if (length(time_num) == 0) {
        return(numeric(0))
    }
    if (anyNA(time_num)) {
        time_num[is.na(time_num)] <- seq_len(sum(is.na(time_num))) + max(time_num, na.rm = TRUE)
    }
    if (length(time_num) == 1) {
        return(0.5)
    }

    unique_times <- sort(unique(time_num))
    positive_gaps <- diff(unique_times)
    positive_gaps <- positive_gaps[positive_gaps > 0]
    min_gap <- if (length(positive_gaps) > 0) min(positive_gaps) else 1
    tie_counts <- table(time_num)
    tie_width <- min_gap / (max(as.integer(tie_counts)) + 1)
    tie_index <- ave(seq_along(time_num), time_num, FUN = seq_along) - 1
    t <- time_num + tie_index * tie_width
    t <- t - min(t)
    if (max(t) == 0) {
        return(seq(0, 1, length.out = length(t)))
    }
    t / max(t)
}

eth_time_numeric <- function(x) {
    if (inherits(x, c("POSIXct", "POSIXt", "Date"))) {
        return(as.numeric(x))
    }
    if (is.numeric(x) || is.integer(x)) {
        return(as.numeric(x))
    }
    parsed <- suppressWarnings(as.POSIXct(
        x,
        tz = "UTC",
        tryFormats = c(
            "%Y-%m-%d %H:%M:%S",
            "%Y-%m-%dT%H:%M:%OS",
            "%Y-%m-%dT%H:%M:%OSZ",
            "%Y-%m-%d"
        )
    ))
    if (sum(!is.na(parsed)) >= length(x) / 2) {
        return(as.numeric(parsed))
    }
    suppressWarnings(as.numeric(x))
}

eth_build_stablecoin_sql <- function(start_date,
                                     end_date,
                                     tokens = eth_stablecoin_registry(),
                                     min_amount = 0,
                                     max_rows = NULL,
                                     table = "`bigquery-public-data.crypto_ethereum.token_transfers`") {
    token_addresses <- paste(sprintf("'%s'", tolower(tokens$token_address)), collapse = ", ")
    limit_clause <- if (!is.null(max_rows)) paste("LIMIT", as.integer(max_rows)) else ""
    paste0(
        "SELECT\n",
        "  block_timestamp,\n",
        "  transaction_hash,\n",
        "  log_index,\n",
        "  block_number,\n",
        "  LOWER(token_address) AS token_address,\n",
        "  LOWER(from_address) AS from_address,\n",
        "  LOWER(to_address) AS to_address,\n",
        "  SAFE_CAST(value AS BIGNUMERIC) AS raw_value\n",
        "FROM ", table, "\n",
        "WHERE block_timestamp >= TIMESTAMP('", start_date, "')\n",
        "  AND block_timestamp < TIMESTAMP('", end_date, "')\n",
        "  AND LOWER(token_address) IN (", token_addresses, ")\n",
        "  AND from_address IS NOT NULL\n",
        "  AND to_address IS NOT NULL\n",
        "  AND from_address != to_address\n",
        "  AND SAFE_CAST(value AS BIGNUMERIC) > ", format(min_amount, scientific = FALSE), "\n",
        "ORDER BY block_timestamp, block_number, log_index\n",
        limit_clause
    )
}

eth_query_bigquery <- function(start_date,
                               end_date,
                               project = Sys.getenv("GCP_PROJECT", unset = ""),
                               billing = project,
                               tokens = eth_stablecoin_registry(),
                               min_amount = 0,
                               max_rows = NULL) {
    if (!requireNamespace("bigrquery", quietly = TRUE)) {
        stop("Install `bigrquery` to query BigQuery, or provide an exported CSV to eth_read_transfers().")
    }
    if (!nzchar(billing)) {
        stop("Set `GCP_PROJECT` or pass `billing` to bill BigQuery public-data queries.")
    }
    sql <- eth_build_stablecoin_sql(
        start_date = start_date,
        end_date = end_date,
        tokens = tokens,
        min_amount = min_amount,
        max_rows = max_rows
    )
    job <- bigrquery::bq_project_query(billing, sql)
    as.data.frame(bigrquery::bq_table_download(job))
}

eth_read_transfers <- function(path) {
    dat <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = "character")
    eth_standardize_transfers(dat)
}

eth_standardize_transfers <- function(dat, tokens = eth_stablecoin_registry()) {
    time_col <- eth_pick_col(dat, c("block_timestamp", "timestamp", "time", "datetime"))
    tx_col <- eth_pick_col(dat, c("transaction_hash", "tx_hash", "hash"), required = FALSE)
    log_col <- eth_pick_col(dat, c("log_index", "evt_index", "event_index"), required = FALSE)
    token_col <- eth_pick_col(dat, c("token_address", "contract_address", "token"))
    from_col <- eth_pick_col(dat, c("from_address", "from", "sender", "src"))
    to_col <- eth_pick_col(dat, c("to_address", "to", "receiver", "dst"))
    value_col <- eth_pick_col(dat, c("raw_value", "value", "amount_raw", "amount"))

    token_map <- tokens
    token_map$token_address <- eth_normalize_address(token_map$token_address)
    dat_token <- eth_normalize_address(dat[[token_col]])
    decimals <- token_map$decimals[match(dat_token, token_map$token_address)]
    symbol <- token_map$symbol[match(dat_token, token_map$token_address)]
    decimals[is.na(decimals)] <- 0

    raw_value <- suppressWarnings(as.numeric(dat[[value_col]]))
    amount <- raw_value / (10 ^ decimals)
    time_num <- eth_time_numeric(dat[[time_col]])
    event_id <- if (!is.null(tx_col)) {
        paste0(dat[[tx_col]], ":", if (!is.null(log_col)) dat[[log_col]] else seq_len(nrow(dat)))
    } else {
        as.character(seq_len(nrow(dat)))
    }

    out <- data.frame(
        event_id = as.character(event_id),
        t_raw = as.character(dat[[time_col]]),
        time_num = time_num,
        from_address = eth_normalize_address(dat[[from_col]]),
        to_address = eth_normalize_address(dat[[to_col]]),
        token_address = dat_token,
        token_symbol = ifelse(is.na(symbol), dat_token, symbol),
        raw_value = raw_value,
        amount = amount,
        transaction_hash = if (!is.null(tx_col)) as.character(dat[[tx_col]]) else NA_character_,
        log_index = if (!is.null(log_col)) suppressWarnings(as.integer(dat[[log_col]])) else seq_len(nrow(dat)),
        stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$from_address) & !is.na(out$to_address) &
                   out$from_address != out$to_address &
                   !is.na(out$time_num) & !is.na(out$amount) & out$amount > 0, , drop = FALSE]
    out <- out[order(out$time_num, out$log_index, seq_len(nrow(out))), , drop = FALSE]
    out$t <- eth_strict_unit_time(out$time_num)
    rownames(out) <- NULL
    out
}

eth_default_label_url <- function() {
    "https://raw.githubusercontent.com/brianleect/etherscan-labels/main/data/etherscan/combined/combinedAllLabels.json"
}

eth_download_labels <- function(dest = "ethereum_labels.json",
                                url = eth_default_label_url()) {
    utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
    dest
}

eth_read_labels <- function(path) {
    if (grepl("\\.json$", path, ignore.case = TRUE)) {
        if (!requireNamespace("jsonlite", quietly = TRUE)) {
            stop("Install `jsonlite` to read JSON label files, or provide a CSV label file.")
        }
        raw <- jsonlite::fromJSON(path, simplifyDataFrame = FALSE)
        addresses <- names(raw)
        entities <- vapply(raw, function(x) {
            if (!is.null(x$name) && nzchar(x$name)) x$name else NA_character_
        }, character(1))
        categories <- vapply(raw, function(x) {
            labels <- x$labels
            if (length(labels) == 0) NA_character_ else paste(labels, collapse = ";")
        }, character(1))
        out <- data.frame(
            address = eth_normalize_address(addresses),
            entity = ifelse(is.na(entities), eth_normalize_address(addresses), entities),
            category = categories,
            stringsAsFactors = FALSE
        )
        out <- out[!is.na(out$address), , drop = FALSE]
        return(out[!duplicated(out$address), , drop = FALSE])
    }

    labels <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                              colClasses = "character")
    address_col <- eth_pick_col(labels, c("address", "Address", "account", "contract"), required = TRUE)
    label_col <- eth_pick_col(labels, c("label", "name", "Label", "Name"), required = FALSE)
    category_col <- eth_pick_col(labels, c("category", "labels", "tag", "type", "LabelType"), required = FALSE)

    out <- data.frame(
        address = eth_normalize_address(labels[[address_col]]),
        entity = if (!is.null(label_col)) as.character(labels[[label_col]]) else eth_normalize_address(labels[[address_col]]),
        category = if (!is.null(category_col)) as.character(labels[[category_col]]) else NA_character_,
        stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$address), , drop = FALSE]
    out[!duplicated(out$address), , drop = FALSE]
}

eth_builtin_entity_labels <- function() {
    data.frame(
        address = eth_normalize_address(c(
            "0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48",
            "0xdAC17F958D2ee523a2206206994597C13D831ec7",
            "0x6B175474E89094C44Da98b954EedeAC495271d0F",
            "0x28c6c06298d514db089934071355e5743bf21d60",
            "0x21a31ee1afc51d94c2efccaa2092ad1028285549",
            "0xdfd5293d8e347dfe59e90efd55b2956a1343963d",
            "0x3f5ce5fbfe3e9af3971dD833D26bA9b5C936f0bE",
            "0x503828976d22510aad0201ac7ec88293211d23da",
            "0xbe0eb53f46cd790cd13851d5eff43d12404d33e8"
        )),
        entity = c(
            "USDC contract",
            "USDT contract",
            "DAI contract",
            "Binance 14",
            "Binance 15",
            "Binance 16",
            "Binance wallet",
            "Coinbase 1",
            "Binance 7"
        ),
        category = c(
            "stablecoin_contract",
            "stablecoin_contract",
            "stablecoin_contract",
            rep("exchange", 6)
        ),
        stringsAsFactors = FALSE
    )
}

eth_apply_entity_labels <- function(transfers,
                                    labels = NULL,
                                    keep_unlabeled = TRUE,
                                    unlabeled_prefix = "addr") {
    labels <- if (is.null(labels)) eth_builtin_entity_labels() else labels
    labels$address <- eth_normalize_address(labels$address)
    labels <- labels[!duplicated(labels$address), , drop = FALSE]

    from_match <- match(transfers$from_address, labels$address)
    to_match <- match(transfers$to_address, labels$address)
    from_entity <- labels$entity[from_match]
    to_entity <- labels$entity[to_match]
    from_category <- labels$category[from_match]
    to_category <- labels$category[to_match]

    if (keep_unlabeled) {
        from_entity[is.na(from_entity)] <- paste0(unlabeled_prefix, ":", transfers$from_address[is.na(from_entity)])
        to_entity[is.na(to_entity)] <- paste0(unlabeled_prefix, ":", transfers$to_address[is.na(to_entity)])
        from_category[is.na(from_category)] <- "unlabeled"
        to_category[is.na(to_category)] <- "unlabeled"
    }

    out <- transfers
    out$from_entity <- from_entity
    out$to_entity <- to_entity
    out$from_category <- from_category
    out$to_category <- to_category
    out <- out[!is.na(out$from_entity) & !is.na(out$to_entity) & out$from_entity != out$to_entity, , drop = FALSE]
    rownames(out) <- NULL
    out
}

eth_filter_active_entities <- function(transfers,
                                       max_entities = 5000,
                                       max_events = NULL,
                                       prefer_labeled = FALSE) {
    if (!is.null(max_entities)) {
        entities <- sort(unique(c(transfers$from_entity, transfers$to_entity)))
        activity <- stats::setNames(rep(0, length(entities)), entities)
        volume <- stats::setNames(rep(0, length(entities)), entities)
        counts <- table(c(transfers$from_entity, transfers$to_entity))
        activity[names(counts)] <- as.numeric(counts)
        out_vol <- tapply(transfers$amount, transfers$from_entity, sum, na.rm = TRUE)
        in_vol <- tapply(transfers$amount, transfers$to_entity, sum, na.rm = TRUE)
        volume[names(out_vol)] <- volume[names(out_vol)] + as.numeric(out_vol)
        volume[names(in_vol)] <- volume[names(in_vol)] + as.numeric(in_vol)
        score <- activity + log1p(volume)
        if (prefer_labeled) {
            categories <- tapply(
                c(transfers$from_category, transfers$to_category),
                c(transfers$from_entity, transfers$to_entity),
                function(x) x[which(!is.na(x))[1]]
            )
            labeled_entities <- names(categories)[!is.na(categories) & categories != "unlabeled"]
            score[labeled_entities] <- score[labeled_entities] + max(score) + 1
        }
        keep <- names(sort(score, decreasing = TRUE))[seq_len(min(max_entities, length(score)))]
        transfers <- transfers[transfers$from_entity %in% keep & transfers$to_entity %in% keep, , drop = FALSE]
    }
    transfers <- transfers[order(transfers$time_num, transfers$log_index), , drop = FALSE]
    if (!is.null(max_events) && nrow(transfers) > max_events) {
        transfers <- transfers[seq_len(max_events), , drop = FALSE]
    }
    transfers$t <- eth_strict_unit_time(transfers$time_num)
    rownames(transfers) <- NULL
    transfers
}

eth_make_repeated_hits <- function(transfers,
                                   weight = c("log_amount", "amount", "count")) {
    weight <- match.arg(weight)
    entities <- sort(unique(c(transfers$from_entity, transfers$to_entity)))
    actor_map <- data.frame(
        actor_id = seq_along(entities),
        entity = entities,
        stringsAsFactors = FALSE
    )
    idx <- stats::setNames(actor_map$actor_id, actor_map$entity)
    hit_weight <- switch(
        weight,
        log_amount = log1p(transfers$amount),
        amount = transfers$amount,
        count = rep(1, nrow(transfers))
    )
    if (weight != "count" && max(hit_weight, na.rm = TRUE) > 0) {
        hit_weight <- hit_weight / max(hit_weight, na.rm = TRUE)
    }
    hits <- data.frame(
        t = transfers$t,
        i = as.integer(idx[transfers$from_entity]),
        j = as.integer(idx[transfers$to_entity]),
        weight = hit_weight,
        amount = transfers$amount,
        token_symbol = transfers$token_symbol,
        from_entity = transfers$from_entity,
        to_entity = transfers$to_entity,
        event_id = transfers$event_id,
        stringsAsFactors = FALSE
    )
    attr(hits, "actor_map") <- actor_map
    hits
}

eth_node_covariates <- function(transfers) {
    entities <- sort(unique(c(transfers$from_entity, transfers$to_entity)))
    out_count <- table(transfers$from_entity)
    in_count <- table(transfers$to_entity)
    out_amount <- tapply(transfers$amount, transfers$from_entity, sum, na.rm = TRUE)
    in_amount <- tapply(transfers$amount, transfers$to_entity, sum, na.rm = TRUE)
    node_entity <- c(transfers$from_entity, transfers$to_entity)
    node_time <- c(transfers$time_num, transfers$time_num)
    first_time <- tapply(node_time, node_entity, min, na.rm = TRUE)
    last_time <- tapply(node_time, node_entity, max, na.rm = TRUE)
    category <- tapply(
        c(transfers$from_category, transfers$to_category),
        node_entity,
        function(x) {
            x <- x[!is.na(x)]
            if (length(x) == 0) NA_character_ else names(sort(table(x), decreasing = TRUE))[1]
        }
    )

    out <- data.frame(
        entity = entities,
        category = unname(category[entities]),
        out_count = as.numeric(out_count[entities]),
        in_count = as.numeric(in_count[entities]),
        out_amount = as.numeric(out_amount[entities]),
        in_amount = as.numeric(in_amount[entities]),
        first_time_num = as.numeric(first_time[entities]),
        last_time_num = as.numeric(last_time[entities]),
        stringsAsFactors = FALSE
    )
    numeric_cols <- c("out_count", "in_count", "out_amount", "in_amount")
    for (col in numeric_cols) out[[col]][is.na(out[[col]])] <- 0
    out$total_count <- out$out_count + out$in_count
    out$total_amount <- out$out_amount + out$in_amount
    out$net_amount <- out$in_amount - out$out_amount
    out$activity_share <- out$total_count / max(1, sum(out$total_count))
    out$active_span_seconds <- out$last_time_num - out$first_time_num
    out
}

eth_temporal_split <- function(hits, train_frac = 0.7) {
    stopifnot(train_frac > 0, train_frac < 1)
    cut <- max(1, floor(nrow(hits) * train_frac))
    split <- rep("test", nrow(hits))
    split[seq_len(cut)] <- "train"
    split
}

eth_default_rhem_params <- function() {
    c(edgeValue = 0,
      recipValue = 0,
      senderValueActivity = 0,
      receiverValueActivity = 0,
      transitiveValue = 0,
      cycleValue = 0,
      commonSourceValue = 0,
      commonTargetValue = 0)
}

eth_default_formula <- function() {
    paste(names(eth_default_rhem_params()), collapse = " + ")
}

eth_timenet_formula <- function(beta = 8,
                                terms = c("edge", "recip", "sender", "receiver")) {
    term_map <- c(
        edge = "decayedEdgeValue",
        recip = "decayedRecipValue",
        sender = "decayedSenderActivity",
        receiver = "decayedReceiverActivity",
        transitive = "decayedTransitiveValue",
        cycle = "decayedCycleValue",
        common_source = "decayedCommonSourceValue",
        common_target = "decayedCommonTargetValue",
        age = "timeSinceLastEdge"
    )
    terms <- match.arg(terms, choices = names(term_map), several.ok = TRUE)
    pieces <- vapply(terms, function(term) {
        name <- unname(term_map[term])
        if (term == "age") {
            paste0(name, "()")
        } else {
            paste0(name, "(beta = ", format(beta, scientific = FALSE), ")")
        }
    }, character(1))
    paste(pieces, collapse = " + ")
}

eth_timenet_formula_grid <- function(betas = c(0, 2, 8, 32)) {
    patterns <- list(
        edge_recip = c("edge", "recip"),
        edge_recip_activity = c("edge", "recip", "sender", "receiver"),
        edge_recip_triads = c("edge", "recip", "transitive", "cycle"),
        edge_activity_shared = c("edge", "sender", "receiver", "common_source", "common_target"),
        full = c("edge", "recip", "sender", "receiver", "transitive", "cycle",
                 "common_source", "common_target")
    )
    rows <- list()
    pos <- 1L
    for (beta in betas) {
        for (pattern in names(patterns)) {
            rows[[pos]] <- data.frame(
                model_id = paste0("timenet_", pattern, "_b", beta),
                beta = beta,
                formula_RHS = eth_timenet_formula(beta, patterns[[pattern]]),
                stringsAsFactors = FALSE
            )
            pos <- pos + 1L
        }
    }
    do.call(rbind, rows)
}

eth_simple_rhem_params <- function() {
    c(repetition = 0,
      reciprocity = 0,
      sender_activity = 0,
      receiver_activity = 0)
}

eth_initial_params <- function(hits) {
    n_events <- max(1, nrow(hits))
    duration <- max(diff(range(hits$t)), 1e-6)
    list(mu = n_events / duration,
         beta_overall = 1,
         K = 0.1,
         beta_edges = 0,
         RHEM_params = eth_default_rhem_params())
}

eth_initial_params_for_formula <- function(hits, formula_RHS = eth_default_formula()) {
    params <- eth_initial_params(hits)
    if (is.null(formula_RHS)) {
        params$RHEM_params <- eth_simple_rhem_params()
    } else {
        params$RHEM_params <- stats::setNames(rep(0, length(eth_formula_param_names(formula_RHS))),
                                              eth_formula_param_names(formula_RHS))
    }
    params
}

eth_formula_param_names <- function(formula_RHS) {
    if (is.null(formula_RHS)) {
        return(names(eth_simple_rhem_params()))
    }
    time_terms <- c("decayedEdgeValue", "decayedRecipValue",
                    "decayedSenderActivity", "decayedReceiverActivity",
                    "decayedTransitiveValue", "decayedCycleValue",
                    "decayedCommonSourceValue", "decayedCommonTargetValue",
                    "edgeAge", "timeSinceLastEdge")
    formula_text <- if (inherits(formula_RHS, "formula")) {
        as.character(formula_RHS)[length(formula_RHS)]
    } else {
        as.character(formula_RHS)
    }
    if (any(vapply(time_terms, function(term) grepl(paste0("\\b", term, "\\b"), formula_text), logical(1)))) {
        if (!requireNamespace("timeNet", quietly = TRUE)) {
            stop("The `timeNet` package is required for temporal Ethereum model formulas.")
        }
        return(timeNet::parse_time_formula(formula_RHS)$term)
    }
    attr(stats::terms(stats::as.formula(paste("~", formula_text))), "term.labels")
}

eth_require_hawkes_deps <- function() {
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

eth_fit_hawkes <- function(hits,
                           actors = sort(unique(c(hits$i, hits$j))),
                           formula_RHS = eth_default_formula(),
                           maxit = 25,
                           trace = 0,
                           cache_max_risk_dyads = 10000,
                           streaming_risk_threshold = 100000,
                           fixed_params = "beta_edges") {
    eth_require_hawkes_deps()
    train_hits <- hits[hits$split == "train", , drop = FALSE]
    if (nrow(train_hits) < 2) {
        stop("Need at least two training hits to fit valued HawkesNet.")
    }

    params_init <- eth_initial_params_for_formula(train_hits, formula_RHS)
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

eth_fit_no_excitation_rhem <- function(hits,
                                       actors = sort(unique(c(hits$i, hits$j))),
                                       formula_RHS = eth_default_formula(),
                                       maxit = 25,
                                       trace = 0,
                                       cache_max_risk_dyads = 10000,
                                       streaming_risk_threshold = 100000,
                                       fixed_params = c("K", "beta_overall", "beta_edges")) {
    eth_require_hawkes_deps()
    train_hits <- hits[hits$split == "train", , drop = FALSE]
    if (nrow(train_hits) < 2) {
        stop("Need at least two training hits to fit no-excitation RHEM.")
    }

    params_init <- eth_initial_params_for_formula(train_hits, formula_RHS)
    params_init$K <- 0
    params_init$beta_overall <- 1
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
        fixed_params = fixed_params,
        trace = trace,
        maxit = maxit,
        get_hessian = FALSE
    )
    fit$params_init <- params_init
    fit$optimized_skeleton <- optimized_skeleton
    fit$formula_RHS <- formula_RHS
    fit$rhem_stats_cache_size <- if (is.null(rhem_stats_cache)) NA_integer_ else length(ls(rhem_stats_cache))
    fit$fixed_params <- fixed_params
    fit$model <- "rhem_no_excitation_hawkesnet"
    fit
}

.eth_risk_set <- function(actors) {
    risk <- expand.grid(i = actors, j = actors,
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
    risk <- risk[risk$i != risk$j, , drop = FALSE]
    rownames(risk) <- NULL
    risk
}

.eth_dyad_key <- function(i, j) {
    paste(i, j, sep = "\r")
}

.eth_logsumexp <- function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
}

.eth_prob_score_row <- function(hits, idx, risk, probs, model) {
    true_idx <- which(risk$i == hits$i[idx] & risk$j == hits$j[idx])
    if (length(true_idx) != 1) {
        true_prob <- NA_real_
        true_rank <- NA_real_
    } else {
        true_prob <- probs[true_idx]
        tie_tol <- sqrt(.Machine$double.eps)
        greater <- sum(probs > true_prob + tie_tol)
        tied <- sum(abs(probs - true_prob) <= tie_tol)
        true_rank <- greater + (tied + 1) / 2
    }
    data.frame(
        event_index = idx,
        event_id = hits$event_id[idx],
        score = -log(true_prob),
        log_mark_density = log(true_prob),
        true_prob = true_prob,
        true_rank = true_rank,
        token_symbol = hits$token_symbol[idx],
        from_entity = hits$from_entity[idx],
        to_entity = hits$to_entity[idx],
        model = model,
        stringsAsFactors = FALSE
    )
}

eth_score_hawkes_marks <- function(hits,
                                   fit,
                                   actors = sort(unique(c(hits$i, hits$j))),
                                   formula_RHS = eth_default_formula(),
                                   streaming_risk_threshold = 100000,
                                   score_type = c("full_intensity", "mark"),
                                   model = "valued_hawkesnet_exact") {
    eth_require_hawkes_deps()
    score_type <- match.arg(score_type)
    estimated <- utils::relist(fit$fit$par, skeleton = fit$optimized_skeleton)
    params <- fit$params_init
    params[names(estimated)] <- estimated
    params$RHEM_params <- stats::setNames(as.numeric(params$RHEM_params), names(fit$params_init$RHEM_params))

    test_idx <- which(hits$split == "test")
    risk <- .eth_risk_set(actors)
    rows <- vector("list", length(test_idx))
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
        if (is.null(pmf$edge_probs)) {
            row <- data.frame(
                event_index = idx,
                event_id = hits$event_id[idx],
                score = -pmf$log_mark_density,
                log_mark_density = pmf$log_mark_density,
                true_prob = exp(pmf$log_mark_density),
                true_rank = NA_integer_,
                token_symbol = hits$token_symbol[idx],
                from_entity = hits$from_entity[idx],
                to_entity = hits$to_entity[idx],
                model = model,
                stringsAsFactors = FALSE
            )
        } else {
            row <- .eth_prob_score_row(hits, idx, risk, pmf$edge_probs, model = model)
        }
        if (score_type == "full_intensity") {
            prior_times <- hits$t[seq_len(idx - 1)]
            decays <- exp(-params$beta_overall * (hits$t[idx] - prior_times))
            ground_intensity <- params$mu + params$K * sum(decays)
            row$score <- -(pmf$log_mark_density + log(ground_intensity))
        }
        rows[[idx_pos]] <- row
    }

    do.call(rbind, rows)
}

eth_summarize_mark_scores <- function(scores) {
    models <- unique(scores$model)
    rows <- lapply(models, function(model) {
        x <- scores[scores$model == model, , drop = FALSE]
        data.frame(
            model = model,
            n = nrow(x),
            mean_negative_log_mark = mean(-x$log_mark_density, na.rm = TRUE),
            median_negative_log_mark = stats::median(-x$log_mark_density, na.rm = TRUE),
            mean_score = mean(x$score, na.rm = TRUE),
            mean_rank = mean(x$true_rank, na.rm = TRUE),
            mrr = mean(1 / x$true_rank, na.rm = TRUE),
            top1 = mean(x$true_rank <= 1, na.rm = TRUE),
            top10 = mean(x$true_rank <= 10, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

eth_fit_score_timenet_candidate <- function(hits,
                                            actors,
                                            candidate,
                                            kind = c("rhem", "hawkes"),
                                            maxit = 10,
                                            streaming_risk_threshold = 100000,
                                            cache_max_risk_dyads = 10000,
                                            trace = 0) {
    kind <- match.arg(kind)
    formula_RHS <- candidate$formula_RHS[1]
    model_id <- paste0(kind, "_", candidate$model_id[1])
    fit_fun <- if (kind == "rhem") eth_fit_no_excitation_rhem else eth_fit_hawkes
    fit <- fit_fun(
        hits,
        actors = actors,
        formula_RHS = formula_RHS,
        maxit = maxit,
        trace = trace,
        cache_max_risk_dyads = cache_max_risk_dyads,
        streaming_risk_threshold = streaming_risk_threshold
    )
    scores <- eth_score_hawkes_marks(
        hits,
        fit,
        actors = actors,
        formula_RHS = formula_RHS,
        streaming_risk_threshold = streaming_risk_threshold,
        score_type = "mark",
        model = model_id
    )
    metrics <- eth_summarize_mark_scores(scores)
    metrics$kind <- kind
    metrics$model_id <- candidate$model_id[1]
    metrics$beta <- candidate$beta[1]
    metrics$formula_RHS <- formula_RHS
    list(fit = fit, scores = scores, metrics = metrics)
}

eth_search_timenet_models <- function(hits,
                                      actors = sort(unique(c(hits$i, hits$j))),
                                      candidates = eth_timenet_formula_grid(),
                                      include_baselines = TRUE,
                                      rhem_maxit = 10,
                                      hawkes_maxit = 10,
                                      hawkes_top_n = 3,
                                      streaming_risk_threshold = 100000,
                                      cache_max_risk_dyads = 10000,
                                      history_maxit = 50,
                                      trace = 0,
                                      output_dir = NULL) {
    if (!requireNamespace("timeNet", quietly = TRUE)) {
        stop("Install the local `timeNet` package before running temporal Ethereum model search.")
    }
    scores <- list()
    metrics <- list()
    fits <- list()
    errors <- list()
    timings <- list()

    if (include_baselines) {
        timings$baselines <- system.time({
            baseline_cmp <- eth_run_baseline_comparisons(
                hits,
                actors = actors,
                history_maxit = history_maxit
            )
        })
        scores$baselines <- baseline_cmp$scores
        metrics$baselines <- baseline_cmp$metrics
        fits$history_softmax <- baseline_cmp$history_fit
    }

    for (row in seq_len(nrow(candidates))) {
        candidate <- candidates[row, , drop = FALSE]
        key <- paste0("rhem_", candidate$model_id)
        timings[[key]] <- system.time({
            result <- tryCatch(
                eth_fit_score_timenet_candidate(
                    hits,
                    actors = actors,
                    candidate = candidate,
                    kind = "rhem",
                    maxit = rhem_maxit,
                    streaming_risk_threshold = streaming_risk_threshold,
                    cache_max_risk_dyads = cache_max_risk_dyads,
                    trace = trace
                ),
                error = function(e) e
            )
        })
        if (inherits(result, "error")) {
            errors[[key]] <- conditionMessage(result)
        } else {
            scores[[key]] <- result$scores
            metrics[[key]] <- result$metrics
            fits[[key]] <- result$fit
        }
    }

    rhem_metrics <- do.call(rbind, metrics[grepl("^rhem_", names(metrics))])
    if (!is.null(rhem_metrics) && nrow(rhem_metrics) > 0 && hawkes_top_n > 0) {
        rhem_metrics <- rhem_metrics[order(rhem_metrics$mean_negative_log_mark,
                                           -rhem_metrics$mrr,
                                           rhem_metrics$mean_rank), , drop = FALSE]
        top_ids <- utils::head(rhem_metrics$model_id, hawkes_top_n)
        for (model_id in top_ids) {
            candidate <- candidates[candidates$model_id == model_id, , drop = FALSE][1, , drop = FALSE]
            key <- paste0("hawkes_", candidate$model_id)
            timings[[key]] <- system.time({
                result <- tryCatch(
                    eth_fit_score_timenet_candidate(
                        hits,
                        actors = actors,
                        candidate = candidate,
                        kind = "hawkes",
                        maxit = hawkes_maxit,
                        streaming_risk_threshold = streaming_risk_threshold,
                        cache_max_risk_dyads = cache_max_risk_dyads,
                        trace = trace
                    ),
                    error = function(e) e
                )
            })
            if (inherits(result, "error")) {
                errors[[key]] <- conditionMessage(result)
            } else {
                scores[[key]] <- result$scores
                metrics[[key]] <- result$metrics
                fits[[key]] <- result$fit
            }
        }
    }

    combined_scores <- if (length(scores) > 0) do.call(rbind, scores) else NULL
    combined_metrics <- if (length(metrics) > 0) do.call(rbind, metrics) else NULL
    if (!is.null(combined_metrics)) {
        combined_metrics <- combined_metrics[order(combined_metrics$mean_negative_log_mark,
                                                   -combined_metrics$mrr,
                                                   combined_metrics$mean_rank), , drop = FALSE]
        rownames(combined_metrics) <- NULL
    }
    out <- list(scores = combined_scores,
                metrics = combined_metrics,
                fits = fits,
                errors = errors,
                timings = timings,
                candidates = candidates)
    if (!is.null(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        if (!is.null(combined_scores)) {
            utils::write.csv(combined_scores, file.path(output_dir, "ethereum_timenet_model_scores.csv"), row.names = FALSE)
        }
        if (!is.null(combined_metrics)) {
            utils::write.csv(combined_metrics, file.path(output_dir, "ethereum_timenet_model_metrics.csv"), row.names = FALSE)
        }
        utils::write.csv(candidates, file.path(output_dir, "ethereum_timenet_candidates.csv"), row.names = FALSE)
        saveRDS(fits, file.path(output_dir, "ethereum_timenet_fits.rds"))
        saveRDS(errors, file.path(output_dir, "ethereum_timenet_errors.rds"))
        saveRDS(timings, file.path(output_dir, "ethereum_timenet_timings.rds"))
    }
    out
}

eth_score_frequency_baseline <- function(hits,
                                         actors = sort(unique(c(hits$i, hits$j))),
                                         alpha = 1,
                                         model = "dyad_frequency") {
    test_idx <- which(hits$split == "test")
    risk <- .eth_risk_set(actors)
    risk_key <- .eth_dyad_key(risk$i, risk$j)
    rows <- vector("list", length(test_idx))
    for (idx_pos in seq_along(test_idx)) {
        idx <- test_idx[idx_pos]
        history <- hits[seq_len(idx - 1), , drop = FALSE]
        counts <- table(.eth_dyad_key(history$i, history$j))
        vals <- rep(alpha, nrow(risk))
        match_idx <- match(names(counts), risk_key, nomatch = 0)
        keep <- match_idx > 0
        vals[match_idx[keep]] <- vals[match_idx[keep]] + as.numeric(counts[keep])
        probs <- vals / sum(vals)
        rows[[idx_pos]] <- .eth_prob_score_row(hits, idx, risk, probs, model)
    }
    do.call(rbind, rows)
}

eth_score_recency_baseline <- function(hits,
                                       actors = sort(unique(c(hits$i, hits$j))),
                                       decay = 10,
                                       alpha = 1e-3,
                                       model = "dyad_recency") {
    test_idx <- which(hits$split == "test")
    risk <- .eth_risk_set(actors)
    risk_key <- .eth_dyad_key(risk$i, risk$j)
    rows <- vector("list", length(test_idx))
    for (idx_pos in seq_along(test_idx)) {
        idx <- test_idx[idx_pos]
        history <- hits[seq_len(idx - 1), , drop = FALSE]
        vals <- rep(alpha, nrow(risk))
        if (nrow(history) > 0) {
            decayed <- history$weight * exp(-decay * (hits$t[idx] - history$t))
            totals <- tapply(decayed, .eth_dyad_key(history$i, history$j), sum)
            match_idx <- match(names(totals), risk_key, nomatch = 0)
            keep <- match_idx > 0
            vals[match_idx[keep]] <- vals[match_idx[keep]] + as.numeric(totals[keep])
        }
        probs <- vals / sum(vals)
        rows[[idx_pos]] <- .eth_prob_score_row(hits, idx, risk, probs, model)
    }
    do.call(rbind, rows)
}

eth_score_activity_baseline <- function(hits,
                                        actors = sort(unique(c(hits$i, hits$j))),
                                        alpha = 1,
                                        model = "sender_receiver_activity") {
    test_idx <- which(hits$split == "test")
    risk <- .eth_risk_set(actors)
    rows <- vector("list", length(test_idx))
    for (idx_pos in seq_along(test_idx)) {
        idx <- test_idx[idx_pos]
        history <- hits[seq_len(idx - 1), , drop = FALSE]
        sender_counts <- table(history$i)
        receiver_counts <- table(history$j)
        sender <- rep(alpha, length(actors))
        receiver <- rep(alpha, length(actors))
        names(sender) <- names(receiver) <- as.character(actors)
        sender[names(sender_counts)] <- sender[names(sender_counts)] + as.numeric(sender_counts)
        receiver[names(receiver_counts)] <- receiver[names(receiver_counts)] + as.numeric(receiver_counts)
        vals <- sender[as.character(risk$i)] * receiver[as.character(risk$j)]
        probs <- vals / sum(vals)
        rows[[idx_pos]] <- .eth_prob_score_row(hits, idx, risk, probs, model)
    }
    do.call(rbind, rows)
}

.eth_history_stats <- function(risk, history, actors) {
    edge_counts <- table(.eth_dyad_key(history$i, history$j))
    recip_counts <- edge_counts[.eth_dyad_key(risk$j, risk$i)]
    dyad_counts <- edge_counts[.eth_dyad_key(risk$i, risk$j)]
    sender_counts <- table(history$i)
    receiver_counts <- table(history$j)
    stats <- cbind(
        repetition = log1p(as.numeric(dyad_counts)),
        reciprocity = log1p(as.numeric(recip_counts)),
        sender_activity = log1p(as.numeric(sender_counts[as.character(risk$i)])),
        receiver_activity = log1p(as.numeric(receiver_counts[as.character(risk$j)]))
    )
    stats[is.na(stats)] <- 0
    stats
}

eth_fit_history_softmax <- function(hits,
                                    actors = sort(unique(c(hits$i, hits$j))),
                                    maxit = 50,
                                    trace = 0) {
    train_idx <- which(hits$split == "train")
    risk <- .eth_risk_set(actors)
    fn <- function(theta) {
        ll <- 0
        for (idx in train_idx) {
            history <- hits[seq_len(idx - 1), , drop = FALSE]
            stats <- .eth_history_stats(risk, history, actors)
            eta <- as.vector(stats %*% theta)
            true_idx <- which(risk$i == hits$i[idx] & risk$j == hits$j[idx])
            ll <- ll + eta[true_idx] - .eth_logsumexp(eta)
        }
        ll
    }
    fit <- stats::optim(
        par = rep(0, length(eth_simple_rhem_params())),
        fn = fn,
        method = "Nelder-Mead",
        control = list(fnscale = -1, trace = trace, maxit = maxit),
        hessian = FALSE
    )
    names(fit$par) <- names(eth_simple_rhem_params())
    list(fit = fit, actors = actors)
}

eth_score_history_softmax <- function(hits,
                                      fit,
                                      actors = fit$actors,
                                      model = "history_softmax") {
    test_idx <- which(hits$split == "test")
    risk <- .eth_risk_set(actors)
    rows <- vector("list", length(test_idx))
    for (idx_pos in seq_along(test_idx)) {
        idx <- test_idx[idx_pos]
        history <- hits[seq_len(idx - 1), , drop = FALSE]
        stats <- .eth_history_stats(risk, history, actors)
        eta <- as.vector(stats %*% fit$fit$par)
        probs <- exp(eta - .eth_logsumexp(eta))
        rows[[idx_pos]] <- .eth_prob_score_row(hits, idx, risk, probs, model)
    }
    do.call(rbind, rows)
}

eth_run_baseline_comparisons <- function(hits,
                                         actors = sort(unique(c(hits$i, hits$j))),
                                         history_maxit = 50) {
    history_fit <- eth_fit_history_softmax(hits, actors = actors, maxit = history_maxit)
    scores <- rbind(
        eth_score_frequency_baseline(hits, actors = actors),
        eth_score_recency_baseline(hits, actors = actors),
        eth_score_activity_baseline(hits, actors = actors),
        eth_score_history_softmax(hits, history_fit, actors = actors)
    )
    list(scores = scores,
         metrics = eth_summarize_mark_scores(scores),
         history_fit = history_fit)
}

eth_prepare_hawkesnet_data <- function(transfers,
                                       labels = NULL,
                                       max_entities = 1000,
                                       max_events = 10000,
                                       weight = "log_amount",
                                       train_frac = 0.7,
                                       prefer_labeled = FALSE) {
    labeled <- eth_apply_entity_labels(transfers, labels = labels, keep_unlabeled = TRUE)
    filtered <- eth_filter_active_entities(
        labeled,
        max_entities = max_entities,
        max_events = max_events,
        prefer_labeled = prefer_labeled
    )
    hits <- eth_make_repeated_hits(filtered, weight = weight)
    hits$split <- eth_temporal_split(hits, train_frac = train_frac)
    covariates <- eth_node_covariates(filtered)
    list(
        transfers = filtered,
        hits = hits,
        actor_map = attr(hits, "actor_map"),
        node_covariates = covariates
    )
}
