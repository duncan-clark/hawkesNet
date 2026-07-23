library(hawkesNet)

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg) > 0) {
        this_file <- normalizePath(sub("^--file=", "", file_arg[1]))
    }
}
this_dir <- if (!is.na(this_file)) dirname(this_file) else file.path(getwd(), "inst", "ethereum_study")
source(file.path(this_dir, "ethereum_utils.R"))

parse_null_int <- function(x) {
    if (!nzchar(x) || tolower(x) %in% c("null", "none", "na")) {
        return(NULL)
    }
    as.integer(x)
}

parse_flag <- function(x) {
    tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_numeric_csv <- function(x, default) {
    if (!nzchar(x)) {
        return(default)
    }
    as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

output_dir <- Sys.getenv("ETH_STUDY_OUTPUT_DIR", unset = file.path(getwd(), "ethereum_study_results"))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

transfers_csv <- Sys.getenv("ETH_TRANSFERS_CSV", unset = "")
labels_csv <- Sys.getenv("ETH_LABELS_CSV", unset = "")
download_labels <- parse_flag(Sys.getenv("ETH_DOWNLOAD_LABELS", unset = "false"))
start_date <- Sys.getenv("ETH_START_DATE", unset = "2024-03-01")
end_date <- Sys.getenv("ETH_END_DATE", unset = "2024-03-08")
max_rows <- parse_null_int(Sys.getenv("ETH_MAX_QUERY_ROWS", unset = ""))
max_entities <- parse_null_int(Sys.getenv("ETH_MAX_ENTITIES", unset = "1000"))
max_events <- parse_null_int(Sys.getenv("ETH_MAX_EVENTS", unset = "10000"))
hawkes_maxit <- as.integer(Sys.getenv("ETH_HAWKES_MAXIT", unset = "10"))
history_maxit <- as.integer(Sys.getenv("ETH_HISTORY_MAXIT", unset = "50"))
run_hawkes <- parse_flag(Sys.getenv("ETH_RUN_HAWKES", unset = "false"))
run_rhem <- parse_flag(Sys.getenv("ETH_RUN_RHEM", unset = "false"))
run_baselines <- parse_flag(Sys.getenv("ETH_RUN_BASELINES", unset = "true"))
run_timenet_search <- parse_flag(Sys.getenv("ETH_RUN_TIMENET_SEARCH", unset = "false"))
timenet_betas <- parse_numeric_csv(Sys.getenv("ETH_TIMENET_BETAS", unset = ""), default = c(0, 2, 8, 32))
timenet_rhem_maxit <- as.integer(Sys.getenv("ETH_TIMENET_RHEM_MAXIT", unset = as.character(hawkes_maxit)))
timenet_hawkes_maxit <- as.integer(Sys.getenv("ETH_TIMENET_HAWKES_MAXIT", unset = as.character(hawkes_maxit)))
timenet_hawkes_top_n <- as.integer(Sys.getenv("ETH_TIMENET_HAWKES_TOP_N", unset = "3"))
prefer_labeled <- parse_flag(Sys.getenv("ETH_PREFER_LABELED", unset = "false"))
rhem_formula_mode <- tolower(Sys.getenv("ETH_RHEM_FORMULA", unset = "valued"))
rhem_formula <- switch(
    rhem_formula_mode,
    valued = eth_default_formula(),
    simple = NULL,
    none = NULL,
    stop("ETH_RHEM_FORMULA must be `valued` or `simple`.")
)

message("Running Ethereum stablecoin HawkesNet study")
message("  output dir:   ", output_dir)
message("  date window:  ", start_date, " to ", end_date)
message("  max entities: ", if (is.null(max_entities)) "none" else max_entities)
message("  max events:   ", if (is.null(max_events)) "none" else max_events)
message("  prefer labels:", prefer_labeled)
message("  RHEM formula: ", rhem_formula_mode)
message("  run baselines:", run_baselines)
message("  run RHEM:     ", run_rhem)
message("  run hawkes:   ", run_hawkes)
message("  timeNet search:", run_timenet_search)

if (nzchar(transfers_csv)) {
    message("Reading exported transfers: ", transfers_csv)
    transfers <- eth_read_transfers(transfers_csv)
} else {
    message("Querying BigQuery public Ethereum token_transfers table")
    transfers <- eth_query_bigquery(
        start_date = start_date,
        end_date = end_date,
        max_rows = max_rows
    )
    transfers <- eth_standardize_transfers(transfers)
    transfers_csv <- file.path(output_dir, "ethereum_stablecoin_transfers.csv")
    utils::write.csv(transfers, transfers_csv, row.names = FALSE)
    message("Saved queried transfers to: ", transfers_csv)
}

labels <- NULL
if (download_labels && !nzchar(labels_csv)) {
    labels_csv <- file.path(output_dir, "etherscan_labels.json")
    message("Downloading public Etherscan labels to: ", labels_csv)
    tryCatch(
        eth_download_labels(labels_csv),
        error = function(e) {
            warning("Could not download labels: ", conditionMessage(e))
            labels_csv <<- ""
        }
    )
}
if (nzchar(labels_csv) && file.exists(labels_csv)) {
    message("Reading labels: ", labels_csv)
    labels <- eth_read_labels(labels_csv)
}

prepared <- eth_prepare_hawkesnet_data(
    transfers = transfers,
    labels = labels,
    max_entities = max_entities,
    max_events = max_events,
    weight = "log_amount",
    prefer_labeled = prefer_labeled
)

saveRDS(prepared, file.path(output_dir, "ethereum_hawkesnet_data.rds"))
utils::write.csv(prepared$transfers, file.path(output_dir, "ethereum_transfers_prepared.csv"), row.names = FALSE)
utils::write.csv(prepared$hits, file.path(output_dir, "ethereum_hits.csv"), row.names = FALSE)
utils::write.csv(prepared$actor_map, file.path(output_dir, "ethereum_actor_map.csv"), row.names = FALSE)
utils::write.csv(prepared$node_covariates, file.path(output_dir, "ethereum_node_covariates.csv"), row.names = FALSE)

summary <- data.frame(
    n_events = nrow(prepared$hits),
    n_actors = nrow(prepared$actor_map),
    n_train = sum(prepared$hits$split == "train"),
    n_test = sum(prepared$hits$split == "test"),
    risk_dyads = nrow(prepared$actor_map) * (nrow(prepared$actor_map) - 1),
    query_start_date = start_date,
    query_end_date = end_date,
    observed_start = min(prepared$transfers$t_raw),
    observed_end = max(prepared$transfers$t_raw),
    stringsAsFactors = FALSE
)
utils::write.csv(summary, file.path(output_dir, "ethereum_summary.csv"), row.names = FALSE)
print(summary)

all_scores <- list()
all_metrics <- list()
timings <- list()
actors <- sort(unique(c(prepared$hits$i, prepared$hits$j)))

if (run_baselines) {
    timings$baselines <- system.time({
        baseline_cmp <- eth_run_baseline_comparisons(
            prepared$hits,
            actors = actors,
            history_maxit = history_maxit
        )
    })
    all_scores$baselines <- baseline_cmp$scores
    all_metrics$baselines <- baseline_cmp$metrics
    saveRDS(baseline_cmp$history_fit, file.path(output_dir, "ethereum_history_softmax_fit.rds"))
    print(baseline_cmp$metrics)
}

if (run_timenet_search) {
    timings$timenet_search <- system.time({
        timenet_candidates <- eth_timenet_formula_grid(betas = timenet_betas)
        timenet_search <- eth_search_timenet_models(
            prepared$hits,
            actors = actors,
            candidates = timenet_candidates,
            include_baselines = !run_baselines,
            rhem_maxit = timenet_rhem_maxit,
            hawkes_maxit = timenet_hawkes_maxit,
            hawkes_top_n = timenet_hawkes_top_n,
            history_maxit = history_maxit,
            output_dir = output_dir
        )
    })
    if (!is.null(timenet_search$scores)) {
        all_scores$timenet_search <- timenet_search$scores
    }
    if (!is.null(timenet_search$metrics)) {
        all_metrics$timenet_search <- timenet_search$metrics
        print(utils::head(timenet_search$metrics, 10))
    }
    if (length(timenet_search$errors) > 0) {
        warning("Some timeNet candidates failed: ", paste(names(timenet_search$errors), collapse = ", "))
    }
}

if (run_rhem) {
    timings$rhem <- system.time({
        rhem_fit <- eth_fit_no_excitation_rhem(
            prepared$hits,
            actors = actors,
            formula_RHS = rhem_formula,
            maxit = hawkes_maxit
        )
        rhem_scores <- eth_score_hawkes_marks(
            prepared$hits,
            rhem_fit,
            actors = actors,
            formula_RHS = rhem_formula,
            score_type = "full_intensity",
            model = paste0("rhem_no_excitation_", rhem_formula_mode)
        )
    })
    rhem_metrics <- eth_summarize_mark_scores(rhem_scores)
    all_scores$rhem <- rhem_scores
    all_metrics$rhem <- rhem_metrics
    saveRDS(rhem_fit, file.path(output_dir, "ethereum_rhem_no_excitation_fit.rds"))
    utils::write.csv(rhem_scores, file.path(output_dir, "ethereum_rhem_no_excitation_scores.csv"), row.names = FALSE)
    utils::write.csv(rhem_metrics, file.path(output_dir, "ethereum_rhem_no_excitation_metrics.csv"), row.names = FALSE)
    print(rhem_metrics)
}

if (run_hawkes) {
    timings$hawkes <- system.time({
        fit <- eth_fit_hawkes(
            prepared$hits,
            actors = actors,
            formula_RHS = rhem_formula,
            maxit = hawkes_maxit
        )
        scores <- eth_score_hawkes_marks(
            prepared$hits,
            fit,
            actors = actors,
            formula_RHS = rhem_formula,
            model = paste0("hawkesnet_", rhem_formula_mode)
        )
    })
    metrics <- eth_summarize_mark_scores(scores)
    all_scores$hawkes <- scores
    all_metrics$hawkes <- metrics

    saveRDS(fit, file.path(output_dir, "ethereum_hawkes_fit.rds"))
    utils::write.csv(scores, file.path(output_dir, "ethereum_hawkes_scores.csv"), row.names = FALSE)
    utils::write.csv(metrics, file.path(output_dir, "ethereum_hawkes_metrics.csv"), row.names = FALSE)
    print(metrics)
}

if (length(all_scores) > 0) {
    combined_scores <- do.call(rbind, all_scores)
    combined_metrics <- do.call(rbind, all_metrics)
    utils::write.csv(combined_scores, file.path(output_dir, "ethereum_model_scores.csv"), row.names = FALSE)
    utils::write.csv(combined_metrics, file.path(output_dir, "ethereum_model_metrics.csv"), row.names = FALSE)
}
if (length(timings) > 0) {
    saveRDS(timings, file.path(output_dir, "ethereum_model_timings.rds"))
    print(timings)
}

message("Saved Ethereum study outputs to: ", output_dir)
