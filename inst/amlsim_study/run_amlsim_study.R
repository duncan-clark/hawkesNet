library(hawkesNet)

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
this_dir <- if (!is.na(this_file)) dirname(this_file) else file.path(getwd(), "inst", "amlsim_study")
source(file.path(this_dir, "amlsim_utils.R"))

resolve_default_amlsim_path <- function() {
    env_path <- Sys.getenv("AMLSIM_TX_PATH", unset = NA_character_)
    if (!is.na(env_path) && nzchar(env_path)) {
        return(env_path)
    }

    env_dir <- Sys.getenv("AMLSIM_DIR", unset = NA_character_)
    if (!is.na(env_dir) && nzchar(env_dir)) {
        return(file.path(env_dir, "outputs", "sample"))
    }

    candidates <- c(
        file.path(getwd(), "external", "AMLSim", "outputs", "sample"),
        file.path(getwd(), "..", "external", "AMLSim", "outputs", "sample"),
        file.path(getwd(), "..", "..", "external", "AMLSim", "outputs", "sample"),
        file.path(getwd(), "external", "AMLSim", "sample", "outputs"),
        file.path(getwd(), "..", "external", "AMLSim", "sample", "outputs"),
        file.path(getwd(), "..", "..", "external", "AMLSim", "sample", "outputs")
    )
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) {
        stop("Set AMLSIM_TX_PATH to a transaction CSV or AMLSIM_DIR to an AMLSim checkout.")
    }
    found[1]
}

tx_path <- resolve_default_amlsim_path()
output_dir <- Sys.getenv("AMLSIM_STUDY_OUTPUT_DIR", unset = file.path(getwd(), "amlsim_study_results"))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

max_actors <- as.integer(Sys.getenv("AMLSIM_MAX_ACTORS", unset = "50"))
max_events <- as.integer(Sys.getenv("AMLSIM_MAX_EVENTS", unset = "500"))
hawkes_maxit <- as.integer(Sys.getenv("AMLSIM_HAWKES_MAXIT", unset = "25"))

message("Running AMLSim SAR study")
message("  transactions: ", tx_path)
message("  max actors:   ", max_actors)
message("  max events:   ", max_events)
message("  Hawkes maxit: ", hawkes_maxit)

study <- amlsim_run_sar_study(
    path = tx_path,
    max_actors = max_actors,
    max_events = max_events,
    hawkes_maxit = hawkes_maxit
)

print(study$summary)
print(study$metrics)
if (length(study$errors) > 0) {
    print(study$errors)
}

saveRDS(study, file.path(output_dir, "amlsim_sar_study.rds"))
utils::write.csv(study$metrics, file.path(output_dir, "amlsim_sar_metrics.csv"), row.names = FALSE)
utils::write.csv(study$predictions, file.path(output_dir, "amlsim_sar_predictions.csv"), row.names = FALSE)

message("Saved AMLSim study results to: ", output_dir)
