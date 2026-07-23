## Sourced by study scripts / SLURM without requiring a full package install.
## Prefer setting HAWKESNET_OUTPUT_DIR on NeSI.

hawkesnet_resolve_output_dir <- function(pkg_root, create = TRUE) {
  env <- Sys.getenv("HAWKESNET_OUTPUT_DIR", unset = "")
  if (nzchar(env)) {
    out <- normalizePath(env, mustWork = FALSE)
  } else {
    cur <- normalizePath(pkg_root, mustWork = FALSE)
    out <- NULL
    for (i in seq_len(6L)) {
      parent <- dirname(cur)
      sibling <- file.path(parent, "cluster_output")
      if (dir.exists(sibling)) {
        out <- normalizePath(sibling, mustWork = FALSE)
        break
      }
      if (identical(parent, cur)) break
      cur <- parent
    }
    if (is.null(out)) {
      out <- normalizePath(file.path(dirname(pkg_root), "cluster_output"), mustWork = FALSE)
    }
  }
  if (isTRUE(create)) {
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    for (sd in c("runs", "logs", "paper_figures", "diagnostics")) {
      dir.create(file.path(out, sd), showWarnings = FALSE, recursive = TRUE)
    }
  }
  out
}
