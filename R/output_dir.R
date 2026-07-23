#' Resolve the hawkesNet run-output directory
#'
#' Cluster studies and SLURM jobs should write RDS, logs, and paper figures
#' outside the package git checkout. Resolution order:
#' \enumerate{
#'   \item Environment variable \code{HAWKESNET_OUTPUT_DIR} (recommended on NeSI)
#'   \item An existing \code{cluster_output/} directory found by walking up from
#'     \code{pkg_root} (typical laptop sync folder)
#'   \item Sibling path \code{file.path(dirname(pkg_root), "cluster_output")}
#'     (typical NeSI layout: \code{.../hawkes_net/hawkesNet} + \code{.../hawkes_net/cluster_output})
#' }
#'
#' @param pkg_root Character path to the package root (directory with DESCRIPTION).
#'   If \code{NULL}, uses \code{getwd()} when it looks like a package root.
#' @param create Logical; create the directory if missing (default \code{TRUE}).
#' @param subdirs Optional character vector of subdirectories to ensure exist
#'   (e.g. \code{c("runs", "logs", "paper_figures", "diagnostics")}).
#' @return Normalized path to the output directory (invisibly creates it).
#' @export
hawkesnet_output_dir <- function(pkg_root = NULL, create = TRUE, subdirs = NULL) {
  if (is.null(pkg_root) || !nzchar(pkg_root)) {
    pkg_root <- getwd()
  }
  pkg_root <- normalizePath(pkg_root, mustWork = FALSE)

  env <- Sys.getenv("HAWKESNET_OUTPUT_DIR", unset = "")
  if (nzchar(env)) {
    out <- normalizePath(env, mustWork = FALSE)
  } else {
    out <- .hawkesnet_find_output_dir(pkg_root)
  }

  if (isTRUE(create)) {
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    if (!is.null(subdirs)) {
      for (sd in subdirs) {
        dir.create(file.path(out, sd), showWarnings = FALSE, recursive = TRUE)
      }
    }
  }
  normalizePath(out, mustWork = FALSE)
}

#' @rdname hawkesnet_output_dir
#' @export
hawkesnet_paper_figures_dir <- function(pkg_root = NULL, create = TRUE) {
  file.path(
    hawkesnet_output_dir(pkg_root = pkg_root, create = create, subdirs = "paper_figures"),
    "paper_figures"
  )
}

#' @noRd
.hawkesnet_find_output_dir <- function(pkg_root) {
  cur <- pkg_root
  for (i in seq_len(6L)) {
    parent <- dirname(cur)
    sibling <- file.path(parent, "cluster_output")
    if (dir.exists(sibling)) {
      return(normalizePath(sibling, mustWork = FALSE))
    }
    if (identical(parent, cur)) break
    cur <- parent
  }
  # Preferred default when nothing exists yet: sibling of package root
  normalizePath(file.path(dirname(pkg_root), "cluster_output"), mustWork = FALSE)
}
