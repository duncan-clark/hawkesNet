## =============================================================================
## Save GOF plots without titles for paper (hypertext results)
## =============================================================================
## Run from package root on the cluster (after hypertext job has completed):
##   Rscript inst/hypertext_conference/save_gof_plots_paper.R
##
## Loads cluster_output/results_hypertext.RDS and writes PDFs to
## cluster_output/paper_outputs/ with same filenames as original but no titles.
## =============================================================================

library(hawkesNet)

PKG_ROOT <- if (nzchar(Sys.getenv("SLURM_SUBMIT_DIR"))) {
  Sys.getenv("SLURM_SUBMIT_DIR")
} else {
  getwd()
}
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) PKG_ROOT <- getwd()

source(file.path(PKG_ROOT, "inst", "resolve_output_dir.R"))
OUTPUT_DIR <- hawkesnet_resolve_output_dir(PKG_ROOT)
RDS_PATH <- file.path(OUTPUT_DIR, "results_hypertext.RDS")
PAPER_DIR <- file.path(OUTPUT_DIR, "paper_figures")

if (!file.exists(RDS_PATH)) {
  stop("Results file not found: ", RDS_PATH, "\nRun the hypertext job first.")
}

dir.create(PAPER_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Loading:", RDS_PATH, "\n")
dat <- readRDS(RDS_PATH)
results <- dat$results

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required. Install with: install.packages('ggplot2')")
}

n_saved <- 0L
for (nm in names(results)) {
  res <- results[[nm]]
  if (is.null(res) || is.null(res$gof) || is.null(res$gof$plots)) next

  safe_label <- gsub("[^a-zA-Z0-9_]", "_", tolower(res$label))
  for (plot_name in names(res$gof$plots)) {
    p <- res$gof$plots[[plot_name]]
    if (is.null(p)) next

    ## Remove title for paper
    p_no_title <- p + ggplot2::labs(title = NULL) +
      ggplot2::theme(plot.title = ggplot2::element_blank())

    fname <- file.path(PAPER_DIR, sprintf("gof_%s_%s.pdf", safe_label, plot_name))
    tryCatch({
      ggplot2::ggsave(filename = fname, plot = p_no_title, width = 10, height = 8)
      cat("  Saved:", fname, "\n")
      n_saved <- n_saved + 1L
    }, error = function(e) {
      cat("  Plot save failed:", fname, e$message, "\n")
    })
  }
}

cat("\nDone. Saved", n_saved, "plots to", PAPER_DIR, "\n")
