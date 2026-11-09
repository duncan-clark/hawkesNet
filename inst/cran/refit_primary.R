# Fit the manuscript's motif-free complete-update laws on an existing cache.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: refit_primary.R CRAN_CACHE OUTPUT_DIRECTORY PERSONNEL_RDS")
suppressPackageStartupMessages(library(hawkesNet))
path <- normalizePath(args[1],mustWork=TRUE)
people <- readRDS(args[3])
output <- args[2]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
if(length(list.files(output))) stop("Use an empty output directory to avoid mixing analyses")
for(law in c("PA","CS-1","CS-pair","TERGM-logdegree")) {
  data <- cran_update_data(path, law = law, people = if (law == "PA") NULL else people)
  label <- if(law=="TERGM-logdegree") "TERGM" else law
  for(window in c("train","full")) {
    checkpoint <- file.path(output,paste0(label,"_",window,"_checkpoint.rds"))
    fit <- fit_cran_updates(data,window=window,starts=2L,checkpoint=checkpoint)
    saveRDS(fit,file.path(output,paste0(label,"_",window,".rds")))
    unlink(checkpoint)
  }
  rm(data)
  gc()
}
