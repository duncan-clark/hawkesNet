# Rscript hawkesNet/inst/cran/refit_shared_consumer.R OUTPUT LAW [train|full|both|check]
# Run CS-1/check once to build the shared cache before starting parallel fits.
args <- commandArgs(TRUE)
if(length(args)<2L)stop("Supply output directory and law")
out <- args[1]
law <- match.arg(args[2],c("CS-1","CS-pair","CS-pair-migration-add","TERGM-logdegree"))
mode <- if(length(args)<3L)"both" else match.arg(args[3],c("train","full","both","check"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pkgload::load_all("hawkesNet",quiet=TRUE,debug=FALSE)
data_path <- "paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
diagnostics <- "paper_1/6_full_rewrite/diagnostics"
archive <- file.path(diagnostics,"pattern_gof_20261107/fits")
m <- cran_update_data(data_path,law,term_specification="shared_consumer")
m$release_fits <- readRDS(file.path(diagnostics,"exponential_rates_20261105/fits/release_fits.rds"))
m$release_values <- lapply(m$release_fits,"[[","values")
fixed <- m$lower==m$upper
expected <- c("CS-1"=31L,"CS-pair"=37L,"CS-pair-migration-add"=46L,"TERGM-logdegree"=31L)[[law]]
stopifnot(sum(!fixed)+15L==expected,m$pattern_cache$consumer_checks>0L)
windows <- if(mode %in% c("both","check"))c("train","full") else mode
for(window in windows){
  name <- paste0(law,"_",window)
  dest <- file.path(out,paste0(name,".rds"))
  if(file.exists(dest)&&mode!="check"){
    stopifnot(identical(readRDS(dest)$term_specification,"shared_consumer"))
    next
  }
  old <- readRDS(file.path(archive,paste0(if(law=="TERGM-logdegree")"TERGM" else law,"_",window,".rds")))
  start <- setNames(numeric(length(m$names)),m$names)
  same <- intersect(names(old$parameters),m$names)
  start[same] <- old$parameters[same]
  # The zero rates are for the nested check only. Start both fits in the interior.
  start[m$meanConsumer] <- -1000
  replay <- loglik_cran_updates(m,start,TRUE)
  null_error <- max(abs(replay$loglik-old$contributions$edge))
  stopifnot(null_error<1e-6)
  start[m$meanConsumer] <- rep(log(.15),length(m$meanConsumer))
  added <- grep("shared_consumer|log_m_consumer",m$names)
  check <- start
  check[grep("shared_consumer",m$names)] <- .2
  direction <- numeric(length(start))
  direction[added] <- sin(seq_along(added))
  h <- 1e-6
  score <- sum(loglik_cran_updates(m,check)$gradient*direction)
  fd <- sum(loglik_cran_updates(m,check+h*direction,TRUE)$loglik-
            loglik_cran_updates(m,check-h*direction,TRUE)$loglik)/(2*h)
  error <- abs(score-fd)/max(1,abs(score),abs(fd))
  stopifnot(error<1e-5)
  record <- list(law=law,window=window,term_specification="shared_consumer",
                 parameter_count=expected,cache=m$pattern_cache,
                 nested_baseline_error=null_error,gradient_error=error)
  message(name,": nested baseline error ",signif(null_error,3),", gradient error ",signif(error,3))
  jsonlite::write_json(record,file.path(out,paste0(name,"_validation.json")),auto_unbox=TRUE,pretty=TRUE,digits=16)
  if(mode=="check")next
  alternate <- start
  alternate[grep("shared_consumer",m$names)] <- .6
  alternate[m$meanConsumer] <- rep(log(.4),length(m$meanConsumer))
  checkpoint <- file.path(out,paste0(name,"_checkpoint.rds"))
  if(file.exists(checkpoint))start <- readRDS(checkpoint)
  fit <- fit_cran_updates(m,window,parameters=list(start,alternate),checkpoint=checkpoint,
                          control=list(polish_tol=2e-6))
  stopifnot(fit$parameter_count==expected,all(fit$parameters[fixed]==m$lower[fixed]),
            fit$edge_loglik>=old$edge_loglik-1e-6)
  saveRDS(fit,dest)
  jsonlite::write_json(c(record,list(loglik=fit$update_loglik,predictive=fit$predictive,
    parameters=as.list(fit$parameters),bounds=fit$bound_parameters,
    projected_gradient=fit$projected_gradient,multistart_loglik=fit$multistart_loglik)),
    file.path(out,paste0(name,".json")),auto_unbox=TRUE,pretty=TRUE,digits=16)
}
