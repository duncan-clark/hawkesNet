# Rscript hawkesNet/inst/cran/refit_channels.R OUTPUT LAW [train|full|both] [base|degree] [poisson|single]
args <- commandArgs(TRUE)
stopifnot(length(args)>=2L)
out <- args[1]; law <- match.arg(args[2],c("CS-1","CS-pair-migration-add"))
mode <- if(length(args)>2L)args[3] else "both"
degree <- length(args)>3L && args[4]=="degree"
pair_proposals <- if(length(args)>4L)args[5]else NULL
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pkgload::load_all("hawkesNet",quiet=TRUE,debug=FALSE)
diagnostics <- "paper_1/6_full_rewrite/diagnostics"
path <- "paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
m <- cran_update_data(path,law,term_specification="pair_channels",migration_degree=degree,
                     pair_proposals=pair_proposals)
# Fail before optimization or expensive normalizers if exact summation exceeds
# the declared computational budget. This is never a reason to drop updates.
if(identical(m$pair_proposals,'poisson'))
  hawkesNet:::cran_pair_union_check(m$pair_union_plan,m$pair_max_states)
m$release_fits <- readRDS(file.path(diagnostics,"exponential_rates_20261105/fits/release_fits.rds"))
m$release_values <- lapply(m$release_fits,"[[","values")
for(window in if(mode=="both")c("train","full")else mode){
  dest <- file.path(out,paste0(law,"_",window,".rds"))
  if(file.exists(dest))next
  old <- readRDS(file.path(diagnostics,"shared_consumer_20261107/fits",paste0(law,"_",window,".rds")))
  start <- setNames(numeric(length(m$names)),m$names)
  start[names(old$parameters)] <- old$parameters
  if(degree){
    base <- file.path(diagnostics,"cs_channels_20261108/fits",paste0(law,"_",window,".rds"))
    if(file.exists(base)){p <- readRDS(base)$parameters;start[names(p)] <- p}
  }
  start <- pmin(m$upper,pmax(m$lower,start))
  free <- which(m$lower<m$upper)
  direction <- numeric(length(start));direction[free]<-sin(free)
  h <- 1e-6
  analytic <- sum(loglik_cran_updates(m,start)$gradient*direction)
  numeric <- (loglik_cran_updates(m,start+h*direction)$loglik-loglik_cran_updates(m,start-h*direction)$loglik)/(2*h)
  error <- abs(analytic-numeric)/max(1,abs(analytic),abs(numeric))
  stopifnot(error<1e-5)
  fit <- fit_cran_updates(m,window,parameters=start,
    checkpoint=sub(".rds$","_checkpoint.rds",dest),control=list(polish_tol=2e-6))
  stopifnot(fit$parameter_count==(if(law=="CS-1")31L else 42L)+as.integer(degree))
  saveRDS(fit,dest)
  jsonlite::write_json(list(law=law,window=window,term_specification=fit$term_specification,
    pair_proposals=fit$pair_proposals,pair_max_states=fit$pair_max_states,
    migration_degree=degree,parameter_count=fit$parameter_count,loglik=fit$update_loglik,predictive=fit$predictive,
    parameters=as.list(fit$parameters),fixed=fit$fixed_parameters,bounds=fit$bound_parameters,
    gradient_check=error,projected_gradient=fit$projected_gradient,multistart=fit$multistart_loglik),
    sub(".rds$",".json",dest),auto_unbox=TRUE,pretty=TRUE,digits=16)
}
