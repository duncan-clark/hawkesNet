# Rscript hawkesNet/inst/cran/refit_diamond.R OUTPUT [train|full|both] [WARM_FITS|add-only|migration-add]
args <- commandArgs(TRUE);stopifnot(length(args)>=1L)
out <- args[1];dir.create(out,recursive=TRUE,showWarnings=FALSE)
windows <- if(length(args)==1L||args[2]=="both")c("train","full")else match.arg(args[2],c("train","full"))
if(dir.exists("hawkesNet"))pkgload::load_all("hawkesNet",quiet=TRUE,debug=FALSE)else library(hawkesNet)
root <- "paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
diagnostics <- "paper_1/6_full_rewrite/diagnostics"
original <- file.path(diagnostics,"exponential_rates_20261105","fits")
triangle <- file.path(diagnostics,"triangle_channel_20261106")
mode <- if(length(args)>2L)args[3] else ""
law <- switch(mode,`add-only`="CS-pair-add",`migration-add`="CS-pair-migration-add","CS-pair-diamond")
message("Preparing exact paired-addition candidate masses for ",law)
m <- cran_update_data(root,law)
m$release_values <- readRDS(file.path(original,"release_values.rds"))
m$release_fits <- readRDS(file.path(original,"release_fits.rds"))
stopifnot(length(m$triangle$group)==58008L,length(m$triangle$diamond$group)==72757L)
message(length(m$triangle$diamond$group)," observed open-diamond pairs; ",m$triangle$normaliser_checks," direct normaliser checks")
for(window in windows){
  name <- paste0(law,"_",window)
  dest <- file.path(out,paste0(name,".rds"))
  if(file.exists(dest))next
  if(law!="CS-pair-diamond"){
    base <- if(law=="CS-pair-add")"CS-1" else "CS-pair"
    old <- readRDS(file.path(original,paste0(base,"_",window,".rds")))
    combined <- readRDS(file.path(diagnostics,"diamond_channel_distinct_20261106",paste0("CS-pair-diamond_",window,".rds")))
    stopifnot(identical(m$names,names(combined$parameters)),all(names(old$parameters) %in% m$names))
    p1 <- pmin(m$upper,pmax(m$lower,combined$parameters))
    single <- match(names(old$parameters),m$names)
    p1[match("migration",m$names)] <- old$parameters[["migration"]]
    p2 <- p1
    p2[single] <- old$parameters
    disabled <- p2
    disabled[c(m$meanT,m$meanDiamond)] <- -1e6
  }else{
    old <- readRDS(file.path(triangle,paste0("CS-pair-triangle_",window,".rds")))
    stopifnot(identical(head(m$names,28),names(old$parameters)))
    disabled <- c(old$parameters,rep(-1e6,3))
    p1 <- c(old$parameters,log(c(.03,.15,.10)))
    p2 <- c(old$parameters,log(c(.18,.60,.35)))
    train <- file.path(out,"CS-pair-diamond_train.rds")
    if(window=="full"&&file.exists(train))p2 <- readRDS(train)$parameters
    if(length(args)>2L){
      warm <- file.path(args[3],paste0(name,".rds"))
      if(file.exists(warm)){
        p1<-p2<-readRDS(warm)$parameters
        p2[m$meanDiamond]<-p2[m$meanDiamond]+c(.3,-.3,.5)
      }
    }
  }
  nested <- loglik_cran_updates(m,disabled,TRUE)
  error <- max(abs(nested$loglik-old$contributions$edge));stopifnot(error<1e-6)
  check <- p1
  v <- sin(seq_along(check))
  h <- if(law=="CS-pair-migration-add")3e-7 else 3e-6
  g <- sum(loglik_cran_updates(m,check)$gradient*v)
  fd <- sum(loglik_cran_updates(m,check+h*v,TRUE)$loglik-loglik_cran_updates(m,check-h*v,TRUE)$loglik)/(2*h)
  rel <- abs(g-fd)/max(1,abs(g),abs(fd))
  message(window,": nesting error ",signif(error,4),"; directional gradient error ",signif(rel,4))
  stopifnot(rel<1e-5)
  fit <- fit_cran_updates(m,window,parameters=list(p1,p2),control=list(polish_tol=2e-6),
                          checkpoint=file.path(out,paste0(name,"_checkpoint.rds")))
  saveRDS(fit,dest)
  jsonlite::write_json(list(nested_error=error,relative_gradient_error=rel,normaliser_checks=m$triangle$normaliser_checks,
    observed_diamond_pairs=length(m$triangle$diamond$group),parameters=as.list(fit$parameters),
    parameter_count=fit$parameter_count,projected_gradient=fit$projected_gradient,starts=fit$multistart_loglik,bounds=fit$bound_parameters,
    loglik=fit$update_loglik,predictive=fit$predictive),file.path(out,paste0(name,".json")),
    pretty=TRUE,auto_unbox=TRUE,digits=16)
  print(fit$parameters[m$meanDiamond])
}
