# Rscript hawkesNet/inst/cran/refit_pattern_preferences.R OUTPUT LAW [train|full|both|check]
# Build the common pattern cache once with CS-1/check before concurrent fits.
args<-commandArgs(TRUE);if(length(args)<2L)stop("Supply output directory and law")
out<-args[1];law<-match.arg(args[2],c("CS-1","CS-pair","CS-pair-migration-add","TERGM-logdegree"))
mode<-if(length(args)<3L)"both"else match.arg(args[3],c("train","full","both","check"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pkgload::load_all("hawkesNet",quiet=TRUE,debug=FALSE)
root<-"paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
diagnostics<-"paper_1/6_full_rewrite/diagnostics"
original<-file.path(diagnostics,"exponential_rates_20261105/fits")
archive<-switch(law,`CS-pair-migration-add`=file.path(diagnostics,"migration_addition_comparison_20261106"),
                `TERGM-logdegree`=file.path(diagnostics,"tergm_logdegree_20261107"),original)
m<-cran_update_data(root,law,term_specification="pattern_preferences")
legacy<-cran_update_data(root,law,term_specification="legacy")
m$release_fits<-readRDS(file.path(original,"release_fits.rds"));m$release_values<-lapply(m$release_fits,"[[","values")
fixed<-m$lower==m$upper
expected<-c(`CS-1`=30L,`CS-pair`=35L,`CS-pair-migration-add`=41L,`TERGM-logdegree`=30L)[[law]]
stopifnot(sum(!fixed)+15L==expected,all(m$lower[grepl("in2star|in3star|in_[23]star",m$names)]==0),
          all(m$upper[grepl("in2star|in3star|in_[23]star",m$names)]==0),
          all(m$lower[fixed]==m$upper[fixed]))
windows<-if(mode %in% c("both","check"))c("train","full")else mode
for(window in windows){
 name<-paste0(law,"_",window)
 old<-readRDS(file.path(archive,paste0(name,".rds")))
 baseline<-old$parameters
 baseline[grepl("in2star|in3star|in_[23]star|cyclic",names(baseline))]<-0
 if(law %in% c("CS-pair","CS-pair-migration-add"))baseline["migration"]<-0
 start<-setNames(numeric(length(m$names)),m$names)
 same<-intersect(names(baseline),names(start));start[same]<-baseline[same]
 if(!is.null(m$pair_preferences))start[m$pair_preferences]<-baseline[c("popularity_add","popularity_delete","age_decay_add")]
 start<-pmin(m$upper,pmax(m$lower,start))
 reference<-loglik_cran_updates(legacy,baseline,TRUE)
 evaluation<-loglik_cran_updates(m,start,TRUE)
 null_error<-max(abs(reference$loglik-evaluation$loglik))
 stopifnot(null_error<1e-6)
 # Check the new diamond coordinate away from zero, along with all free terms.
 check<-start;check[match("addition_diamond",m$names)]<-.002
 if(!is.null(m$ip))check[match("pair_diamond",m$names)]<--.003
 direction<-sin(seq_along(check)) * (!fixed) / pmax(1,abs(check))
 h<-1e-7
 score<-sum(loglik_cran_updates(m,check)$gradient*direction)
 fd<-sum(loglik_cran_updates(m,check+h*direction,TRUE)$loglik-
         loglik_cran_updates(m,check-h*direction,TRUE)$loglik)/(2*h)
 error<-abs(score-fd)/max(1,abs(score),abs(fd))
 message(name,": ",expected," parameters; nested baseline error ",signif(null_error,3),"; gradient error ",signif(error,3))
 stopifnot(error<1e-5)
 record<-list(law=law,window=window,term_specification="pattern_preferences",parameter_count=expected,
   fixed_parameters=m$names[fixed],single_preferences=c("migration","transitivity","focal open diamond"),
   independent_pair_parameters=if(is.null(m$pair_preferences))character()else m$names[c(m$ip[3:4],m$pair_preferences)],
   cache=m$pattern_cache,nested_baseline_error=null_error,gradient_error=error,status="validated")
 jsonlite::write_json(record,file.path(out,paste0(name,"_validation.json")),auto_unbox=TRUE,pretty=TRUE,digits=16)
 if(mode=="check")next
 dest<-file.path(out,paste0(name,".rds"))
 if(file.exists(dest)){fit<-readRDS(dest);stopifnot(identical(fit$term_specification,"pattern_preferences"));next}
 checkpoint<-file.path(out,paste0(name,"_checkpoint.rds"))
 alternative<-start;alternative[match("addition_diamond",m$names)]<-.01
 if(!is.null(m$ip))alternative[match("pair_diamond",m$names)]<-.01
 if(file.exists(checkpoint))start<-readRDS(checkpoint)
 fit<-fit_cran_updates(m,window,parameters=list(start,alternative),checkpoint=checkpoint,control=list(polish_tol=2e-6))
 stopifnot(fit$parameter_count==expected,all(fit$parameters[fixed]==m$lower[fixed]))
 saveRDS(fit,dest);record$status<-"fitted"
 jsonlite::write_json(c(record,list(loglik=fit$update_loglik,predictive=fit$predictive,
   parameters=as.list(fit$parameters),bounds=fit$bound_parameters,projected_gradient=fit$projected_gradient,
   multistart_loglik=fit$multistart_loglik)),file.path(out,paste0(name,".json")),auto_unbox=TRUE,pretty=TRUE,digits=16)
}
