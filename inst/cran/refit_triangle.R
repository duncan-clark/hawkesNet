# Run from the project root: Rscript hawkesNet/inst/cran/refit_triangle.R OUTPUT [train|full|both]
args <- commandArgs(TRUE);stopifnot(length(args)>=1L)
out <- args[1];dir.create(out,recursive=TRUE,showWarnings=FALSE)
windows <- if(length(args)==1L || args[2]=="both")c("train","full")else match.arg(args[2],c("train","full"))
if(dir.exists("hawkesNet"))pkgload::load_all("hawkesNet",quiet=TRUE,debug=FALSE)else library(hawkesNet)
root <- "paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
archive <- "paper_1/6_full_rewrite/diagnostics/exponential_rates_20261105/fits"
message("Preparing eligible pairs and legacy edit risks")
m <- cran_update_data(root,"CS-pair-triangle")
m$release_values <- readRDS(file.path(archive,"release_values.rds"))
m$release_fits <- readRDS(file.path(archive,"release_fits.rds"))
message(length(m$triangle$group)," observed eligible pairs in ",length(unique(m$triangle$group))," releasing lines; ",
        m$triangle$normaliser_checks," direct checks of full candidate sums")
for(window in windows){
 name <- paste0("CS-pair-triangle_",window)
 dest <- file.path(out,paste0(name,".rds"))
 if(file.exists(dest))next
 old <- readRDS(file.path(archive,paste0("CS-pair_",window,".rds")))
 stopifnot(identical(head(m$names,24),names(old$parameters)))
 fixed <- c(old$parameters,rep(-1e6,3),0)
 nested <- loglik_cran_updates(m,fixed,TRUE)
 error <- max(abs(nested$loglik-old$contributions$edge))
 stopifnot(error<1e-6)
 p1 <- c(old$parameters,log(c(.03,.3,.5)),4)
 p2 <- c(old$parameters,log(c(.15,1,1.2)),8)
 train <- file.path(out,"CS-pair-triangle_train.rds")
 if(window=="full"&&file.exists(train)){p2<-p1;p1<-readRDS(train)$parameters}
 v <- sin(seq_along(p1));h <- 3e-6
 analytic <- sum(loglik_cran_updates(m,p1)$gradient*v)
 numerical <- (loglik_cran_updates(m,p1+h*v)$loglik-loglik_cran_updates(m,p1-h*v)$loglik)/(2*h)
 rel <- abs(analytic-numerical)/max(1,abs(analytic),abs(numerical))
 stopifnot(rel<1e-5)
 message(window,": nesting max error ",format(error,digits=4),"; directional gradient relative error ",format(rel,digits=4))
 fit <- fit_cran_updates(m,window,parameters=list(p1,p2),
             checkpoint=file.path(out,paste0(name,"_checkpoint.rds")))
 saveRDS(fit,dest)
 jsonlite::write_json(list(nested_error=error,relative_gradient_error=rel,normaliser_checks=m$triangle$normaliser_checks,
   observed_eligible_pairs=length(m$triangle$group),parameters=as.list(fit$parameters),
   projected_gradient=fit$projected_gradient,starts=fit$multistart_loglik,bounds=fit$bound_parameters,
   loglik=fit$update_loglik,predictive=fit$predictive),file.path(out,paste0(name,".json")),
   pretty=TRUE,auto_unbox=TRUE,digits=16)
 print(fit$parameters[c(m$meanT,m$copyT)])
}
