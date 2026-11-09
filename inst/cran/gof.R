# Rscript hawkesNet/inst/cran/gof.R DATA_DIRECTORY FIT_DIRECTORY OUTPUT_DIRECTORY
args<-commandArgs(TRUE);if(length(args)<3L)stop('Supply data, fits, output; optionally a model name')
if(dir.exists('hawkesNet'))pkgload::load_all('hawkesNet',quiet=TRUE,debug=FALSE)else library(hawkesNet)
path<-normalizePath(args[1]);fits<-normalizePath(args[2]);output<-args[3];dir.create(output,recursive=TRUE,showWarnings=FALSE)
laws<-if(length(args)>3L)args[4]else c('PA','CS-1','CS-pair','TERGM')
reps<-if(length(args)>4L)as.integer(args[5])else 20L
stopifnot(length(reps)==1L,!is.na(reps),reps>=2L)
portable<-file.path(output,'graph_2025.rds')
if(!file.exists(portable))saveRDS(cran_graph_data(path,'2025-01-01'),portable)
d<-readRDS(portable)
times<-(hawkesNet:::cran_read_npy(file.path(path,'mark_forward_gof_20260925/times.npy'))-1546300800)/86400
checkpoints<-as.numeric(as.Date(c('2025-04-01','2025-07-01','2025-10-01','2026-01-01'))-as.Date('2019-01-01'))
for(law in laws){
 fit<-readRDS(file.path(fits,paste0(law,"_train.rds")))
 if(fit$projected_gradient>2e-5)stop("Unfinished fit: ",law)
 run<-function(i){
  dest<-file.path(output,sprintf("%s_%03d.rds",law,i));if(file.exists(dest))return(dest)
  law_index<-match(law,c("PA","CS-1","CS-1U","CS-pair","TERGM","CS-pair-migration-add","TERGM-logdegree"))
  if(is.na(law_index))stop("No reproducible seed allocation for ",law)
  began<-proc.time()[3];x<-simulate_cran_updates(d$state,fit,times,checkpoints,seed=20261105L+i*97L+law_index*10000L)
  x$term_specification<-fit$term_specification
  saveRDS(x,dest);message(law," replicate ",i," elapsed ",round(proc.time()[3]-began,1)," seconds");dest
 }
 cores<-if(.Platform$OS.type=="windows")1L else as.integer(Sys.getenv("HAWKESNET_GOF_CORES","4"))
 # These workers only call the native graph simulator; no BLAS calls occur after fork.
 result<-parallel::mclapply(seq_len(reps),run,mc.cores=cores,mc.preschedule=FALSE)
 if(any(vapply(result,inherits,FALSE,"try-error")))stop("Failed GOF path: ",law)
}
