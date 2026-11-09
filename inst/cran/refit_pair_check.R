# Additional independent full-window starts for auditing CS-pair local maxima.
# Rscript hawkesNet/inst/cran/refit_pair_check.R DATA_DIRECTORY FIT_DIRECTORY
args<-commandArgs(TRUE);if(length(args)!=2L)stop('Supply data and fit directories')
if(dir.exists('hawkesNet'))pkgload::load_all('hawkesNet',quiet=TRUE,debug=FALSE)else library(hawkesNet)
path<-args[1];output<-args[2];m<-cran_update_data(path,'CS-pair')
m$release_values<-readRDS(file.path(output,'release_values.rds'))
p<-readRDS(file.path(output,'CS-1_full.rds'))$parameters
one<-setNames(numeric(length(m$names)),m$names);one[names(p)]<-p;one[m$meanP]<-log(.2);one<-pmin(m$upper,pmax(m$lower,one))
two<-one;two[which(m$names=='migration')]<-0;two[m$meanP]<-log(.5)
f<-fit_cran_updates(m,'full',parameters=list(one,two),checkpoint=file.path(output,'CS-pair_parallel_full_checkpoint.rds'))
saveRDS(f,file.path(output,'CS-pair_parallel_full.rds'))
jsonlite::write_json(unclass(f),file.path(output,'CS-pair_parallel_full.json'),auto_unbox=TRUE,pretty=TRUE,digits=16,na='null')
