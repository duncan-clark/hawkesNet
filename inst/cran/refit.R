# Run from the project root after installing hawkesNet:
# Rscript hawkesNet/inst/cran/refit.R DATA_DIRECTORY OUTPUT_DIRECTORY
args<-commandArgs(TRUE)
if(length(args)<2L)stop('Supply the CRAN cache directory and output directory; optionally a model name')
if(dir.exists('hawkesNet'))pkgload::load_all('hawkesNet',quiet=TRUE,debug=FALSE)else library(hawkesNet)
data_root<-normalizePath(args[1]);output<-args[2];dir.create(output,recursive=TRUE,showWarnings=FALSE)
laws<-if(length(args)>2L)args[3] else c('PA','CS-1','CS-pair','TERGM')
# All laws share the same release model. Evaluate it once in each estimation window.
release_cache<-file.path(output,'release_values.rds')
if(!file.exists(release_cache)){
 m<-cran_update_data(data_root,'PA');rel<-list()
 for(window in c('train','full')){
  p<-hawkesNet:::.cran_release_parameters(m,window)
  rel[[window]]<-hawkesNet:::cran_release_evaluate(m$release_data,p,m$families)
 };saveRDS(rel,release_cache);rm(m);gc()
}
rel<-readRDS(release_cache)
release_fit_file<-file.path(output,"release_fits.rds")
if(!file.exists(release_fit_file)){m<-cran_update_data(data_root,"PA");m$release_values<-rel;rf<-setNames(lapply(c("train","full"),function(w)fit_cran_release(m,w)),c("train","full"));saveRDS(rf,release_fit_file)}
release_fits<-readRDS(release_fit_file)
rel<-lapply(release_fits,`[[`,"values")
run<-function(law){
 m<-cran_update_data(data_root,law);m$release_values<-rel;m$release_fits<-release_fits
 record<-list();p0<-hawkesNet:::.cran_old_start(m,'full')
 if(law!='TERGM'){
  nested<-p0;if(law=='CS-pair')nested[m$meanP]<--1e6
  v<-loglik_cran_updates(m,nested,TRUE);old<-hawkesNet:::.cran_old_fit(m,'full')
  error<-sum(v$loglik)+sum(rel$full$loglik)-old$summary$mark_loglik
  if(m$proposal_means=='constant'&&!m$exponential&&abs(error)>1e-6)stop('Shared-coefficient nesting failed: ',law,' ',error)
  v<-loglik_cran_updates(m,p0,TRUE)
  direction<-sin(seq_along(p0));h<-1e-6
  fd<-(loglik_cran_updates(m,p0+h*direction)$loglik-loglik_cran_updates(m,p0-h*direction)$loglik)/(2*h)
  analytic<-sum(colSums(v$gradient)*direction)
  if(abs(fd-analytic)/max(1,abs(fd),abs(analytic))>1e-5)stop('Directional gradient check failed: ',law)
  record$validation<-list(directional_gradient_error=fd-analytic)
 }
 for(window in c('train','full')){
  dest<-file.path(output,paste0(law,'_',window,'.rds'))
  if(file.exists(dest)){previous<-readRDS(dest);if(previous$projected_gradient<=2e-5){record[[window]]<-previous;next}}
  start<-if(window=='full')record$train$parameters else NULL
  if(law=='CS-1'&&window=='train'){
   pa<-readRDS(file.path(output,'PA_train.rds'))$parameters
   start<-setNames(numeric(length(m$names)),m$names);common<-intersect(names(pa),m$names);start[common]<-pa[common]
   start['popularity_add']<-pa['pa_add'];start['popularity_delete']<--pa['pa_delete']
  }
  if(law=='CS-pair'){
   single<-readRDS(file.path(output,paste0('CS-1_',window,'.rds')))$parameters
   embed<-setNames(numeric(length(m$names)),m$names);embed[names(single)]<-single;embed[m$meanP]<-log(.2);embed<-pmin(m$upper,pmax(m$lower,embed))
   alternative<-embed;alternative[which(m$names=='migration')]<-0;alternative[m$meanP]<-log(.5)
   start<-if(window=='train')list(embed,alternative)else list(start,alternative)
  }
  if(exists('previous')&&previous$window==window&&previous$projected_gradient>2e-5)start<-previous$parameters
  full_path<-file.path(output,paste0(law,"_full.rds"))
  if(window=="train"&&!is.null(start)&&file.exists(full_path))start<-list(start,readRDS(full_path)$parameters)
  checkpoint<-file.path(output,paste0(law,'_',window,'_checkpoint.rds'))
  if(is.null(start)&&file.exists(checkpoint))start<-readRDS(checkpoint)
  # A previously fitted model is a valid starting value, never a replacement result.
  fit<-fit_cran_updates(m,window,parameters=start,checkpoint=checkpoint)
  saveRDS(fit,dest);jsonlite::write_json(unclass(fit),file.path(output,paste0(law,'_',window,'.json')),auto_unbox=TRUE,pretty=TRUE,digits=16,na='null')
  record[[window]]<-fit
 };saveRDS(record,file.path(output,paste0(law,'_record.rds')));law
}
# Independent models can be fitted in separate R processes; each native sum uses four workers.
results<-lapply(laws,run)
if(any(vapply(results,inherits,FALSE,'try-error')))stop('One or more model fits failed')
manifest<-list(runtime=R.version.string,package_version=as.character(utils::packageVersion('hawkesNet')),data=data_root,laws=laws,training='2020--2023',full='2020--2025',prediction=c(2024,2025),weights='exponential CS proposal weights; directed single and joint pair motifs',proposal_means='addition release class and log1p outdegree; deletion log1p outdegree',age='separate continuous age decay by operation',penalty=FALSE,completed=as.character(Sys.time()))
jsonlite::write_json(manifest,file.path(output,paste0(laws[1],'-manifest.json')),auto_unbox=TRUE,pretty=TRUE)
