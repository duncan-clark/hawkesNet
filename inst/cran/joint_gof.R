# Rscript hawkesNet/inst/cran/joint_gof.R DATA FITS OUTPUT TIMING_JSON CLIPPED_M2_JSON [SPEC]
args<-commandArgs(TRUE);if(length(args)<5L)stop('Supply data, fits, output, timing fit JSON and clipped M2 JSON')
if(dir.exists('hawkesNet'))pkgload::load_all('hawkesNet',quiet=TRUE,debug=FALSE)else library(hawkesNet)
path<-normalizePath(args[1]);fits<-normalizePath(args[2]);output<-args[3];dir.create(output,recursive=TRUE,showWarnings=FALSE)
archive<-jsonlite::fromJSON(args[4],simplifyVector=FALSE)$fits
clipped<-jsonlite::fromJSON(args[5],simplifyVector=FALSE)$results$two_exp_M2_clipped
physical<-function(model){
 if(model=='M2'){
  p<-unlist(clipped$parameters);# dictionary_engine ordering: baseline, masses, gamma_B, gamma_Delta, delta, log means
  return(setNames(c(exp(p[1:3]),exp(p[7:8]),p[6],p[4:5]),c('mu','alpha_short','alpha_long','tau_short','tau_long','delta','gamma_B','gamma_Delta')))
 }
 x<-archive[[model]];p<-setNames(unlist(x$parameters),unlist(x$parameter_names))
 setNames(c(exp(p[c('log_mu','log_alpha_short','log_alpha_long','log_tau_short','log_tau_long')]),if(model=='M0')c(0,0,0)else c(p['delta'],p['gamma_birth'],0)),c('mu','alpha_short','alpha_long','tau_short','tau_long','delta','gamma_B','gamma_Delta'))
}
calendar<-utils::read.csv(file.path(path,'current_fits/calendar_values.csv'))$value
specs<-data.frame(label=c('M0_CS-pair','M1_CS-pair','M2_CS-pair','M2_PA','M0_TERGM'),law=c(rep('CS-pair',3),'PA','TERGM'),timing=c('M0','M1','M2','M2','M0'))
origin_file<-file.path(output,'origins.rds')
if(!file.exists(origin_file)){
 d<-cran_graph_data(path,'2025-01-01');g<-hawkesNet:::cran_graph_create(d$state);prior<-d$records[vapply(d$records,function(r)r$time<d$origin,FALSE)]
 # Rewind the observed pre-origin changes, then replay them to recover every
 # post-update degree. The pre-origin labels remain isolated until their birth.
 reverse<-lapply(rev(prior),function(r){a<-r$add;r$add<-r$delete;r$delete<-a;r$births<-character();r$clocks<-numeric();r$sources<-character();r})
 hawkesNet:::cran_graph_replay(g,reverse)
 history<-hawkesNet:::cran_graph_replay(g,prior)
 stopifnot(identical(hawkesNet:::cran_graph_snapshot(g)$edges,d$state$edges))
 origins<-as.numeric(as.Date(unlist(lapply(1:12,function(m)sprintf('2025-%02d-%02d',m,c(1,15)))))-as.Date('2019-01-01'))
 records<-d$records[vapply(d$records,function(r)r$time>=d$origin,FALSE)];ci<-1L;result<-vector('list',24)
 for(i in seq_along(origins)){
  take<-which(vapply(records,function(r)r$time<origins[i],FALSE));take<-take[take>=ci]
  if(length(take)){history<-c(history,hawkesNet:::cran_graph_replay(g,records[take]));ci<-max(take)+1L}
  result[[i]]<-list(origin=origins[i],state=hawkesNet:::cran_graph_snapshot(g),history=history)
 };saveRDS(result,origin_file)
}
origins<-readRDS(origin_file)
if(length(args)>5L)specs<-specs[specs$label==args[6],,drop=FALSE]
for(j in seq_len(nrow(specs))){s<-specs[j,];fit<-readRDS(file.path(fits,paste0(s$law,'_train.rds')));tm<-physical(s$timing)
 tasks<-expand.grid(rep=seq_len(32L),origin=seq_along(origins))
 run<-function(task){i<-tasks$origin[task];rep<-tasks$rep[task]
  dest<-file.path(output,sprintf('%s_%02d_%02d.rds',s$label,i,rep));if(file.exists(dest))return(dest)
  o<-origins[[i]];seed<-20261105L+100000L*match(s$label,c('M0_CS-pair','M1_CS-pair','M2_CS-pair','M2_PA','M0_TERGM'))+100L*i+rep
  x<-simulate_cran_joint(o$state,o$history,fit,tm,calendar,o$origin,model=s$timing,seed=seed);saveRDS(x,dest)
  if(rep==32L)message(s$label,' origin ',i,' complete')
 }
 cores<-if(.Platform$OS.type=='windows')1L else as.integer(Sys.getenv('HAWKESNET_JOINT_CORES','4'))
 # Each child calls only the native simulator, with an explicit path seed.
 result<-parallel::mclapply(seq_len(nrow(tasks)),run,mc.cores=cores,mc.preschedule=FALSE)
 if(any(vapply(result,inherits,FALSE,'try-error')))stop('Failed joint path: ',s$label)
}
