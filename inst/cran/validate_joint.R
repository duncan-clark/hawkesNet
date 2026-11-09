# Verify simulation and likelihood event rates at every forecast origin.
if(dir.exists("hawkesNet"))pkgload::load_all("hawkesNet",quiet=TRUE,debug=FALSE)else library(hawkesNet)
args<-commandArgs(TRUE);if(length(args)<5L)stop("Supply DATA, FITS, JOINT, TIMING_JSON and CLIPPED_M2_JSON")
path<-args[1];out<-args[2]
origins<-readRDS(file.path(args[3],"origins.rds"))
archive<-jsonlite::fromJSON(args[4],simplifyVector=FALSE)$fits
clipped<-jsonlite::fromJSON(args[5],simplifyVector=FALSE)$results$two_exp_M2_clipped
fit<-readRDS(file.path(out,'PA_train.rds'));d<-cran_timing_data(path);result<-list()
for(model in c('M0','M1','M2')){
 if(model=='M2'){q<-unlist(clipped$parameters);p<-c(q[1:3],q[7:8],q[4:6])}else{a<-archive[[model]];q<-setNames(unlist(a$parameters),unlist(a$parameter_names));p<-c(q[c('log_mu','log_alpha_short','log_alpha_long','log_tau_short','log_tau_long')],if(model=='M1')q[c('gamma_birth','delta')])}
 amp<-switch(model,M0=c(0,0,0),M1=c(p[7],p[6],0),M2=c(p[8],p[6:7]));tm<-setNames(c(exp(p[1:5]),amp),c('mu','alpha_short','alpha_long','tau_short','tau_long','delta','gamma_B','gamma_Delta'))
 rates<-vapply(origins,function(o){a<-simulate_cran_joint(o$state,o$history,fit,tm,d$calendar,o$origin,horizon=0,model=model,seed=1)$initial_intensity;b<-loglik_cran_timing(d,p,model,end=o$origin)$end_intensity;c(simulation=a,likelihood=b)},numeric(2))
 result[[model]]<-rates;cat(model,'maximum absolute discrepancy',max(abs(rates[1,]-rates[2,])),'\n');flush.console();stopifnot(max(abs(rates[1,]-rates[2,]))<1e-6)
}
saveRDS(result,file.path(out,'joint_rate_validation.rds'))
