# Independent full-candidate and gradient checks for exponential CS refits.
# Rscript hawkesNet/inst/cran/validate_exponential.R DATA FIT_DIRECTORY
args<-commandArgs(TRUE);stopifnot(length(args)>=2L)
pkgload::load_all('hawkesNet',quiet=TRUE,debug=FALSE)
m<-cran_update_data(args[1],'CS-pair');fit<-readRDS(file.path(args[2],'CS-pair_train.rds'));p<-fit$parameters
v<-loglik_cran_updates(m,p);free<-which(m$lower<m$upper)
fd<-vapply(free,function(j){h<-1e-5*max(1,abs(p[j]));d<-numeric(length(p));d[j]<-h;(loglik_cran_updates(m,p+d)$loglik-loglik_cran_updates(m,p-d)$loglik)/(2*h)},0)
gradient_error<-max(abs(fd-v$gradient[free])/pmax(1,abs(fd),abs(v$gradient[free])))
stopifnot(gradient_error<2e-5)
normal<-hawkesNet:::.cran_add_sum(m,p[m$ia],p[m$ageA])
d<-cran_graph_data(args[1],'2020-01-01');g<-hawkesNet:::cran_graph_create(d$state)
records<-d$records[vapply(d$records,function(r)r$time>=d$origin,FALSE)]
groups<-jsonlite::fromJSON(file.path(args[1],'cs2_paired_20260924/groups.json'))
errors<-numeric();checks<-list();gi<-1L
for(b in 1:4){
 r<-records[[b]];birth<-r;birth$add<-birth$delete<-matrix(character(),0,2)
 hawkesNet:::cran_graph_replay(g,list(birth));state<-hawkesNet:::cran_graph_snapshot(g)
 n<-length(state$names);family<-sub('@[^@]+$','',state$names);major<-as.integer(sub('.*@','',state$names))
 ed<-state$edges+1L;out<-split(ed[,2],factor(ed[,1],levels=1:n));inn<-split(ed[,1],factor(ed[,2],levels=1:n));din<-lengths(inn)
 while(gi<=m$ng&&m$gb[gi]==b-1L){
  u<-match(groups$source[gi],state$names);tr<-cy<-numeric(n)
  bump<-function(x,i){if(length(i))x<-x+tabulate(i,nbins=n);x}
  for(w in out[[u]]){tr<-bump(tr,inn[[w]]);tr<-bump(tr,out[[w]])}
  for(w in inn[[u]]){tr<-bump(tr,out[[w]]);cy<-bump(cy,inn[[w]])}
  migration<-vapply(seq_len(n),function(v)any(family[out[[u]]]==family[v]&major[out[[u]]]<major[v]),FALSE)
  predecessor<-which(family==family[u]&major<major[u]);copy<-rep(FALSE,n)
  if(!length(out[[u]])&&length(predecessor)){pr<-predecessor[which.max(major[predecessor])];copy[out[[pr]]]<-TRUE}
  x<-hawkesNet:::.cran_features(din,length(out[[u]]),tr,cy,flags=cbind(migration,copy),exponential=TRUE)
  keep<-setdiff(seq_len(n),c(u,out[[u]]));bt<-(state$birth-365)/365.25
  z<-drop(x[keep,,drop=FALSE]%*%p[m$ia])+p[m$ageA]*bt[keep]
  mx<-max(z);w<-exp(z-mx);s<-sum(w);direct<-mx+log(s)
  score<-colSums(x[keep,,drop=FALSE]*w)/s;age_score<-sum(w*bt[keep])/s
  error<-max(abs(c(direct-normal$logsum[gi],score-normal$score[gi,],age_score-normal$birth[gi])))
  errors<-c(errors,error);checks[[length(checks)+1L]]<-list(group=gi,source=state$names[u],candidates=length(keep),max_error=error)
  gi<-gi+1L
 }
 r$births<-character();r$clocks<-numeric();hawkesNet:::cran_graph_replay(g,list(r))
}
stopifnot(max(errors)<1e-7)
single<-cran_update_data(args[1],'CS-1');sp<-readRDS(file.path(args[2],'CS-1_train.rds'))$parameters
embed<-setNames(numeric(length(m$names)),m$names);embed[names(sp)]<-sp;embed[m$meanP]<--1000
nesting_error<-loglik_cran_updates(m,embed)$loglik-loglik_cran_updates(single,sp)$loglik
constant<-cran_update_data(args[1],'CS-pair',proposal_means='constant');q<-p;q[tail(seq_along(q),4)]<-0
constant_error<-loglik_cran_updates(m,q)$loglik-loglik_cran_updates(constant,p[seq_along(constant$names)])$loglik
stopifnot(abs(nesting_error)<1e-7,abs(constant_error)<1e-7)
legacy_errors<-NULL
if(length(args)>2L){legacy_errors<-vapply(c('PA','CS-1','CS-pair'),function(law){d<-cran_update_data(args[1],law,weights='logistic',proposal_means='constant');f<-readRDS(file.path(args[3],paste0(law,'_full.rds')));loglik_cran_updates(d,f$parameters)$loglik-f$edge_loglik},0);stopifnot(max(abs(legacy_errors))<1e-7)}
result<-list(gradient_max_relative_error=gradient_error,direct_candidate_max_error=max(errors),pair_zero_rate_nesting_error=nesting_error,constant_rate_nesting_error=constant_error,legacy_errors=legacy_errors,contexts=checks,release_class_counts=table(m$release_class),rate_design_rank=c(addition=qr(m$ZA)$rank,deletion=qr(m$ZD[m$exposure[,2]>0,,drop=FALSE])$rank))
saveRDS(result,file.path(args[2],'exponential_validation.rds'));print(result)
