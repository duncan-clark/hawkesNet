# Frozen-state prediction for one hypothetical new package, with no copying.
# This is a diagnostic of fitted weights, not another fitted specification.
# Rscript hawkesNet/inst/cran/diagnose_exponential.R FIT_DIRECTORY GOF_DIRECTORY
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
pkgload::load_all('hawkesNet',quiet=TRUE,debug=FALSE)
d<-readRDS(file.path(args[2],'graph_2025.rds'));g<-hawkesNet:::cran_graph_create(d$state)
result<-list();i<-0L
for(date in c('2025-01-01','2026-01-01')){
 if(date=='2026-01-01')hawkesNet:::cran_graph_replay(g,d$records[vapply(d$records,function(r)r$time>=d$origin,FALSE)])
 state<-hawkesNet:::cran_graph_snapshot(g);degree<-tabulate(state$edges[,2]+1L,nbins=length(state$names));t<-as.numeric(as.Date(date)-as.Date('2019-01-01'))
 for(law in c('PA','CS-1','CS-pair')){
  p<-readRDS(file.path(args[1],paste0(law,'_train.rds')))$parameters
  z<-(if(law=='PA')p['pa_add']else p['popularity_add'])*log1p(degree)-p['age_decay_add']*(t-state$birth)/365.25
  if(law!='PA')z<-z+p['addition_in2star']*degree/100+p['addition_in3star']*degree*(degree-1)/20000
  w<-exp(z-max(z));pr<-w/sum(w);mean<-exp(p['log_m_add']+p['mean_add_new_package'])
  i<-i+1L;result[[i]]<-data.frame(date=date,model=law,proposal_mean=unname(mean),largest_target_probability=max(pr),expected_distinct_additions=sum(-expm1(-mean*pr)))
 }
}
ans<-do.call(rbind,result);utils::write.csv(ans,file.path(args[1],'concentration_diagnostic.csv'),row.names=FALSE);print(ans,row.names=FALSE)
