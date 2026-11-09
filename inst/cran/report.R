# Rscript hawkesNet/inst/cran/report.R DATA FITS GOF JOINT FIGURES
args<-commandArgs(TRUE);if(length(args)<5L)stop('Supply data, fits, mark GOF, joint GOF and figure directory')
if(dir.exists('hawkesNet'))pkgload::load_all('hawkesNet',quiet=TRUE,debug=FALSE)else library(hawkesNet)
library(ggplot2)
path<-normalizePath(args[1]);fits<-args[2];gof<-args[3];joint<-args[4];figures<-args[5];dir.create(figures,recursive=TRUE,showWarnings=FALSE)
laws<-if(length(args)>6L)strsplit(args[7],',',fixed=TRUE)[[1]]else c('PA','CS-1','CS-pair','TERGM')
cols<-c(PA='#0072B2',`CS-1`='#E69F00',`CS-pair`='#D55E00',`CS-pair-migration-add`='#009E73',TERGM='#CC79A7')
model_labels<-c(PA='PA',`CS-1`='CS-1',`CS-pair`='CS-pair-mig.',`CS-pair-migration-add`='CS-pair-mig.+add.',TERGM='TERGM')
main_pair<-if('CS-pair-migration-add' %in% laws)'CS-pair-migration-add'else 'CS-pair'
mode<-if(length(args)>5L)args[6]else'all';stopifnot(mode %in% c('all','marks','joint'))
base<-theme_bw(base_size=10)+theme(legend.position='bottom',panel.grid.minor=element_blank())
saveplot<-function(p,name,w,h)ggsave(file.path(figures,paste0(name,'.pdf')),p,width=w,height=h,device='pdf')
fitrows<-lapply(laws,function(law){train<-readRDS(file.path(fits,paste0(law,'_train.rds')));full<-readRDS(file.path(fits,paste0(law,'_full.rds')));if(max(train$projected_gradient,full$projected_gradient)>2e-5)stop('Unfinished fit: ',law);data.frame(model=law,parameters=full$parameter_count,full_loglik=full$update_loglik,predictive_2024=train$predictive[['2024']]$update_loglik,predictive_2025=train$predictive[['2025']]$update_loglik)})
utils::write.csv(do.call(rbind,fitrows),file.path(fits,'comparison.csv'),row.names=FALSE)
if(mode %in% c('all','marks')){
d<-readRDS(file.path(gof,'graph_2025.rds'));g<-hawkesNet:::cran_graph_create(d$state);checkpoints<-as.numeric(as.Date(c('2025-04-01','2025-07-01','2025-10-01','2026-01-01'))-as.Date('2019-01-01'));observed<-vector('list',4);ci<-1L
for(k in seq_along(checkpoints)){
 use<-which(vapply(d$records,function(r)r$time>=d$origin&&r$time<checkpoints[k],FALSE));use<-use[use>=ci]
 if(length(use)){hawkesNet:::cran_graph_replay(g,d$records[use]);ci<-max(use)+1L}
 observed[[k]]<-hawkesNet:::cran_graph_summary(g)
}
obs_events<-jsonlite::fromJSON(file.path(path,'mark_forward_gof_20260925/observed.json'))$events
cdf<-function(h,x){ix<-as.numeric(names(h));vapply(x,function(z)sum(h[ix<=z])/sum(h),0)}
band<-function(x)stats::quantile(x,c(.025,.5,.975),names=FALSE)
distributions<-trajectories<-events<-distances<-list();di<-ti<-ei<-ki<-0L
for(law in laws){
 files<-file.path(gof,sprintf('%s_%03d.rds',law,1:20));if(!all(file.exists(files)))stop('Incomplete mark simulations: ',law);paths<-lapply(files,readRDS)
 for(q in 1:4){for(metric in names(observed[[q]]$metrics)){
  if(metric=='betweenness_landmarks')next
  values<-vapply(paths,function(x)x$snapshots[[q]]$metrics[[metric]],0);b<-band(values);ti<-ti+1L;trajectories[[ti]]<-data.frame(model=law,date=as.Date('2019-01-01')+checkpoints[q],metric=metric,observed=observed[[q]]$metrics[[metric]],lo=b[1],median=b[2],hi=b[3])
 }}
 for(metric in c('indegree','outdegree','esp')){
  h<-observed[[4]]$distributions[[metric]];thresholds<-if(metric=='esp')c(0:10,15,20,30,50,100)else c(0:10,15,20,30,50,75,100,150,200,300,500,1000,2000,4000,8000)
  actual<-cdf(h,thresholds);mat<-sapply(paths,function(x)cdf(x$snapshots[[4]]$distributions[[metric]],thresholds));bs<-t(apply(mat,1,band));di<-di+1L;distributions[[di]]<-data.frame(model=law,metric=metric,threshold=thresholds,observed=actual,lo=bs[,1],median=bs[,2],hi=bs[,3])
  dd<-vapply(paths,function(x){s<-x$snapshots[[4]]$distributions[[metric]];z<-sort(unique(c(as.numeric(names(h)),as.numeric(names(s)))));max(abs(cdf(s,z)-cdf(h,z)))},0);ki<-ki+1L;distances[[ki]]<-data.frame(model=law,metric=metric,median_distance=median(dd),lo=band(dd)[1],hi=band(dd)[3])
 }
 for(metric in c('births','additions','deletions','releasing_lines')){
  key<-c(births='nodes_added',additions='edges_added',deletions='edges_deleted',releasing_lines='actors')[[metric]];v<-vapply(paths,function(x)mean(x$sizes[,metric]),0);b<-band(v);ei<-ei+1L;events[[ei]]<-data.frame(model=law,metric=metric,observed=mean(obs_events[[key]]),lo=b[1],median=b[2],hi=b[3])
 }
}
D<-do.call(rbind,distributions);T<-do.call(rbind,trajectories);E<-do.call(rbind,events);K<-do.call(rbind,distances)
for(n in c('D','T','E','K'))utils::write.csv(get(n),file.path(gof,paste0(n,'.csv')),row.names=FALSE)
D$model<-factor(D$model,laws);T$model<-factor(T$model,laws);E$model<-factor(E$model,laws)
D$metric<-factor(D$metric,c('indegree','outdegree','esp'))
distribution_labels<-as_labeller(c(indegree='In-degree',outdegree='Out-degree',esp='Shared partners'))
count_scale<-function()scale_x_continuous(trans='log1p',breaks=c(0,1,10,100,1000),labels=function(x)format(x,trim=TRUE,scientific=FALSE))
base<-theme_bw(base_size=10)+theme(legend.position='bottom',panel.grid.minor=element_blank())
saveplot<-function(p,name,w,h)ggsave(file.path(figures,paste0(name,'.pdf')),p,width=w,height=h,device='pdf')
p<-ggplot(subset(D,model %in% c('PA',main_pair,'TERGM')),aes(threshold,median-observed,color=model,fill=model))+geom_hline(yintercept=0,color='grey50')+geom_ribbon(aes(ymin=lo-observed,ymax=hi-observed),alpha=.2,color=NA)+geom_line()+facet_grid(metric~model,scales='free_y',labeller=labeller(metric=distribution_labels,model=as_labeller(model_labels)))+count_scale()+scale_color_manual(values=cols,labels=model_labels)+scale_fill_manual(values=cols,labels=model_labels)+labs(x='Count',y='Simulated minus observed CDF')+base+theme(legend.position='none')
saveplot(p,'cran_mark_gof',7.4,6.2)
p<-ggplot(D,aes(threshold,1-median,color=model,fill=model))+geom_ribbon(aes(ymin=1-hi,ymax=1-lo),alpha=.2,color=NA)+geom_line()+geom_line(aes(y=1-observed),color='black')+facet_grid(model~metric,scales='free_y',labeller=labeller(metric=distribution_labels,model=as_labeller(model_labels)))+count_scale()+scale_color_manual(values=cols,labels=model_labels)+scale_fill_manual(values=cols,labels=model_labels)+labs(x='Count',y='Survival function')+base+theme(legend.position='none')
saveplot(p,'cran_mark_distributions_all',7.4,if(length(laws)>4)8.7 else 8)
keys<-c('edges','transitivity','max_outdegree','indegree_centralization','pagerank_gini','mean_betweenness')
p<-ggplot(subset(T,metric %in% keys),aes(date,median,color=model,fill=model))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.12,color=NA)+geom_line()+geom_line(aes(y=observed,group=metric),color='black')+facet_wrap(~metric,scales='free_y',ncol=2,labeller=as_labeller(c(edges='Edges',transitivity='Transitivity',max_outdegree='Maximum out-degree',indegree_centralization='In-degree centralisation',pagerank_gini='PageRank Gini',mean_betweenness='Mean betweenness')))+scale_color_manual(values=cols,labels=model_labels)+scale_fill_manual(values=cols,labels=model_labels)+labs(x=NULL,y=NULL,color=NULL,fill=NULL)+base
saveplot(p,'cran_mark_trajectories',7.4,7)
p<-ggplot(E,aes(model,median,color=model))+geom_hline(aes(yintercept=observed),linetype=2,color='grey50')+geom_linerange(aes(ymin=lo,ymax=hi))+geom_point()+facet_wrap(~metric,scales='free_y',ncol=2,labeller=as_labeller(c(births='Births',additions='Edges added',deletions='Edges deleted',releasing_lines='Releasing lines')))+scale_color_manual(values=cols,labels=model_labels)+scale_x_discrete(labels=model_labels)+labs(x=NULL,y='Mean per update')+base+theme(legend.position='none',axis.text.x=element_text(angle=25,hjust=1))
saveplot(p,'cran_event_sizes',7.4,4.8)
}
if(mode %in% c('all','joint')&&dir.exists(joint)){
 origin_data<-readRDS(file.path(joint,'origins.rds'));obs<-utils::read.csv(file.path(path,'../../../application/cran/data/joint_forecast/observed.csv'))
 specs<-c('M0_CS-pair','M1_CS-pair','M2_CS-pair','M2_PA','M0_TERGM');labels<-setNames(c('M0 / CS-pair','M1b / CS-pair','M2 / CS-pair','M2 / PA','M0 / TERGM'),specs)
 kinds<-c('Migration','Triangle','Two-star','Three-star');old_kind<-c('migration','triangle','star2','star3');rows<-curves<-origin_rows<-list();ri<-si<-oi<-0L
 for(spec in specs){values<-array(NA_real_,c(24,32,4));for(i in 1:24)for(rep in 1:32){dest<-file.path(joint,sprintf('%s_%02d_%02d.rds',spec,i,rep));if(!file.exists(dest))stop('Incomplete joint simulations: ',dest);values[i,rep,]<-readRDS(dest)$wait}
  for(k in 1:4){observed_wait<-obs$days[obs$kind==old_kind[k]];stopifnot(length(observed_wait)==24);p<-rowMeans(values[,,k]<=1);interval<-t(apply(values[,,k],1,stats::quantile,c(.05,.5,.95)))
   ri<-ri+1L;rows[[ri]]<-data.frame(model=spec,target=kinds[k],observed_RM=mean(observed_wait),predicted_RM=mean(values[,,k]),Brier=mean((p-as.integer(observed_wait<=1))^2),coverage=sum(observed_wait>=interval[,1]&observed_wait<=interval[,3]),within_one_day=mean(p))
   for(i in 1:24){oi<-oi+1L;origin_rows[[oi]]<-data.frame(model=spec,target=kinds[k],date=as.Date('2019-01-01')+origin_data[[i]]$origin,observed=observed_wait[i],lo=interval[i,1],median=interval[i,2],hi=interval[i,3])}
   grid<-seq(0,14,length.out=281);si<-si+1L;curves[[si]]<-data.frame(model=spec,target=kinds[k],day=grid,observed=vapply(grid,function(t)mean(observed_wait>t),0),predicted=vapply(grid,function(t)mean(values[,,k]>t),0))
  }
 }
 J<-do.call(rbind,rows);C<-do.call(rbind,curves);O<-do.call(rbind,origin_rows);for(n in c('J','C','O'))utils::write.csv(get(n),file.path(joint,paste0(n,'.csv')),row.names=FALSE)
 C$model<-factor(C$model,specs);O$model<-factor(O$model,specs)
 jcols<-setNames(c('#0072B2','#009E73','#D55E00','#E69F00','#CC79A7'),specs)
 p<-ggplot(subset(C,target=='Migration'),aes(day,predicted,color=model))+geom_step()+geom_step(aes(y=observed),color='black')+scale_color_manual(values=jcols,labels=labels)+labs(x='Days to next migration',y='Survival probability',color=NULL)+base
 curve<-p
 intervals<-ggplot(subset(O,target=="Migration"&model %in% c("M2_CS-pair","M0_TERGM")),aes(date,median,color=model))+geom_linerange(aes(ymin=lo,ymax=hi),position=position_dodge(width=2))+geom_point(position=position_dodge(width=2),size=1)+geom_point(aes(y=observed),shape=4,color="black")+scale_color_manual(values=jcols,labels=labels)+labs(x=NULL,y="Days to next migration",color=NULL)+base+theme(axis.text.x=element_text(angle=45,hjust=1),legend.text=element_text(size=7))
 legend_plot<-ggplotGrob(curve+guides(color=guide_legend(nrow=1))+theme(legend.text=element_text(size=8)))
 legend_index<-which(legend_plot$layout$name=='guide-box-bottom');if(!length(legend_index))legend_index<-which(legend_plot$layout$name=='guide-box')
 shared_legend<-legend_plot$grobs[[legend_index[1]]]
 grDevices::pdf(file.path(figures,"cran_migration_gof.pdf"),width=7.4,height=4.4)
 grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2,heights=grid::unit(c(1,.34),c('null','inches')))))
 print(curve+theme(legend.position='none'),vp=grid::viewport(layout.pos.row=1,layout.pos.col=1));print(intervals+theme(legend.position='none'),vp=grid::viewport(layout.pos.row=1,layout.pos.col=2))
 grid::pushViewport(grid::viewport(layout.pos.row=2,layout.pos.col=1:2));grid::grid.draw(shared_legend);grid::popViewport();grDevices::dev.off()
 p<-ggplot(subset(C,target!='Migration'&day<=2),aes(day,predicted,color=model))+geom_step()+geom_step(aes(y=observed),color='black')+facet_wrap(~target,nrow=1)+scale_color_manual(values=jcols,labels=labels)+labs(x='Days to next motif creation',y='Survival probability',color=NULL)+base
 saveplot(p,'cran_motif_waits',7.4,3.7)
}
