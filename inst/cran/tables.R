# Rscript hawkesNet/inst/cran/tables.R FIT_DIRECTORY MAIN_TEX [MIGRATION_ADDITION_FIT_DIRECTORY]
args<-commandArgs(TRUE);if(!length(args) %in% 2:3)stop('Supply fits and existing main.tex, with an optional combined-fit directory')
fits<-args[1];main<-args[2];combined<-'CS-pair-migration-add'
s<-paste(readLines(main,warn=FALSE),collapse='\n')
include_combined<-length(args)==3L
if(!include_combined&&grepl(paste0(combined,' &'),s,fixed=TRUE))stop('Supply the combined-fit directory to refresh the current manuscript comparison')
laws<-c('PA','CS-1','CS-pair',if(include_combined)combined,'TERGM')
fit_path<-function(law,window)file.path(if(law==combined)args[3]else fits,paste0(law,'_',window,'.rds'))
report_dir<-if(include_combined)args[3]else fits
f<-setNames(lapply(laws,function(law)readRDS(fit_path(law,'full'))),laws)
replace_table<-function(label,new){
 at<-regexpr(paste0('\\label{',label,'}'),s,fixed=TRUE)[1];stopifnot(at>0)
 starts<-gregexpr('\\begin{table}',substr(s,1,at),fixed=TRUE)[[1]];a<-max(starts)
 end<-regexpr('\\end{table}',substr(s,at,nchar(s)),fixed=TRUE)[1]+at-1L
 s<<-paste0(substr(s,1,a-1L),new,substr(s,end+nchar('\\end{table}'),nchar(s)))
}
fmt<-function(x,digits=4L){if(is.na(x))return('--');v<-formatC(abs(x),format='f',digits=digits);paste0(if(x<0)'$-$'else'',v)}
value<-function(fit,key,exponent=FALSE,fixed=FALSE){if(fixed)return('0$^{\\rm f}$');if(!key %in% names(fit$parameters))return('--');v<-fit$parameters[[key]];if(exponent)v<-exp(v);paste0(fmt(v),if(key %in% fit$bound_parameters)'$^\\dagger$'else'')}
row<-function(label,values)paste0(label,' & ',paste(values,collapse=' & '),'\\\\')
# Exact complete-update comparison, full-window fit and held-out train fit.
comparison<-do.call(rbind,lapply(laws,function(law){tr<-readRDS(fit_path(law,'train'));fu<-f[[law]];stopifnot(max(tr$projected_gradient,fu$projected_gradient)<=2e-5);data.frame(model=law,parameters=fu$parameter_count,full_loglik=fu$update_loglik,predictive_2024=tr$predictive[['2024']]$update_loglik,predictive_2025=tr$predictive[['2025']]$update_loglik)}))
utils::write.csv(comparison,file.path(report_dir,'comparison.csv'),row.names=FALSE)
pretty_ll<-function(v)paste0(if(v<0)'$-$'else'',formatC(abs(v),format='f',digits=1,big.mark=','))
rows<-vapply(seq_len(nrow(comparison)),function(i){x<-comparison[i,];row(if(include_combined&&x$model=='CS-pair')'CS-pair-migration'else x$model,c(x$parameters,pretty_ll(x$full_loglik),pretty_ll(x$predictive_2024),pretty_ll(x$predictive_2025)))},'')
replace_table('tab:cran_marks',paste(c(r"{\begin{table}[t]
\centering\small
\caption{Update-law log-likelihoods for 2020--2025 and predictive log-likelihoods for 2024 and 2025, evaluated with parameters estimated from 2020--2023; larger is better. All laws share the release/birth law and include popularity, migration, copying and separate continuous addition/deletion age effects. PA and the CS laws use release-type and out-degree effects on single-edit proposal means. The CS laws use directed motifs. TERGM uses a lagged-covariate transition law (Supplement~\ref{sec:supp_cran_tergm}); every entry evaluates normalised probabilities of the same complete updates. Parameter counts include the release/birth law and exclude cancelling or fixed coefficients.}\label{tab:cran_marks}
\begin{tabular}{lrrrr}
\toprule
Edit law & \shortstack{No. of\\parameters} & Log-likelihood & \shortstack[r]{2024 predictive\\log-likelihood} & \shortstack[r]{2025 predictive\\log-likelihood}\\
\midrule}",rows,r"{\bottomrule
\end{tabular}
\end{table}}"),collapse='\n'))
pa<-f[['PA']]
parows<-c(row('Reference proposal mean',c(value(pa,'log_m_add',TRUE),value(pa,'log_m_delete',TRUE))),row('Popularity exponent',c(value(pa,'pa_add'),value(pa,'pa_delete'))),row('Age decay (per year)',c(value(pa,'age_decay_add'),value(pa,'age_decay_delete'))),row('Migration',c(value(pa,'migration'),'--')),row('Copying',c(value(pa,'copy'),'--')))
patable<-paste(c(r"{\begin{table}[htbp]
\centering\small
\caption{Full-window PA edge-law estimates. Reference proposal means are the exponentiated intercepts in \eqref{eq:cran_proposal_means}; rate coefficients are in Table~\ref{tab:cran_rate_est}. $\dagger$ marks a numerical bound.}\label{tab:cran_edge_pa}
\begin{tabular}{lrr}
\toprule
Parameter & Addition & Deletion\\
\midrule}",parows,r"{\bottomrule
\end{tabular}
\end{table}}"),collapse='\n')
rate_keys<-c('mean_add_new_package','mean_add_new_major','mean_add_log_outdegree','mean_delete_log_outdegree')
rate_labels<-c('Addition: new package','Addition: new major','Addition: log one plus out-degree','Deletion: log one plus out-degree')
rate_rows<-vapply(seq_along(rate_keys),function(i)row(rate_labels[i],vapply(f[c('PA','CS-1','CS-pair')],value,'',key=rate_keys[i])),'')
rate_table<-paste(c(r"{\begin{table}[htbp]
\centering\small
\caption{Full-window coefficients for the log proposal means in \eqref{eq:cran_proposal_means}. Existing lines are the reference release class; out-degree is measured before the update. These coefficients govern proposal counts, separately from target selection.}\label{tab:cran_rate_est}
\begin{tabular}{lrrr}
\toprule
Coefficient & PA & CS-1 & CS-pair\\
\midrule}",rate_rows,r"{\bottomrule
\end{tabular}
\end{table}}"),collapse='\n')
if(grepl('\\label{tab:cran_rate_est}',s,fixed=TRUE))replace_table('tab:cran_rate_est',rate_table)else patable<-paste(patable,rate_table,sep='\n\n')
labels<-c('Intercept','In-two-stars','Out-two-stars','In-three-stars','Out-three-stars','Transitive triples','Cyclic triples','Popularity','Age decay (per year)','Migration','Copying','Reference proposal mean')
keysA<-c('addition_intercept','addition_in2star','addition_out2star','addition_in3star','addition_out3star','addition_transitive','addition_cyclic','popularity_add','age_decay_add','migration','copy','log_m_add')
keysD<-c('deletion_intercept','deletion_in2star','deletion_out2star','deletion_in3star','deletion_out3star','deletion_transitive',NA,'popularity_delete','age_decay_delete',NA,NA,'log_m_delete')
tergmA<-c('addition_intercept','addition_in_2star','addition_out_2star','addition_in_3star','addition_out_3star','addition_transitive','addition_cyclic','addition_popularity','addition_age_decay','addition_migration','addition_copy',NA)
tergmD<-c('deletion_intercept','deletion_in_2star','deletion_out_2star','deletion_in_3star','deletion_out_3star','deletion_transitive',NA,'deletion_popularity','deletion_age_decay',NA,NA,NA)
get<-function(fit,key,i){if(is.na(key))return(if(i==7L)'0$^{\\rm f}$'else'--');value(fit,key,i==12L)}
csrows<-vapply(seq_along(labels),function(i)row(labels[i],c(get(f[['CS-1']],keysA[i],i),get(f[['CS-1']],keysD[i],i),get(f[['CS-pair']],keysA[i],i),get(f[['CS-pair']],keysD[i],i),get(f[['TERGM']],tergmA[i],i),get(f[['TERGM']],tergmD[i],i))),'')
cstable<-paste(c(r"{\begin{table}[htbp]
\centering\small\setlength{\tabcolsep}{4pt}
\caption{Full-window single-edit coefficients for CS-1, CS-pair-migration (CS-pair in the table) and TERGM. Additions and deletions have separate coefficients, including ageing. Intercepts and out-star terms cancel from the exponential CS selection weights and are omitted (--). Reference proposal means are as in \eqref{eq:cran_proposal_means}. $\dagger$ marks a numerical bound and $\mathrm f$ a coefficient fixed at zero; point estimates are reported.}\label{tab:cran_edge_directed}
\begin{tabular}{lrrrrrr}
\toprule
 & \multicolumn{2}{c}{CS-1} & \multicolumn{2}{c}{CS-pair} & \multicolumn{2}{c}{TERGM}\\
Parameter & Addition & Deletion & Addition & Deletion & Addition & Deletion\\
\midrule}",csrows,r"{\bottomrule
\end{tabular}
\end{table}}"),collapse='\n')
if(grepl('\\label{tab:cran_edge_pa}',s,fixed=TRUE)){
 replace_table('tab:cran_edge_pa',patable);replace_table('tab:cran_edge_directed',cstable)
}else replace_table('tab:cran_edge_directed',paste(patable,cstable,sep='\n\n'))
pair<-f[['CS-pair']]
pairrows<-c(row(r"{$\omega_P$}",value(pair,'log_m_migration',TRUE)),row('In-two-stars',value(pair,'pair_in2star')),row('In-three-stars',value(pair,'pair_in3star')),row('Transitive triples',value(pair,'pair_transitive')),row('Cyclic triples',value(pair,'pair_cyclic',fixed=TRUE)))
pair_label<-if(grepl('\\label{tab:cran_edge_undirected}',s,fixed=TRUE))'tab:cran_edge_undirected'else'tab:cran_edge_pairs'
replace_table(pair_label,paste(c(r"{\begin{table}[htbp]
\centering\small
\caption{Full-window directed paired-migration estimates for the migration-only law, with exponential pair weights. Popularity coefficients and the new target's age effect are inherited from its single-edit channels in Table~\ref{tab:cran_edge_directed}. A constant pair intercept cancels; both out-star changes vanish and are omitted. The cyclic coefficient is fixed at zero because all candidate scores vanish in the data.}\label{tab:cran_edge_pairs}
\begin{tabular}{lr}
\toprule
Parameter & Estimate\\
\midrule}",pairrows,r"{\bottomrule
\end{tabular}
\end{table}}"),collapse='\n'))
s<-gsub('tab:cran_edge_undirected','tab:cran_edge_pairs',s,fixed=TRUE)
writeLines(s,main)
# Record independently reproducible point estimates.
estimates<-do.call(rbind,lapply(laws,function(law)data.frame(model=law,parameter=names(f[[law]]$parameters),estimate=unname(f[[law]]$parameters),bound=names(f[[law]]$parameters)%in%f[[law]]$bound_parameters)))
utils::write.csv(estimates,file.path(report_dir,'estimates.csv'),row.names=FALSE)
