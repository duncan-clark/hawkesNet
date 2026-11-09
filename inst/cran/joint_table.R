# Rscript hawkesNet/inst/cran/joint_table.R JOINT_DIRECTORY MAIN_TEX
args<-commandArgs(TRUE);if(length(args)!=2L)stop('Supply joint output and existing main.tex')
x<-utils::read.csv(file.path(args[1],'J.csv'));main<-args[2];s<-paste(readLines(main,warn=FALSE),collapse='\n')
models<-c('M0_CS-pair','M1_CS-pair','M2_CS-pair','M2_PA','M0_TERGM');targets<-c('Migration','Triangle','Two-star','Three-star')
labels<-setNames(c('M0 / CS-pair','M1 / CS-pair','M2 / CS-pair','M2 / PA','M0 / TERGM'),models)
rows<-character();for(target in targets)for(model in models){z<-x[x$target==target&x$model==model,];stopifnot(nrow(z)==1L);rows<-c(rows,paste0(target,' & ',labels[[model]],' & ',sprintf('%.3f',z$observed_RM),' & ',sprintf('%.3f',z$predicted_RM),' & ',sprintf('%.3f',z$Brier),' & ',z$coverage,'/24\\\\'))}
text<-paste(c(r"{\begin{table}[htbp]
\centering\small
\caption{Joint forecasts from 24 origins in 2025, with 32 paths per model and origin. RM is the mean restricted wait $\E\{\min(W,14)\}$ in days; the Brier loss concerns occurrence within one day. Coverage counts observed restricted waits inside pointwise central 90\% predictive intervals and is descriptive. TERGM uses the same two-exponential calendar Hawkes M0 as the M0/CS-pair comparison.}\label{tab:cran_joint}
\begin{tabular}{llrrrr}
\toprule
Target & Model & Observed RM & Predicted RM & Brier & Coverage\\
\midrule}",rows,r"{\bottomrule
\end{tabular}
\end{table}}"),collapse='\n')
position<-regexpr('\\label{tab:cran_joint}',s,fixed=TRUE)[1];stopifnot(position>0)
a<-max(gregexpr('\\begin{table}',substr(s,1,position),fixed=TRUE)[[1]]);b<-regexpr('\\end{table}',substr(s,position,nchar(s)),fixed=TRUE)[1]+position-1L
s<-paste0(substr(s,1,a-1L),text,substr(s,b+nchar('\\end{table}'),nchar(s)))
writeLines(s,main)
