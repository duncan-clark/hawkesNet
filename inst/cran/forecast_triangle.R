# Run from the project root. The two-exponential adjusted M0 clock is common
# to all mark laws. Its RNG is separate, so marks cannot shift the clock.
args <- commandArgs(TRUE);stopifnot(length(args)>=1L)
output <- args[1];dir.create(output,recursive=TRUE,showWarnings=FALSE)
fit_directory <- if(length(args)>4L)normalizePath(args[5],mustWork=TRUE)else NULL
law_names <- if(is.null(fit_directory))c("CS-pair","TERGM","CS-pair-triangle")else c("CS-pair","TERGM","CS-pair-migration-add")
laws <- if(length(args)>1L)match.arg(args[2],c(law_names,"CS-pair-diamond","CS-pair-add","CS-pair-migration","CS-pair-migration-add","all"))else "all"
if(laws=="CS-pair-migration")laws<-"CS-pair"
if(laws=="all")laws<-law_names
diamond_run <- "CS-pair-diamond" %in% laws
first_only <- length(args)>3L && identical(args[4],"first")
if(diamond_run)law_names<-c(law_names,"CS-pair-diamond")
additions_law <- intersect(laws,c("CS-pair-add","CS-pair-migration-add"))
additions_run <- length(additions_law)>0L
 if(additions_run)law_names<-c("CS-pair","TERGM",additions_law)
check_diamonds <- first_only||diamond_run||additions_run
reps <- if(length(args)>2L)as.integer(args[3])else 32L
stopifnot(reps>=2L)
if(dir.exists("hawkesNet"))pkgload::load_all("hawkesNet",quiet=TRUE,debug=FALSE)else library(hawkesNet)
diagnostics <- "paper_1/6_full_rewrite/diagnostics"
archive <- file.path(diagnostics,"exponential_rates_20261105")
triangle_archive <- file.path(diagnostics,"triangle_channel_20261106")
adjusted <- readRDS(file.path(diagnostics,"calendar_gap_pilot_20261106","fit.rds"))
origins <- readRDS(file.path(archive,"joint","origins.rds"))
p <- unname(exp(adjusted$adjusted[1:5])); calendar <- adjusted$calendar
batches <- readRDS(file.path(archive,"gof","graph_2025.rds"))$records
observed_times <- vapply(batches,function(r)r$time,0.)
triangle_counts <- read.csv(file.path(diagnostics,"transitive_triples_20261106","by_update.csv"))
observed_triples <- triangle_counts$time[triangle_counts$created>0]
if(check_diamonds){
 diamond_counts<-read.csv(file.path(diagnostics,"pair_motif_audit_20261106","by_update.csv"))
 stopifnot(identical(triangle_counts$date,diamond_counts$date),identical(triangle_counts$additions,diamond_counts$additions))
}
observed <- read.csv("paper_1/5_redevelopment/application/cran/data/joint_forecast/observed.csv")
H <- 14
make_clock <- function(o,r,i){
 set.seed(20261106L+1000L*i+r)
 t <- o; q <- vapply(p[4:5],function(tau)sum(exp(-(t-observed_times[observed_times<t])/tau)),0.)
 c <- p[2:3]/p[4:5]; times <- numeric()
 while(t<o+H){
  day <- floor(t); boundary <- min(day+1,o+H)
  base <- p[1]*calendar[day+1]; bound <- base+sum(c*q)
  next_time <- min(boundary,t+rexp(1,bound))
  q <- q*exp(-(next_time-t)/p[4:5]); t <- next_time
  if(t<boundary && runif(1)<(base+sum(c*q))/bound) {times<-c(times,t); q<-q+1}
  if(length(times)>10000L)stop("Clock event guard reached")
 }
 times
}
cores <- as.integer(Sys.getenv("HAWKESNET_JOINT_CORES","4"))
targets <- c("Migration","Triangle","Transitive")
dates <- as.Date("2019-01-01")+vapply(origins,function(o)o$origin,0.)
for(law in laws){
 path <- switch(law,`CS-pair-triangle`=if(file.exists(file.path(output,"CS-pair-triangle_train.rds")))output else triangle_archive,
                `CS-pair-diamond`=output,`CS-pair-add`=output,`CS-pair-migration-add`=output,file.path(archive,"fits"))
 if(!is.null(fit_directory))path<-fit_directory
 fit <- readRDS(file.path(path,paste0(law,"_train.rds")))
 hawkesNet:::.cran_check_fit(fit)
 location <- file.path(output,paste0(if(first_only)"first_"else"paths_",law,"_",reps))
 dir.create(location,showWarnings=FALSE)
 run_origin <- function(i){
  file <- file.path(location,sprintf("origin_%02d.rds",i))
  if(file.exists(file))return(readRDS(file))
  o <- origins[[i]]; t <- o$origin
  if(check_diamonds){
   j<-which(triangle_counts$time>=t & triangle_counts$time<t+H)[1L]
   observed_diamond<-if(is.na(j))0L else diamond_counts$unlinked_shared_dependency[j]
  }
  obs <- observed[observed$origin==i-1L,]
  w <- c(obs$days[match(c("migration","triangle"),obs$kind)],
         min(c(H,observed_triples[observed_triples>=t & observed_triples<t+H]-t)))
  actual_rank <- vapply(w,function(z)if(z<H)sum(observed_times>=t & observed_times<=t+z+1e-9)else NA_integer_,0L)
  rows <- lapply(seq_len(reps),function(r) {
   clock <- make_clock(t,r,i)
   set.seed(20261106L+1000000L+1000L*i+r)
   sim <- hawkesNet:::cran_shared_clock_waits(o$state,if(first_only)head(clock,1L)else clock,t,H,hawkesNet:::.cran_sampler_law(fit),fit$parameters,fit$release_parameters)
   if(first_only)return(data.frame(law=law,origin=i,rep=r,observed=observed_diamond,predicted=sim$first_diamond))
   z<-data.frame(law=law,origin=i,date=dates[i],rep=r,target=targets,
     observed=w,wait=unname(sim$wait),observed_rank=actual_rank,rank=unname(sim$rank),
     first_update=min(c(H,clock-t)),n_clock=length(clock))
   if(check_diamonds){z$observed_diamond<-observed_diamond;z$first_diamond<-sim$first_diamond}
   z
  })
  z <- do.call(rbind,rows)
  saveRDS(z,file);message(law," origin ",i," complete");z
 }
 all <- parallel::mclapply(seq_along(origins),run_origin,mc.cores=cores,mc.preschedule=FALSE)
 if(any(vapply(all,function(z)!is.data.frame(z),FALSE)))stop("One or more forecast origins failed")
 z<-do.call(rbind,all)
 if(check_diamonds){
   first<-if(first_only)z else data.frame(law=law,origin=z$origin[z$target=="Triangle"],rep=z$rep[z$target=="Triangle"],
           observed=z$observed_diamond[z$target=="Triangle"],predicted=z$first_diamond[z$target=="Triangle"])
   write.csv(first,file.path(output,paste0("first_diamond_",law,"_",reps,".csv")),row.names=FALSE)
 }
 if(!first_only)write.csv(z[,setdiff(names(z),c("observed_diamond","first_diamond"))],
       file.path(output,paste0("paths_",law,"_",reps,".csv")),row.names=FALSE)
}
if(first_only)quit(save="no",status=0L)
paths <- paste0("paths_",law_names,"_",reps,".csv")
files <- file.path(output,paths)
if(is.null(fit_directory)&&(diamond_run||additions_run)){missing<-!file.exists(files);files[missing]<-file.path(triangle_archive,paths[missing])}
files<-files[file.exists(files)]
raw <- do.call(rbind,lapply(files,read.csv))
groups <- split(raw,list(raw$law,raw$target),drop=TRUE)
cmp <- do.call(rbind,lapply(groups,function(z){
 by_origin <- lapply(split(z,z$origin),function(a) {
  q <- quantile(a$wait,c(.05,.95));y<-a$observed[1];p<-mean(a$wait<=1)
  c(loss=(p-(y<=1))^2,covered=y>=q[1]&&y<=q[2],var=var(a$wait)/nrow(a))
 })
 ss <- do.call(rbind,by_origin)
 data.frame(law=z$law[1],target=z$target[1],observed_RM=mean(z$observed),
    predicted_RM=mean(z$wait),RM_MC_SE=sqrt(sum(ss[,"var"]))/nrow(ss),
    Brier=mean(ss[,"loss"]),coverage=sum(ss[,"covered"]),
    observed_within_day=mean(z$observed<=1),predicted_within_day=mean(z$wait<=1),
    observed_first_update=mean(!is.na(z$observed_rank)&z$observed_rank==1),
    predicted_first_update=mean(!is.na(z$rank)&z$rank==1))
}))
write.csv(cmp,file.path(output,paste0("forecast_",reps,".csv")),row.names=FALSE)
print(cmp,row.names=FALSE)
if(all(law_names %in% raw$law)){
  sort_paths <- function(z)z[order(z$origin,z$rep,z$target),]
  new_law <- if(additions_run)additions_law else if(diamond_run)"CS-pair-diamond"else"CS-pair-triangle"
  new <- sort_paths(raw[raw$law==new_law,])
  paired <- lapply(setdiff(law_names,new_law),function(law){
    old <- sort_paths(raw[raw$law==law,])
    stopifnot(identical(new$origin,old$origin),identical(new$rep,old$rep),
              identical(new$target,old$target),identical(new$observed,old$observed),
              identical(new$first_update,old$first_update),identical(new$n_clock,old$n_clock))
    do.call(rbind,lapply(targets,function(target){
      i <- new$target==target; difference<-new$wait[i]-old$wait[i]
      variation <- tapply(difference,new$origin[i],var)/reps
      data.frame(benchmark=law,target=target,mean_difference=mean(difference),
                 MC_SE=sqrt(sum(variation))/length(variation))
    }))
  })
  write.csv(do.call(rbind,paired),file.path(output,paste0("paired_differences_",reps,".csv")),row.names=FALSE)
}
