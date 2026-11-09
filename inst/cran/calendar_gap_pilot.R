# Exploratory calendar sensitivity, separate from the M2 topology pilots.
# Run from workspace root: Rscript hawkesNet/inst/cran/calendar_gap_pilot.R OUT REPS
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
out <- args[1]; reps <- as.integer(args[2]); stopifnot(reps>=1L,reps<=32L)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pkgload::load_all("hawkesNet",quiet=TRUE,debug=FALSE)
paper <- "paper_1/6_full_rewrite"
archive <- file.path(paper,"diagnostics/exponential_rates_20261105")
path <- "paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
d <- cran_timing_data(path)
origin <- as.Date("2019-01-01")
date <- format(origin+seq_along(d$calendar)-1L,"%m-%d")
D <- cbind(year_end=date>="12-24"|date<="01-06",summer=date>="07-25"&date<="08-15")
edges <- as.numeric(as.Date(c("2020-01-01","2024-01-01","2025-01-01","2026-01-01"))-origin)
N <- sum(d$times>=edges[1]&d$times<edges[2])
p0 <- fit_cran_timing(d,c(log(3),log(.45),log(.25),log(.03),log(2)),"M0",edges[1],edges[2])$parameters
calendar_data <- function(w) { x <- d; x$calendar <- d$calendar*exp(drop(D%*%w)); x }
p_last <- val <- NULL
evaluate <- function(p){
 if(identical(p,p_last)) return(val)
 x <- calendar_data(p[6:7])
 z <- loglik_cran_timing(x,p[1:5],"M0",edges[1],edges[2],TRUE)
 ids <- which(d$times>=edges[1]&d$times<edges[2]); day <- floor(d$times[ids])+1L
 e <- seq.int(edges[1]+1L,edges[2]); mu <- exp(p[1])
 wgrad <- colSums(D[day,,drop=FALSE]*mu*x$calendar[day]/z$intensity[ids]) -
   colSums(D[e,,drop=FALSE]*mu*x$calendar[e])
 val <<- list(loglik=z$loglik,gradient=c(z$gradient,wgrad))
 p_last <<- p; val
}
fit <- optim(c(p0,-2,-1),function(p)-evaluate(p)$loglik/N,function(p)-evaluate(p)$gradient/N,
             method="L-BFGS-B",lower=c(log(1e-8),log(1e-8),log(1e-8),log(1/144),0,-8,-8),
             upper=c(log(1500),log(1000),log(1000),0,log(100),0,0),
             control=list(factr=100,pgtol=1e-8,maxit=500))
stopifnot(fit$convergence==0, max(abs(evaluate(fit$par)$gradient))/N<1e-6)
x <- calendar_data(fit$par[6:7])
ll <- function(x,p) vapply(1:3,function(i)loglik_cran_timing(x,p,"M0",edges[i],edges[i+1L])$loglik,0.)
cmp <- rbind(original=ll(d,p0),gap_calendar=ll(x,fit$par[1:5])); colnames(cmp)<-c("train","LL2024","LL2025")
print(cmp); print(setNames(exp(fit$par[6:7]),colnames(D)))
write.csv(cmp,file.path(out,"timing.csv"))
saveRDS(list(original=p0,adjusted=fit$par,calendar=x$calendar),file.path(out,"fit.rds"))
origins <- readRDS(file.path(archive,"joint","origins.rds"))
first_wait <- function(x,p,model) {
 p <- unname(exp(p)); csum <- c(0,cumsum(x$calendar))
 clock <- function(t) {
  j <- floor(t)+1L; csum[j]+(t-floor(t))*x$calendar[pmin(j,length(x$calendar))]
 }
 do.call(rbind,lapply(origins,function(o) {
  t <- o$origin; lag <- t-x$times[x$times<t]
  a <- p[2:3]*vapply(p[4:5],function(tau)sum(exp(-lag/tau)),0.)
  surv <- function(h) exp(-p[1]*(clock(t+h)-clock(t))+
                        a[1]*expm1(-h/p[4])+a[2]*expm1(-h/p[5]))
  pred <- sum(vapply(0:13,function(j)integrate(surv,j,j+1,rel.tol=1e-9)$value,0.))
  obs <- min(14,min(x$times[x$times>=t])-t); prob <- 1-surv(1)
  data.frame(model=model,date=origin+t,observed=obs,predicted=pred,
             within_day=prob,Brier=(prob-(obs<=1))^2)
 }))
}
write.csv(rbind(first_wait(d,p0,"original"),first_wait(x,fit$par[1:5],"gap_calendar")),
          file.path(out,"first_update_exact.csv"),row.names=FALSE)
observed <- read.csv("paper_1/5_redevelopment/application/cran/data/joint_forecast/observed.csv")
tasks <- expand.grid(rep=seq_len(reps),origin=seq_along(origins))
metrics <- c(first="first",migration="Migration",triangle="Triangle")
outcome <- function(i) {
 t <- origins[[i]]$origin
 w <- observed[observed$origin==i-1L,]
 c(first=min(14,d$times[which(d$times>=t)[1]]-t),
   migration=w$days[w$kind=="migration"],triangle=w$days[w$kind=="triangle"])
}
run <- function(k,law,adjusted){
 i <- tasks$origin[k]; r <- tasks$rep[k]; o <- origins[[i]]
 if(adjusted){
  u <- fits[[law]]; seed <- 20261105L+if(law=="TERGM")500000L else 100000L
  set.seed(seed+100L*i+r)
  v <- cran_joint_path(o$state,o$history,o$origin,14,u$law,u$parameters,u$release_parameters,
                       c(exp(fit$par[1:5]),0,0,0),x$calendar,"M0")
 }else{
  v <- readRDS(file.path(archive,"joint",sprintf("M0_%s_%02d_%02d.rds",law,i,r)))
 }
 waits <- c(first=if(length(v$times))v$times[1]-o$origin else 14,
            migration=v$wait[["Migration"]],triangle=v$wait[["Triangle"]])
 data.frame(law=law,calendar=if(adjusted)"gap_calendar"else"original",
            origin=i,rep=r,target=names(metrics),pred=unname(waits),obs=unname(outcome(i)))
}
fits <- setNames(lapply(c("CS-pair","TERGM"),function(law)readRDS(file.path(archive,"fits",paste0(law,"_train.rds")))),c("CS-pair","TERGM"))
z <- do.call(rbind,lapply(names(fits),function(law)do.call(rbind,lapply(c(FALSE,TRUE),function(adjusted){
 a <- parallel::mclapply(seq_len(nrow(tasks)),run,law=law,adjusted=adjusted,mc.cores=4L)
 if(any(vapply(a,inherits,FALSE,"try-error")))stop("Calendar path failed")
 message("Forward paths completed: ",law," / adjusted = ",adjusted)
 do.call(rbind,a)
}))))
write.csv(z,file.path(out,"paths.csv"),row.names=FALSE)
compare <- do.call(rbind,lapply(split(z,list(z$law,z$calendar,z$target),drop=TRUE),function(a){
 b <- aggregate(cbind(p=a$pred<=1,y=a$obs<=1)~origin,a,mean)
 data.frame(law=a$law[1],calendar=a$calendar[1],target=a$target[1],observed_RM=mean(a$obs),
            predicted_RM=mean(a$pred),Brier=mean((b$p-b$y)^2),predicted_within_day=mean(b$p),observed_within_day=mean(b$y))
}))
write.csv(compare,file.path(out,"forward.csv"),row.names=FALSE)
print(compare,row.names=FALSE)
