# Run from the workspace root. Uses the existing calendar, update law, and paths.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
output <- args[1]
reps <- as.integer(args[2])
stopifnot(reps >= 1L, reps <= 32L)
dir.create(output, recursive=TRUE, showWarnings=FALSE)
pkgload::load_all("hawkesNet", quiet=TRUE, debug=FALSE)
paper <- "paper_1/6_full_rewrite"
experiments <- file.path(paper, "diagnostics")
saved <- file.path(experiments, "exponential_rates_20261105")
path <- "paper_1/5_redevelopment/application_exploration/cran_all_20260922/release_selection_10min_20260924"
time <- cran_timing_data(path)
bounds <- as.numeric(as.Date(c("2020-01-01","2024-01-01","2025-01-01","2026-01-01"))-as.Date("2019-01-01"))
d <- readRDS(file.path(saved,"gof","graph_2025.rds"))
stopifnot(max(abs(time$times-vapply(d$records, `[[`, 0., "time")))<1e-10)
g <- cran_graph_create(d$state)
prior <- d$records[time$times<d$origin]
reverse <- lapply(rev(prior), function(r) {
  a <- r$add; r$add <- r$delete; r$delete <- a
  r$births <- r$sources <- character(); r$clocks <- numeric(); r
})
invisible(cran_graph_replay(g,reverse,TRUE))
history <- cran_graph_replay(g,d$records,TRUE)
legacy <- jsonlite::fromJSON(file.path(experiments,"timing_kernels_20261105","fits.json"))$results
p1 <- legacy$two_exp_M1$parameters[c(1:3,6:7,4:5)]
p2 <- legacy$two_exp_M2_clipped$parameters[c(1:3,7:8,4:6)]
names_new <- c("M2_major_in","M2_family_out")
features <- setNames(lapply(seq_along(names_new), function(k) {
  x <- lapply(history,function(h) h$topology[[k]])
  get <- function(key) unlist(lapply(x, `[[`, key), use.names=FALSE)
  ptr <- function(key) as.integer(c(0,cumsum(lengths(lapply(x,`[[`,key)))))
  ans <- time
  ans$tp <- ptr("id"); ans$cp <- ptr("cu"); ans$tu <- get("id"); ans$cu <- get("cu")
  ans$tz <- get("summary"); ans$cz <- get("cz")
  ans$initial <- attr(history,"timing_initial")[[k]]
  length(ans$initial) <- max(length(ans$initial),1L+max(ans$tu,ans$cu))
  ans$initial[is.na(ans$initial)] <- 0
  ans$frozen <- unlist(lapply(history,`[[`,"summary"), use.names=FALSE)
  count <- diff(ans$tp)
  ans$tx <- rep(log(count),count); ans$tw <- rep(1/count,count)
  ans
}), names_new)
times <- time$times
fit_m1 <- loglik_cran_timing(time,p1,"M1",bounds[1],bounds[2],TRUE)
stopifnot(abs(fit_m1$loglik-legacy$two_exp_M1$training_loglik)<1e-6)
nested <- c(p1[1:6],0,p1[7])
for(x in features) {
  fit0 <- loglik_cran_timing(x,nested,"M2",bounds[1],bounds[2],TRUE)
  stopifnot(abs(fit0$loglik-fit_m1$loglik)<1e-6,
            max(abs(fit0$intensity-fit_m1$intensity))<1e-8)
}
message("Observed topology replay and M1b nesting verified")
# Descriptive screen: changing anchor summaries while the relevant
# kernel retains at least exp(-1) of its value at the release time.
opportunities <- function(x) {
  last <- rep(-Inf,length(x$initial)); z <- x$initial
  nshort <- nlong <- total <- positive <- negative <- 0L
  for(i in which(times<bounds[2])) {
    if(x$cp[i]<x$cp[i+1]) {
      j <- seq.int(x$cp[i]+1L,x$cp[i+1]); u <- x$cu[j]+1L
      nz <- x$cz[j]!=z[u]
      lag <- times[i]-last[u]
      if(times[i]>=bounds[1]) {
        total <- total+sum(nz)
        nshort <- nshort+sum(nz & lag<=exp(p1[4]))
        nlong <- nlong+sum(nz & lag<=exp(p1[5]))
        positive <- positive+sum(nz & lag<=exp(p1[5]) & x$cz[j]>z[u])
        negative <- negative+sum(nz & lag<=exp(p1[5]) & x$cz[j]<z[u])
      }
      z[u] <- x$cz[j]
    }
    last[x$tu[seq.int(x$tp[i]+1L,x$tp[i+1])]+1L] <- times[i]
  }
  c(changes=total,short=nshort,long=nlong,long_gain=positive,long_loss=negative)
}
# Project the optional score variation against all M1b nuisance parameters.
screen <- function(x) {
  out <- function(p) loglik_cran_timing(x,p,"M2",bounds[1],bounds[2],TRUE)
  U <- out(nested)$gradient[7]
  keep <- times>=bounds[1] & times<bounds[2]
  D <- vapply(seq_along(nested), function(j) {
    a <- b <- nested; a[j] <- a[j]+1e-5; b[j] <- b[j]-1e-5
    (log(out(a)$intensity[keep])-log(out(b)$intensity[keep]))/(2e-5)
  }, numeric(sum(keep)))
  J <- sum(lm.fit(D[,-7],D[,7])$residuals^2)
  c(opportunities(x),score=unname(U),predicted_SE=1/sqrt(J),standardized_score=unname(U/sqrt(J)))
}
screens <- as.data.frame(t(vapply(c(list(M2_family_in=time),features),screen,numeric(8))))
write.csv(screens,file.path(output,"topology_screen.csv"))
print(screens)

data_list <- c(list(M0=time,M1b=time,M2_family_in=time),features)
p0 <- fit_cran_timing(time,p1[1:5],"M0",bounds[1],bounds[2])$parameters
params <- list(M0=p0,M1b=p1,M2_family_in=p2)
fits <- list()
for(model in names_new) {
  start <- nested; start[7] <- sign(screens[model,"score"])*2
  all <- lapply(list(nested,start), function(p) {
    fit_cran_timing(features[[model]],p,"M2",bounds[1],bounds[2])
  })
  best <- all[[which.max(vapply(all,`[[`,0.,"loglik"))]]
  grad <- best$gradient
  at_bound <- abs(best$parameters[7])>=10-1e-7
  if(at_bound && grad[7]*best$parameters[7]>0) grad[7] <- 0
  stopifnot(best$convergence==0, max(abs(grad))/sum(times>=bounds[1]&times<bounds[2])<2e-6)
  params[[model]] <- best$parameters
  fits[[model]] <- all
  message(model," fit: gamma_Delta = ",round(best$parameters[7],3),
          ", start difference ",format(diff(vapply(all,`[[`,0.,"loglik")),digits=3))
}
fixed_poisson <- function(mu,a,b){
  idx <- times>=a & times<b
  sum(log(mu*time$calendar[floor(times[idx])+1L]))-mu*sum(time$calendar[seq.int(a+1L,b)])
}
mu <- sum(times>=bounds[1]&times<bounds[2])/sum(time$calendar[seq.int(bounds[1]+1L,bounds[2])])
results <- rbind(data.frame(model="Calendar Poisson",gamma_Delta=NA,
  train=fixed_poisson(mu,bounds[1],bounds[2]),ll2024=fixed_poisson(mu,bounds[2],bounds[3]),
  ll2025=fixed_poisson(mu,bounds[3],bounds[4])),
  do.call(rbind,lapply(names(params),function(m) {
    p <- params[[m]]
    type <- if(m=="M0")"M0"else if(m=="M1b")"M1"else"M2"
    ll <- vapply(1:3,function(j) loglik_cran_timing(data_list[[m]],p,type,bounds[j],bounds[j+1L])$loglik,0.)
    data.frame(model=m,gamma_Delta=if(type=="M2")p[7]else NA,train=ll[1],ll2024=ll[2],ll2025=ll[3])
  })))
write.csv(results,file.path(output,"timing.csv"),row.names=FALSE)
saveRDS(list(params=params,fits=fits,screens=screens),file.path(output,"fits.rds"))
print(results)

origins <- readRDS(file.path(saved,"joint","origins.rds"))
fit <- readRDS(file.path(saved,"fits","CS-pair_train.rds"))
physical <- function(p) c(exp(p[1:5]),p[8],p[6:7])
history_at <- lapply(origins,function(o) history[times<o$origin])
# Validate that the same graph and history produce the fitted intensity at
# every simulation origin, independent of the thinning recurrence.
for(m in names_new){
  p <- params[[m]]
  discrepancies <- vapply(seq_along(origins),function(i){
    o <- origins[[i]]
    x <- cran_joint_path(o$state,history_at[[i]],o$origin,0,fit$law,
         fit$parameters,fit$release_parameters,physical(p),time$calendar,m)
    y <- loglik_cran_timing(data_list[[m]],p,"M2",bounds[1],o$origin)$end_intensity
    abs(x$initial_intensity-y)
  },0.)
  stopifnot(max(discrepancies)<1e-8)
  message(m,": max simulation origin intensity error ",format(max(discrepancies),digits=3))
}

tasks <- expand.grid(rep=seq_len(reps),origin=seq_along(origins))
observed <- read.csv("paper_1/5_redevelopment/application/cran/data/joint_forecast/observed.csv")
use <- c(Migration="migration",Triangle="triangle")
read_path <- function(m,i,r) readRDS(file.path(saved,"joint",sprintf("%s_%02d_%02d.rds",m,i,r)))
one <- function(i,r,model,x) {
  o <- origins[[i]]$origin
  obs <- observed[observed$origin==i-1L,]
  w <- obs$days[match(use,obs$kind)]
  # A forecast may stop at the last diagnostic hit, so use counts up to
  # each first hit only and distinguish censoring at the horizon.
  data.frame(model=model,origin=i,rep=r,target=names(use),observed=w,
             simulated=unname(x$wait[names(use)]),
             observed_updates=vapply(w,function(v)sum(times>=o & times<=o+v+1e-9),0L),
             simulated_updates=vapply(x$wait[names(use)],function(v)sum(x$times<=o+v+1e-9),0L))
}
all <- lapply(c("M0_CS-pair","M1_CS-pair","M2_CS-pair","M0_TERGM",names_new),function(m) {
  # Pilot variants share the archived M1b seed for each origin/replicate.
  run <- function(j) {
    i <- tasks$origin[j]; r <- tasks$rep[j]
    if(m %in% names_new){
      o <- origins[[i]]; set.seed(20261105L+200000L+100L*i+r)
      x <- cran_joint_path(o$state,history_at[[i]],o$origin,14,fit$law,
           fit$parameters,fit$release_parameters,physical(params[[m]]),time$calendar,m)
    }else x <- read_path(m,i,r)
    one(i,r,m,x)
  }
  ans <- parallel::mclapply(seq_len(nrow(tasks)),run,mc.cores=4L)
  if(any(vapply(ans,inherits,FALSE,"try-error"))) stop("Pilot failed: ",m)
  message("Forward paths completed: ",m)
  do.call(rbind,ans)
})
paths <- do.call(rbind,all)
write.csv(paths,file.path(output,"paths.csv"),row.names=FALSE)
summary <- do.call(rbind,lapply(split(paths,list(paths$model,paths$target),drop=TRUE),function(x) {
  p <- aggregate(cbind(pred=x$simulated<=1,obs=x$observed<=1)~origin,x,mean)
  data.frame(model=x$model[1],target=x$target[1],observed_RM=mean(x$observed),
     predicted_RM=mean(x$simulated),Brier=mean((p$pred-p$obs)^2),
     pred_within_day=mean(p$pred),observed_within_day=mean(p$obs),
     predicted_first_update=mean(x$simulated_updates==1 & x$simulated<14),
     observed_first_update=mean(x$observed_updates==1 & x$observed<14),
     median_simulated_updates=median(x$simulated_updates),
     median_observed_updates=median(x$observed_updates))
}))
write.csv(summary,file.path(output,"forward.csv"),row.names=FALSE)
print(summary,row.names=FALSE)
