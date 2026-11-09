#' Calendar Hawkes timing for complete CRAN updates
#'
#' M0 is an ordinary two-exponential calendar Hawkes model, batch adds update
#' size, M1 adds topology at occurrence and M2 revises amplitudes as topology
#' changes. The temporal kernel is held fixed across these specifications.
#' @param path Existing CRAN application cache directory.
#' @param data Prepared timing data, including daily `calendar` values from
#'   1 January 2019. For calendar effects also supply integer `calendar_phase`
#'   of the same length: 0 for ordinary days, 1 for breaks, 2 for reopening.
#' @param parameters Timing parameters on the optimisation scale: log baseline,
#'   log masses, log mean delays, then the active amplitude coefficients, and
#'   finally `log_break`, `log_reopening` when `calendar_effects=TRUE`. Fully
#'   named vectors are reordered by name. In fitting, NULL supplies a default start.
#' @param model `M0` uses two exponentials and unit amplitudes; `batch` (or
#'   `M1a`) adds batch size; `M1` adds initial topology; `M2` adds later topology.
#'   The manuscript aliases `M1b` and `M1c` automatically use separate
#'   release-type effects for both kernels, without and with initial topology.
#'   For the manuscript M2 also set `release_types=TRUE, type_by_kernel=TRUE`.
#'   `Poisson` and `one_exponential` provide unmarked timing controls.
#' @param release_types Use class-specific excitation for each releasing line.
#'   In the likelihood and fit, append `log_amp_new_package` and
#'   `log_amp_new_major` to the amplitude parameters. The existing-line
#'   multiplier is fixed at one. M0 always uses unit amplitudes.
#' @param type_by_kernel With `release_types=TRUE`, fit four log multipliers
#'   instead of two. Append the new-package and new-major coefficients for
#'   the short kernel first, then for the long kernel, using `_short` and
#'   `_long` suffixes in their names.
#' @param start,end Window boundaries in days since 1 January 2019.
#' @param details Return event contributions, interval compensators and rates.
#' @param clip Bound for later log-degree changes in M2.
#' @param control Additional L-BFGS-B optimisation controls.
#' @param calendar_effects Estimate the publication-break and reopening
#'   multipliers. During breaks the fixed calendar is multiplied by
#'   `exp(log_break)`; reopening days have level `exp(log_reopening)` relative
#'   to the baseline scale, replacing the ordinary calendar profile.
#' @return Prepared timing data; a list containing `loglik`, `gradient`,
#'   `compensator` and `end_intensity`; or a fit with parameters, bounds,
#'   convergence code and projected gradient. Fits also contain `timing`
#'   (physical parameters), `calendar` (the fitted daily profile) and
#'   `simulation_model`, ready for [simulate_cran_joint()]. The calendar
#'   conditions on all supplied phase labels, including any future break dates.
#'   The endpoint intensity is NA
#'   when the calendar ends exactly at the window boundary. Search bounds
#'   follow the manuscript and are numerical constraints, not penalties.
#' @section Data format:
#' `times` are strictly increasing days since 2019. `initial` holds initial
#' family log-degree summaries. `tz`, `tw`, `tx`, `tu` hold each releasing
#' line's summary at occurrence, weight (one divided by batch size), log batch
#' size and zero-based family ID. `tp` is a zero-based cumulative pointer of
#' length `length(times)+1` into those arrays. `cz`, `cu` and `cp` similarly
#' store post-update summary changes and their family IDs. For release-type
#' amplitudes, `release_class` labels each trigger 0 (existing line), 1 (new
#' package) or 2 (new major). All history before `start` contributes excitation;
#' its event densities are excluded from the fitting window.
#' @examples
#' d <- list(times=c(.2,.8), initial=0, tz=c(0,0), tw=c(1,1),
#'           tx=c(0,0), tu=c(0L,0L), tp=0:2, cz=numeric(), cu=integer(),
#'           cp=c(0L,0L,0L), calendar=1)
#' loglik_cran_timing(d, log(c(2,.4,.2,.1,2)), "M0", start=0, end=1)$loglik
#' @name cran_timing
NULL
#' @rdname cran_timing
#' @export
cran_timing_data<-function(path,release_types=FALSE){
 root<-file.path(path,'excitation_exploration_20260923/family_source_package')
 rd<-function(n)utils::read.csv(gzfile(file.path(root,paste0(n,'.csv.gz'))))
 b<-utils::read.csv(gzfile(file.path(path,'batches.csv.gz')));tr<-rd('triggers');ch<-rd('degree_changes');packages<-rd('packages')
 n<-nrow(b);logn<-log(tabulate(tr$batch+1L,n));tr$log_count<-logn[tr$batch+1L]
 ans<-list(times=b$time,initial=log1p(packages$initial_degree),tz=log1p(tr$degree),tw=tr$weight,tx=tr$log_count,tu=as.integer(tr$package),tp=as.integer(c(0,cumsum(tabulate(tr$batch+1L,n)))),cz=log1p(ch$degree),cu=as.integer(ch$package),cp=as.integer(c(0,cumsum(tabulate(ch$batch+1L,n)))),calendar=utils::read.csv(file.path(path,'current_fits/calendar_values.csv'))$value)
 if(release_types){
   members<-utils::read.csv(gzfile(file.path(path,'members.csv.gz')))
   member_key<-paste(members$batch,sub('@[^@]+$','',members$package))
   trigger_key<-paste(tr$batch,packages$name[match(tr$package,packages$package)])
   occurrence_key<-function(key)paste(key,ave(seq_along(key),key,FUN=seq_along))
   idx<-match(occurrence_key(trigger_key),occurrence_key(member_key))
   if(anyNA(idx)||length(idx)!=nrow(members))stop('Release members disagree with timing triggers')
   ans$release_class<-.cran_release_class(members$release_type[idx])
 }
 ans
}
.cran_release_class<-function(type){
 code<-match(type,c('dependency_update','new_package','new_major'))-1L
 if(anyNA(code))stop('Unknown CRAN release class')
 code
}
.cran_type_amplitudes<-c('log_amp_new_package','log_amp_new_major')
.cran_kernel_type_amplitudes<-c(paste0(.cran_type_amplitudes,'_short'),paste0(.cran_type_amplitudes,'_long'))
.cran_timing_layout <- function(model, release_types=FALSE, type_by_kernel=FALSE,
                                calendar_effects=FALSE) {
  model <- match.arg(model, c("M0", "batch", "M1", "M2", "M1a", "M1b", "M1c",
                            "Poisson", "one_exponential"))
  if(model == "M1a") {
    model <- "batch"
    release_types <- type_by_kernel <- FALSE
  }
  if(model %in% c("M1b", "M1c")) {
    release_types <- type_by_kernel <- TRUE
    model <- if(model == "M1b") "batch" else "M1"
  }
  k <- if(model == "Poisson") 0L else if(model == "one_exponential") 1L else 2L
  typed <- release_types && model %in% c("batch", "M1", "M2")
  amplitude <- c(switch(model, batch="delta", M1=c("gamma_B", "delta"),
                        M2=c("gamma_B", "gamma_Delta", "delta")),
                 if(typed) if(type_by_kernel) .cran_kernel_type_amplitudes else .cran_type_amplitudes)
  names <- c("log_mu", c("log_alpha_short", "log_alpha_long")[seq_len(k)],
             c("log_tau_short", "log_tau_long")[seq_len(k)], amplitude,
             if(calendar_effects) c("log_break", "log_reopening"))
  list(model=model, k=k, typed=typed, type_by_kernel=type_by_kernel,
       names=names, alpha=1L+seq_len(k), tau=1L+k+seq_len(k),
       amplitude=1L+2L*k+seq_along(amplitude),
       phase=if(calendar_effects) length(names)-1:0 else integer())
}

.cran_timing_names <- function(model, release_types=FALSE, type_by_kernel=FALSE,
                               calendar_effects=FALSE) {
  .cran_timing_layout(model,release_types,type_by_kernel,calendar_effects)$names
}

.cran_timing_parameters <- function(parameters, layout) {
  .cran_parameters(parameters,layout$names,"timing parameters")
}

.cran_validate_timing <- function(data, start, end, clip, typed, calendar_effects) {
  finite <- function(x) is.numeric(x) && all(is.finite(x))
  if(length(start)!=1L || length(end)!=1L || !finite(c(start,end)) || start<0 || end<=start)
    stop("Use finite window boundaries with 0 <= start < end")
  if(length(clip)!=1L || !is.finite(clip) || clip<0) stop("clip must be finite and nonnegative")
  if(!finite(data$times) || any(data$times<0) || is.unsorted(data$times,strictly=TRUE))
    stop("Timing times must be finite, nonnegative and strictly increasing")
  if(!finite(data$calendar) || length(data$calendar)<ceiling(end) || any(data$calendar<0))
    stop("The nonnegative daily calendar must cover the complete window")
  if(!finite(data$initial)) stop("Invalid initial topology summaries")
  n <- length(data$times)
  for(channel in c("t","c")) {
    value <- data[[if(channel=="t") "tz" else "cz"]]
    id <- data[[paste0(channel,"u")]]
    pointer <- data[[paste0(channel,"p")]]
    if(!finite(value) || !finite(id) || length(id)!=length(value) ||
       any(id!=floor(id) | id<0 | id>=length(data$initial)))
      stop("Invalid timing ",channel,"-channel summaries or family indices")
    if(!finite(pointer) || length(pointer)!=n+1L || pointer[1]!=0 ||
       any(pointer!=floor(pointer)) || is.unsorted(pointer) || tail(pointer,1)!=length(value))
      stop("Invalid timing ",channel,"-channel pointers")
  }
  if(!finite(data$tw) || !finite(data$tx) || length(data$tw)!=length(data$tz) ||
     length(data$tx)!=length(data$tz) || any(data$tw<0))
    stop("Invalid trigger weights or log batch sizes")
  mass <- diff(c(0,cumsum(data$tw))[data$tp+1L])
  if(any(abs(mass-1)>1e-8)) stop("Trigger weights must sum to one within every update")
  if(!is.null(data$frozen) && (!finite(data$frozen) || length(data$frozen)!=length(data$tz)))
    stop("Invalid frozen topology summaries")
  if(typed && is.null(data$release_class)) stop("Missing timing release classes")
  if(typed && (length(data$release_class)!=length(data$tz) || anyNA(data$release_class) ||
               any(!data$release_class %in% 0:2))) stop("Invalid timing release classes")
  if(calendar_effects && (length(data$calendar_phase)!=length(data$calendar) ||
                         anyNA(data$calendar_phase) || any(!data$calendar_phase %in% 0:2)))
    stop("calendar_phase must label each day 0 (ordinary), 1 (break), or 2 (reopening)")
  invisible(TRUE)
}

#' @rdname cran_timing
#' @export
loglik_cran_timing <- function(data, parameters, model=c("M0","batch","M1","M2"),
                               start=365, end=2557, details=FALSE, clip=log(4),
                               release_types=FALSE, type_by_kernel=FALSE,
                               calendar_effects=FALSE) {
  if(missing(model)) model <- "M0"
  if(length(model)!=1L) stop("model must name one timing specification")
  layout <- .cran_timing_layout(model,release_types,type_by_kernel,calendar_effects)
  p <- .cran_timing_parameters(parameters,layout)
  .cran_validate_timing(data,start,end,clip,layout$typed,calendar_effects)
  n <- sum(data$times<end)
  d <- data
  d$times <- data$times[seq_len(n)]
  d$tp <- data$tp[seq_len(n+1L)]
  d$cp <- data$cp[seq_len(n+1L)]
  base <- data$calendar
  phase <- if(calendar_effects) data$calendar_phase else rep(0L,length(base))
  if(calendar_effects) base[phase==2L] <- 1
  weights <- c(1,if(calendar_effects) exp(p[layout$phase]) else c(1,1))
  calendar <- base*weights[phase+1L]
  integral <- function(profile,t) {
    day <- floor(t)
    c(0,cumsum(profile))[day+1L]+(t-day)*profile[pmin(day+1L,length(profile))]
  }
  left <- pmax(c(0,d$times),start)
  right <- pmax(left,pmin(c(d$times,end),end))
  exposure <- integral(calendar,right)-integral(calendar,left)
  z <- unname(p[layout$amplitude])
  amplitude <- switch(layout$model,batch=c(0,z),M1=z,M2=z,c(0,0))
  core <- c(p[1],p[layout$alpha],amplitude)
  endpoint <- if(floor(end)<length(calendar)) calendar[floor(end)+1L] else NA_real_
  evaluate <- function(delay,detail=FALSE) {
    cran_timing_native(d,core,exp(-delay),rep(1L,layout$k),
      if(layout$model=="M2") 2L else 1L,FALSE,FALSE,clip,start,end,
      calendar[floor(d$times)+1L],exposure,detail,endpoint)
  }
  out <- evaluate(p[layout$tau], details || calendar_effects)
  g <- numeric(length(p))
  first <- seq_len(1L+layout$k)
  g[first] <- out$gradient[first]
  if(length(layout$amplitude)) g[layout$amplitude] <-
    if(layout$model=="batch") out$gradient[-seq_len(2L+layout$k)] else out$gradient[-first]
  for(j in seq_along(layout$tau)) {
    plus <- minus <- p[layout$tau]
    plus[j] <- plus[j]+1e-5
    minus[j] <- minus[j]-1e-5
    g[layout$tau[j]] <- (evaluate(plus)$loglik-evaluate(minus)$loglik)/2e-5
  }
  if(calendar_effects) {
    day <- floor(d$times)+1L
    for(j in 1:2) {
      take <- d$times>=start & phase[day]==j
      profile <- base*(phase==j)
      mass <- sum(integral(profile,right)-integral(profile,left))
      g[layout$phase[j]] <- exp(p[1]+p[layout$phase[j]])*
        (sum(base[day[take]]/out$intensity[take])-mass)
    }
  }
  if(!details) out[c("contributions","interval_compensator","intensity")] <- NULL
  out$gradient <- setNames(g,layout$names)
  out
}

#' @rdname cran_timing
#' @export
fit_cran_timing <- function(data, parameters=NULL, model=c("M0","batch","M1","M2"),
                            start=365, end=1826, clip=log(4), control=list(),
                            release_types=FALSE, type_by_kernel=FALSE,
                            calendar_effects=FALSE) {
  if(missing(model)) model <- "M0"
  if(length(model)!=1L) stop("model must name one timing specification")
  layout <- .cran_timing_layout(model,release_types,type_by_kernel,calendar_effects)
  .cran_validate_timing(data,start,end,clip,layout$typed,calendar_effects)
  n <- sum(data$times>=start & data$times<end)
  if(n==0) stop("No events in estimation window")
  if(is.null(parameters)) {
    parameters <- setNames(numeric(length(layout$names)),layout$names)
    parameters[1] <- log(max(.001,min(50,.5*n/(end-start))))
    parameters[layout$alpha] <- log(.2)
    parameters[layout$tau] <- log(c(.05,2)[seq_len(layout$k)])
  }
  parameters <- .cran_timing_parameters(parameters,layout)
  lower <- upper <- setNames(numeric(length(parameters)),layout$names)
  lower[1] <- log(.001); upper[1] <- log(50)
  lower[layout$alpha] <- log(1e-8); upper[layout$alpha] <- log(5)
  lower[layout$tau] <- c(log(1/144),0)[seq_len(layout$k)]
  upper[layout$tau] <- if(layout$k==1L) log(100) else c(0,log(100))[seq_len(layout$k)]
  for(name in intersect(layout$names,c("gamma_B","gamma_Delta","delta"))) {
    lower[name] <- c(gamma_B=-2,gamma_Delta=-10,delta=-1)[name]
    upper[name] <- c(gamma_B=2,gamma_Delta=10,delta=2)[name]
  }
  type <- grep("^log_amp_",layout$names)
  lower[type] <- -8; upper[type] <- 3
  lower[layout$phase] <- log(1e-4); upper[layout$phase] <- log(50)
  last <- value <- NULL
  evaluate <- function(p) {
    if(!identical(last,p)) {
      last <<- p
      value <<- loglik_cran_timing(data,p,model,start,end,clip=clip,
        release_types=release_types,type_by_kernel=type_by_kernel,calendar_effects=calendar_effects)
    }
    value
  }
  fit <- stats::optim(parameters,function(p) -evaluate(p)$loglik/n,
    function(p) -evaluate(p)$gradient/n,method="L-BFGS-B",lower=lower,upper=upper,
    control=utils::modifyList(list(maxit=700,factr=100,pgtol=1e-8),control))
  names(fit$par) <- layout$names
  score <- evaluate(fit$par)$gradient
  projected <- .cran_projected(fit$par,-score/n,lower,upper)
  physical <- c(mu=exp(unname(fit$par[1])),alpha_short=0,alpha_long=0,
                tau_short=.1,tau_long=1,delta=0,gamma_B=0,gamma_Delta=0)
  if(layout$k) {
    physical[c("alpha_short","alpha_long")[seq_len(layout$k)]] <- exp(fit$par[layout$alpha])
    physical[c("tau_short","tau_long")[seq_len(layout$k)]] <- exp(fit$par[layout$tau])
  }
  coefficients <- intersect(c("delta","gamma_B","gamma_Delta"),names(fit$par))
  physical[coefficients] <- fit$par[coefficients]
  physical <- c(physical,fit$par[grep("^log_amp_",names(fit$par))])
  attr(physical,"clip") <- clip
  calendar <- data$calendar
  if(calendar_effects) {
    calendar[data$calendar_phase==2L] <- 1
    calendar <- calendar*c(1,exp(fit$par[layout$phase]))[data$calendar_phase+1L]
  }
  simulation_model <- if(layout$model=="M2") "M2" else
    if(layout$model %in% c("batch","M1")) "M1" else "M0"
  list(model=model,parameters=fit$par,loglik=-fit$value*n,convergence=fit$convergence,
       window=c(start,end),clip=clip,gradient=score,projected_gradient=max(abs(projected)),
       lower=lower,upper=upper,release_types=layout$typed,type_by_kernel=layout$type_by_kernel,
       calendar_effects=calendar_effects,timing=physical,calendar=calendar,
       simulation_model=simulation_model)
}
