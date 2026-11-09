# Two normalized component laws share structural preferences. The latent
# component is drawn once for the complete update, so release and edge rates
# are fitted jointly and the observed component label is never imputed.
.cran_mixture_layout <- function(m) {
  edge_count <- length(m$names)
  release_names <- c("log_new_package_proposals", "log_major_proposals", "log_update_proposals",
    paste0(rep(c("major", "update"), each=6), "_", rep(c("age", "age_squared", "in1", "out1", "in2", "out2"), 2)))
  edge_shifts <- list(mixture_add_shift=if(m$law %in% c("TERGM", "TERGM-logdegree"))1L else m$rateA[1],
    mixture_delete_shift=if(m$law %in% c("TERGM", "TERGM-logdegree"))11L else m$rateD[1])
  if(!is.null(m$meanP) && m$lower[m$meanP]<m$upper[m$meanP])
    edge_shifts$mixture_migration_shift <- m$meanP
  pair_rates <- c(m$meanT,m$meanDiamond,m$meanConsumer)
  if(length(pair_rates))edge_shifts$mixture_pair_shift <- pair_rates
  shift_names <- c("mixture_new_shift", "mixture_major_shift", "mixture_update_shift", names(edge_shifts), "mixture_logit")
  joint_names <- c(m$names, paste0("release_",release_names), shift_names)
  # Component 2 has at least as many new-package proposals, fixing label order.
  lower <- c(m$lower, rep(-15,3), rep(-Inf,12), 0, rep(-4,length(shift_names)-2), -6)
  upper <- c(m$upper, rep(5,3), rep(Inf,12), rep(4,length(shift_names)-1), 6)
  list(names=joint_names, lower=lower, upper=upper, edge=seq_len(edge_count),
    release=edge_count+seq_len(15L), release_names=release_names,
    shift=setNames(edge_count+15L+seq_along(shift_names),shift_names), edge_shifts=edge_shifts)
}

.cran_mixture_components <- function(m,parameters,layout=.cran_mixture_layout(m)) {
  p <- unname(.cran_parameters(parameters,layout$names,"joint mixture parameters"))
  edge1 <- edge2 <- p[layout$edge]
  release1 <- release2 <- p[layout$release]
  release2[1:3] <- release2[1:3]+p[layout$shift[c("mixture_new_shift","mixture_major_shift","mixture_update_shift")]]
  for(name in names(layout$edge_shifts))edge2[layout$edge_shifts[[name]]] <-
    edge2[layout$edge_shifts[[name]]]+p[layout$shift[name]]
  list(edge=list(edge1,edge2),release=list(release1,release2),logit=p[layout$shift["mixture_logit"]])
}

.cran_logmix <- function(first,second,logit) {
  log_weight2 <- stats::plogis(logit,log.p=TRUE)
  log_weight1 <- stats::plogis(-logit,log.p=TRUE)
  a <- first+log_weight1
  b <- second+log_weight2
  largest <- pmax(a,b)
  value <- largest+log(exp(a-largest)+exp(b-largest))
  list(loglik=value,responsibility=exp(b-value))
}

.cran_mixture_release <- function(m,parts) {
  prefix <- if(is.null(m$mixture_upto))0L else m$mixture_upto
  coefficients <- as.numeric(parts$release[[1]][4:15])
  cached <- m$mixture_release_cache
  if(is.null(cached)||!identical(cached$coefficients,coefficients)||cached$prefix!=prefix){
    prepared <- cran_release_mixture_cache(m$release_data,parts$release[[1]],m$families,prefix)
    cached <- list(coefficients=coefficients,prefix=prefix,data=prepared)
    m$mixture_release_cache <- cached
  }
  cran_release_mixture_cached(cached$data,parts$release[[1]],parts$release[[2]])
}

.cran_mixture_loglik <- function(m,parameters,per_update=FALSE) {
  if(inherits(parameters,"cran_update_fit"))parameters <- parameters$mixture_parameters
  layout <- .cran_mixture_layout(m)
  parts <- .cran_mixture_components(m,parameters,layout)
  edges <- lapply(parts$edge,function(p).cran_single_loglik(m,p,TRUE))
  releases <- .cran_mixture_release(m,parts)
  mixture <- .cran_logmix(edges[[1]]$loglik+releases[[1]]$loglik,
    edges[[2]]$loglik+releases[[2]]$loglik,parts$logit)
  weight2 <- mixture$responsibility
  weight1 <- 1-weight2
  gradient <- matrix(0,m$n,length(layout$names),dimnames=list(NULL,layout$names))
  gradient[,layout$edge] <- weight1*edges[[1]]$gradient+weight2*edges[[2]]$gradient
  gradient[,layout$release] <- weight1*releases[[1]]$gradient+weight2*releases[[2]]$gradient
  gradient[,layout$shift[c("mixture_new_shift","mixture_major_shift","mixture_update_shift")]] <- weight2*releases[[2]]$gradient[,1:3,drop=FALSE]
  for(name in names(layout$edge_shifts))gradient[,layout$shift[name]] <-
    weight2*rowSums(edges[[2]]$gradient[,layout$edge_shifts[[name]],drop=FALSE])
  gradient[,layout$shift["mixture_logit"]] <- weight2-stats::plogis(parts$logit)
  release_marginal <- .cran_logmix(releases[[1]]$loglik,releases[[2]]$loglik,parts$logit)$loglik
  if(per_update)return(list(loglik=mixture$loglik,gradient=gradient,release=release_marginal,
    edge=mixture$loglik-release_marginal,responsibility=weight2))
  list(loglik=sum(mixture$loglik),gradient=colSums(gradient))
}

.cran_mixture_start <- function(m,window,previous=NULL,batch_probability=.25,rate_ratio=4) {
  layout <- .cran_mixture_layout(m)
  if(inherits(previous,"cran_update_fit")&&!is.null(previous$mixture_parameters))
    return(.cran_parameters(previous$mixture_parameters,layout$names,"mixture warm start"))
  if(is.numeric(previous)&&length(previous)==length(layout$names))
    return(.cran_parameters(previous,layout$names,"mixture warm start"))
  edge <- if(inherits(previous,"cran_update_fit"))previous$parameters else
    if(is.numeric(previous))previous else .cran_old_start(m,window)
  release <- if(inherits(previous,"cran_update_fit"))previous$release_parameters else .cran_release_parameters(m,window)
  if(length(edge)!=length(m$names)||length(release)!=15L)stop("Mixture warm start must supply the base edge and release parameters")
  p <- setNames(numeric(length(layout$names)),layout$names)
  p[layout$edge] <- edge
  p[layout$release] <- release
  p[layout$shift[c("mixture_new_shift","mixture_major_shift","mixture_update_shift")]] <- log(rate_ratio)
  p[layout$release[1:3]] <- release[1:3]-log(1+batch_probability*(rate_ratio-1))
  for(name in names(layout$edge_shifts)){
    ratio <- if(name=="mixture_add_shift")2.5 else 1.5
    p[layout$shift[name]] <- log(ratio)
    index <- layout$edge_shifts[[name]]
    p[index] <- p[index]-log(1+batch_probability*(ratio-1))
  }
  p[layout$shift["mixture_logit"]] <- stats::qlogis(batch_probability)
  pmin(layout$upper,pmax(layout$lower,p))
}

.cran_fit_mixture <- function(m,window,starts,control,parameters,checkpoint) {
  if(identical(m$pair_proposals,"poisson"))cran_pair_union_check(m$pair_union_plan,m$pair_max_states)
  layout <- .cran_mixture_layout(m)
  take <- if(window=="train")m$year<=2023 else rep(TRUE,m$n)
  count <- sum(take)
  m$mixture_upto <- max(which(take))
  on.exit(m$mixture_upto <- NULL,add=TRUE)
  supplied <- is.list(parameters)&&!inherits(parameters,"cran_update_fit")
  initial <- .cran_mixture_start(m,window,if(supplied)parameters[[1]]else parameters)
  initial <- pmin(layout$upper,pmax(layout$lower,initial))
  tolerance <- if(is.null(control$polish_tol))2e-6 else control$polish_tol
  control$polish_tol <- NULL
  if(!is.finite(tolerance)||tolerance<=0||tolerance>2e-5)stop("Invalid mixture fitting tolerance")
  last <- NULL
  result <- NULL
  calls <- 0L
  evaluate <- function(p){
    if(!identical(as.numeric(p),last)){
      raw <- .cran_mixture_loglik(m,p,TRUE)
      result <<- list(value=-mean(raw$loglik[take]),gradient=-colMeans(raw$gradient[take,,drop=FALSE]))
      last <<- as.numeric(p)
      calls <<- calls+1L
      if(calls%%50L==0L)message(m$law," ",window," mixture evaluation ",calls,": ",format(-result$value*count,digits=12))
    }
    result
  }
  objective <- function(p)evaluate(p)$value
  gradient <- function(p)evaluate(p)$gradient
  raw <- .cran_mixture_loglik(m,initial,TRUE)
  scales <- pmin(100,pmax(.001,1/sqrt(colMeans(raw$gradient[take,,drop=FALSE]^2)+1e-10)))
  projected <- function(p)max(abs(.cran_projected(p,gradient(p),layout$lower,layout$upper)))
  release_preferences <- layout$release[4:15]
  optimize_rates <- function(start) {
    lower <- layout$lower
    upper <- layout$upper
    lower[release_preferences] <- upper[release_preferences] <- start[release_preferences]
    fit <- stats::optim(start,objective,gradient,method="L-BFGS-B",lower=lower,upper=upper,
      control=utils::modifyList(list(maxit=700L,factr=100,pgtol=1e-8,lmm=30L,parscale=scales),control))
    block_score <- max(abs(.cran_projected(fit$par,gradient(fit$par),lower,upper)))
    if(block_score>tolerance){
      refined <- stats::nlminb(fit$par,objective,gradient,lower=lower,upper=upper,scale=1/scales,
        control=list(eval.max=1200L,iter.max=600L,rel.tol=1e-13,x.tol=1e-9))
      if(refined$objective<=objective(fit$par)+1e-10)fit$par <- refined$par
    }
    fit
  }
  improve_preferences <- function(p,exact_hessian=FALSE) {
    raw <- .cran_mixture_loglik(m,p,TRUE)
    score <- -colMeans(raw$gradient[take,release_preferences,drop=FALSE])
    if(max(abs(score))<=tolerance)return(p)
    if(exact_hessian){
      information <- matrix(0,12L,12L)
      for(j in seq_along(release_preferences)){
        coordinate <- release_preferences[j]
        step <- 1e-5*max(1,abs(p[coordinate]))
        delta <- numeric(length(p));delta[coordinate] <- step
        information[,j] <- (gradient(p+delta)-gradient(p-delta))[release_preferences]/(2*step)
      }
      information <- (information+t(information))/2
    }else{
      information <- crossprod(raw$gradient[take,release_preferences,drop=FALSE])/count
    }
    scaling <- 1/sqrt(pmax(diag(information),1e-12))
    decomposition <- eigen(information*outer(scaling,scaling),symmetric=TRUE)
    direction <- -scaling*drop(decomposition$vectors %*%
      (drop(crossprod(decomposition$vectors,score*scaling))/pmax(decomposition$values,1e-6)))
    previous <- objective(p)
    step <- 1
    while(step>=1/128){
      proposal <- p
      proposal[release_preferences] <- p[release_preferences]+step*direction
      if(objective(proposal)<=previous+1e-4*step*sum(score*direction))return(proposal)
      step <- step/2
    }
    p
  }
  trials <- vector("list",starts)
  for(i in seq_len(starts)){
    start <- if(supplied&&i<=length(parameters)).cran_mixture_start(m,window,parameters[[i]])else initial
    if(i>1L&&!supplied){
      start <- .cran_mixture_start(m,window,parameters,batch_probability=.4,rate_ratio=2.5)
      if((inherits(parameters,"cran_update_fit")&&!is.null(parameters$mixture_parameters))||
         (is.numeric(parameters)&&length(parameters)==length(layout$names)))
        start <- initial+.05*scales*sin(seq_along(initial)*i)
    }
    start <- pmin(layout$upper,pmax(layout$lower,start))
    fit <- list(par=start,convergence=0L)
    for(cycle in seq_len(10L)){
      fit <- optimize_rates(fit$par)
      if(!is.null(checkpoint))saveRDS(setNames(fit$par,layout$names),checkpoint)
      if(projected(fit$par)<=tolerance)break
      updated <- improve_preferences(fit$par,exact_hessian=cycle>=7L)
      if(identical(updated,fit$par)&&cycle<7L)updated <- improve_preferences(fit$par,exact_hessian=TRUE)
      fit$par <- updated
      message(m$law," ",window," mixture block ",cycle,"; joint projected score ",format(projected(fit$par),digits=3))
    }
    fit$value <- objective(fit$par)
    fit$projected_gradient <- projected(fit$par)
    trials[[i]] <- fit
    if(!is.null(checkpoint))saveRDS(setNames(fit$par,layout$names),checkpoint)
    message(m$law," ",window," mixture start ",i,": ",format(-fit$value*count,digits=13),
      "; projected gradient ",format(fit$projected_gradient,digits=3))
  }
  best <- trials[[which.min(vapply(trials,`[[`,0.,"value"))]]
  if(best$projected_gradient>2e-5)stop("Mixture projected gradient remains above tolerance; checkpoint retained")
  parameters <- setNames(best$par,layout$names)
  parts <- .cran_mixture_components(m,parameters,layout)
  m$mixture_upto <- NULL
  raw <- .cran_mixture_loglik(m,parameters,TRUE)
  native_shifts <- setNames(numeric(4),c("mixture_add_shift","mixture_delete_shift","mixture_migration_shift","mixture_pair_shift"))
  native_shifts[names(layout$edge_shifts)] <- parameters[layout$shift[names(layout$edge_shifts)]]
  native_release_shifts <- parameters[layout$shift[c("mixture_new_shift","mixture_major_shift","mixture_update_shift","mixture_logit")]]
  prediction <- lapply(c(2024L,2025L),function(year){
    use <- m$year==year
    list(year=year,updates=sum(use),edge_loglik=sum(raw$edge[use]),release_loglik=sum(raw$release[use]),update_loglik=sum(raw$loglik[use]))
  })
  names(prediction) <- c("2024","2025")
  bound <- abs(parameters-layout$lower)<1e-5 | abs(parameters-layout$upper)<1e-5
  structure(list(law=m$law,window=window,
    parameters=c(setNames(parts$edge[[1]],m$names),native_shifts),
    release_parameters=c(setNames(parts$release[[1]],layout$release_names),native_release_shifts),
    mixture_parameters=parameters,mixture="two_component",mixture_probability=stats::plogis(parts$logit),
    component_parameters=list(edge=lapply(parts$edge,setNames,m$names),release=lapply(parts$release,setNames,layout$release_names)),
    parameter_count=sum(layout$lower<layout$upper),edge_loglik=sum(raw$edge[take]),update_loglik=sum(raw$loglik[take]),
    predictive=prediction,standard_error=setNames(rep(NA_real_,length(parameters)),layout$names),
    bound_parameters=layout$names[bound & layout$lower<layout$upper],fixed_parameters=layout$names[layout$lower==layout$upper],
    projected_gradient=best$projected_gradient,multistart_loglik=-count*vapply(trials,`[[`,0.,"value"),
    convergence=vapply(trials,`[[`,0L,"convergence"),contributions=list(year=m$year,edge=raw$edge,release=raw$release,responsibility=raw$responsibility),
    weights=if(m$law %in% c("TERGM","TERGM-logdegree"))"logit"else if(m$pa)"exponential"else m$weights,
    proposal_means=m$proposal_means,term_specification=m$term_specification,migration_degree=isTRUE(m$migration_degree),
    pair_proposals=m$pair_proposals,pair_max_states=m$pair_max_states,
    shared_people=!is.null(m$shared_people),
    specification="Two-component complete-update mixture; shared structural preferences; component-specific release and edit rates; exact likelihood"),
    class="cran_update_fit")
}
