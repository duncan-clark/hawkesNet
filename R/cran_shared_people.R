# Keep archived zero slots so the likelihood and sampler share one layout.
.cran_fix_motif_terms <- function(m) {
  index <- grep("star|transitive|diamond|cyclic|shared_consumer", m$names)
  m$lower[index] <- m$upper[index] <- 0
  invisible(m)
}

.cran_validate_people <- function(people) {
  integer_values <- function(x, positive=FALSE) is.numeric(x) &&
    all(is.finite(x)) && all(x==floor(x)) && all(abs(x)<=.Machine$integer.max) &&
    (!positive || all(x>0))
  if(!is.list(people) || !all(c("bin","package","people") %in% names(people)) ||
     !integer_values(people$bin) || !is.character(people$package) ||
     anyNA(people$package) || any(!nzchar(people$package)) || !is.list(people$people) ||
     length(people$bin)!=length(people$package) || length(people$bin)!=length(people$people) ||
     is.unsorted(people$bin) || anyDuplicated(data.frame(bin=people$bin,package=people$package)) ||
     !all(vapply(people$people,integer_values,logical(1),positive=TRUE)))
    stop("people must contain sorted integer bins, package names and positive integer ID lists, with one row per package/bin")
  people$bin <- as.integer(people$bin)
  people$people <- lapply(people$people,function(x)unique(as.integer(x)))
  people
}

.cran_attach_people <- function(m, people, sources) {
  tergm <- m$law=="TERGM-logdegree"
  stopifnot(tergm || (m$law %in% c("CS-1","CS-pair") && m$exponential &&
                      m$proposal_means=="covariate"),
            m$term_specification %in% c("shared_consumer", "motif_free"),
            m$mixture == "two_component", !m$migration_degree, is.null(m$triangle))
  people <- .cran_validate_people(people)
  fixed <- c("addition_transitive","addition_diamond","addition_shared_consumer","deletion_transitive",
             "migration","pair_transitive","pair_diamond","pair_shared_consumer")
  if(m$law!="CS-pair")fixed <- setdiff(fixed,"migration")
  fixed <- match(intersect(fixed,m$names),m$names)
  m$lower[fixed] <- m$upper[fixed] <- 0
  a <- m$a
  data <- cran_graph_data(m$path,"2020-01-01")
  records <- Filter(function(r)r$time>=data$origin,data$records)
  stopifnot(length(records)==m$n)
  cache <- cran_shared_people_cache(data$state,records,sources,as.integer(m$gb),people)
  add <- cache$observed_add
  deletion <- cache$deletion
  pairs <- cache$pairs
  pop <- match("popularity_add",m$namesA)
  migration <- match("migration",m$namesA)
  copy <- match("copy",m$namesA)
  # Compare identities in the replay order used by the pattern and migration caches.
  stopifnot(isTRUE(all.equal(unname(add[,1:5]),
    unname(cbind(m$ag,m$XA[,pop],a$added_birth,m$XA[,migration],m$XA[,copy])),tolerance=1e-12)))
  # Identical baseline rows may tie. Use the actual paired deletion's replay
  # index below so that the tied rows retain the correct joint event.
  old <- cbind(m$eg,drop(a$existing_degree),a$existing_birth,a$existing_outcome)
  order_deletions <- function(z)order(z[,1],z[,2],round(z[,3]*365.25*144),z[,4])
  old_order <- order_deletions(old)
  replay_order <- order_deletions(deletion)
  stopifnot(isTRUE(all.equal(unname(old[old_order,]),unname(deletion[replay_order,1:4]),tolerance=1e-12)))
  to_old <- integer(nrow(deletion))
  to_old[replay_order] <- old_order
  deleted_people <- numeric(nrow(deletion))
  deleted_people[to_old] <- deletion[,5]
  index <- length(m$names)+1:3
  if(tergm){
    shared <- cache$shared_add
    g <- shared[,1]+1L
    d <- which(deleted_people==1)
    dg <- m$eg[d]+1L
    age <- function(group,birth)birth-m$time[m$gb[group]+1L]
    m$shared_people <- list(index=index, addition_group=g-1L, deletion_group=dg-1L,
      addition=cbind(1,shared[,2],shared[,4:5,drop=FALSE],age(g,shared[,3]),
                     m$ZA[g,-1,drop=FALSE]),
      deletion=cbind(1,-log(pmax(drop(a$existing_degree)[d],1)),age(dg,a$existing_birth[d]),
                     m$ZD[dg,2]),
      addition_index=c(1L,8:10,18L,20:22),deletion_index=c(11L,17L,19L,23L),
      pair=pairs[,5],pair_group=pairs[,1],observed_pairs=pairs[,6]>=0)
    m$obs <- cbind(m$obs,.cran_sum(m$ag,add[,6],m$ng),
                   .cran_sum(m$sg,deleted_people[m$sel],m$ng),
                   .cran_sum(m$ag,add[,4]*add[,6],m$ng))
  } else {
    if(m$law=="CS-pair"){
      stopifnot(isTRUE(all.equal(unname(pairs[,1:4]),unname(cbind(m$pg,m$Pb,m$Ppop)),tolerance=1e-12)),
                identical(as.integer(which(pairs[,6]>=0)),as.integer(m$oi)),
                identical(as.integer(pairs[m$oi,7]+1L),as.integer(m$ai)))
      m$di <- nrow(add)+match(to_old[pairs[m$oi,6]+1L],m$sel)
      stopifnot(!anyNA(m$di),!anyDuplicated(m$di))
    }
    shared <- cache$shared_add
    base <- matrix(0,nrow(shared),m$na)
    base[,c(pop,migration,copy)] <- shared[,c(2,4,5)]
    additive <- function(mig, flag)
      if(m$law=="CS-1") cbind(flag,flag*mig) else matrix(flag,ncol=1L)
    m$shared_people <- list(index=index, group=shared[,1], base=base,
      birth=shared[,3], design=additive(shared[,4],rep(1,nrow(shared))),
      observed_add=additive(add[,4],add[,6]), deletion=deleted_people, pair=pairs[,5],
      pair_group=pairs[,1], observed_pairs=pairs[,6]>=0)
  }
  m$names <- c(m$names,"addition_shared_people","deletion_shared_people","migration_shared_people")
  m$lower <- c(m$lower,rep(-8,3L))
  m$upper <- c(m$upper,rep(8,3L))
  invisible(m)
}

# Modify the Bernoulli normalizer only at candidates with shared people.
.cran_people_tergm_sum <- function(m,H,dh,p) {
  x <- m$shared_people
  eta_a <- drop(x$addition%*%p[x$addition_index])
  eta_d <- drop(x$deletion%*%p[x$deletion_index])
  eta_a_shared <- eta_a+p[x$index[1]]+x$addition[,3]*p[x$index[3]]
  eta_d_shared <- eta_d+p[x$index[2]]
  H <- H+.cran_sum(x$addition_group,.cran_softplus(eta_a_shared)-.cran_softplus(eta_a),m$ng)+
    .cran_sum(x$deletion_group,.cran_softplus(eta_d_shared)-.cran_softplus(eta_d),m$ng)
  pa <- stats::plogis(eta_a_shared)
  pd <- stats::plogis(eta_d_shared)
  dh[,x$addition_index] <- dh[,x$addition_index,drop=FALSE]+
    .cran_sum(x$addition_group,x$addition*(pa-stats::plogis(eta_a)),m$ng)
  dh[,x$deletion_index] <- dh[,x$deletion_index,drop=FALSE]+
    .cran_sum(x$deletion_group,x$deletion*(pd-stats::plogis(eta_d)),m$ng)
  list(H=H,dh=cbind(dh,.cran_sum(x$addition_group,pa,m$ng),
    .cran_sum(x$deletion_group,pd,m$ng),.cran_sum(x$addition_group,pa*x$addition[,3],m$ng)))
}

# Tilt only the shared-person subset of the full addition normalizer.
.cran_people_add_sum <- function(m, normalizer, p, index) {
  x <- m$shared_people
  weight <- exp(drop(x$base%*%p[m$ia])+p[m$ageA]*x$birth-normalizer$logsum[x$group+1L])
  mass <- .cran_sum(x$group,weight,m$ng)
  if(any(mass>1+1e-8))stop("Shared-person subset exceeds complete addition normalizer")
  tilted <- weight*exp(drop(x$design%*%p[index]))
  denominator <- pmax(0,1-mass)+.cran_sum(x$group,tilted,m$ng)
  change <- tilted-weight
  score <- normalizer$score+.cran_sum(x$group,x$base*change,m$ng)
  score <- cbind(score,.cran_sum(x$group,x$design*tilted,m$ng))
  list(logsum=normalizer$logsum+log(denominator),score=score/denominator,
       birth=(normalizer$birth+.cran_sum(x$group,change*x$birth,m$ng))/denominator)
}
