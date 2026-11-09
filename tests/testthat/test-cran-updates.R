test_that('native CRAN samplers preserve collapsed proposal and transition laws', {
  state<-list(names=c('A@1','B@1','B@2','C@1'),birth=rep(0,4),clock=rep(0,4),edges=matrix(c(0L,1L),ncol=2))
  rp<-c(-100,-100,0,rep(0,12));set.seed(9307)
  for(law in c('PA','CS-1','CS-1U','CS-pair','TERGM')){
    p<-switch(law,PA=c(log(1.2),log(.7),rep(0,6)),`CS-1`=c(log(1.2),log(.7),rep(0,19)),`CS-1U`=c(log(1.2),log(.7),rep(0,14)),`CS-pair`=c(log(1.2),log(.7),rep(0,24),log(.6)),TERGM=c(rep(-2,1),rep(0,9),.3,rep(0,6),0,0))
    g<-cran_graph_create(state);trials<-0L;add<-del<-migration<-empty<-0L
    for(i in seq_len(15000L)){
      x<-cran_graph_step(g,100,law,p,rp,apply=FALSE)
      if(!0L %in% x$sources)next
      trials<-trials+1L;a<-x$add[x$add[,1]==0,2];d<-x$delete[x$delete[,1]==0,2]
      empty<-empty+as.integer(length(a)+length(d)==0)
      add<-add+as.integer(2L %in% a);del<-del+as.integer(1L %in% d);migration<-migration+as.integer(2L %in% a&&1L %in% d)
    }
    if(law=='TERGM'){
      ra<-plogis(-2);rd<-plogis(.3);z<-1-(1-ra)^2*(1-rd);expected<-c(ra/z,rd/z,ra*rd/z)
    }else if(law=='CS-pair'){
      z<-1-exp(-2.5);ra<-1-exp(-.6);rd<-1-exp(-.7);rpair<-1-exp(-.6)
      expected<-c((1-(1-ra)*(1-rpair))/z,(1-(1-rd)*(1-rpair))/z,(rpair+(1-rpair)*ra*rd)/z)
    }else{z<-1-exp(-1.9);expected<-c((1-exp(-.6))/z,(1-exp(-.7))/z,(1-exp(-.6))*(1-exp(-.7))/z)}
    actual<-c(add,del,migration)/trials
    expect_equal(empty,0L)
    expect_true(all(abs(actual-expected)<6*sqrt(expected*(1-expected)/trials)+.002),info=law)
  }
})

test_that('directed pair changes equal full motif recomputation on four vertices', {
  statistics<-function(y)c(sum(choose(colSums(y),2))/100,sum(choose(rowSums(y),2))/100,sum(choose(colSums(y),3))/10000,sum(choose(rowSums(y),3))/10000,sum((y%*%y)*y)/10,sum(diag(y%*%y%*%y))/30)
  dyads<-which(row(matrix(0,4,4))!=col(matrix(0,4,4)),arr.ind=TRUE)
  available<-dyads[!(dyads[,1]==1&dyads[,2] %in% c(2,3)),,drop=FALSE]
  error<-0
  for(bits in 0:1023){
    y<-matrix(0,4,4);y[1,2]<-1;y[available]<-as.integer(intToBits(bits))[1:10]
    edges<-which(y==1,arr.ind=TRUE)-1L
    state<-list(names=paste0(LETTERS[1:4],'@1'),birth=rep(0,4),clock=rep(0,4),edges=edges)
    actual<-cran_directed_pair_change(cran_graph_create(state),0L,1L,2L)
    after<-y;after[1,2]<-0;after[1,3]<-1
    error<-max(error,abs(actual-(statistics(after)-statistics(y))))
  }
  expect_lt(error,1e-12)
})

test_that('directed paired sampler agrees with weighted union probabilities', {
 for(exponential in c(FALSE,TRUE)){
 y<-matrix(0,6,6);y[rbind(c(1,2),c(1,5),c(6,5),c(5,3))]<-1
 birth<-c(0,0,20,50,5,0);state<-list(names=c('A@1','B@1','B@2','B@3','C@1','D@1'),birth=birth,clock=c(100,rep(0,5)),edges=which(y==1,arr.ind=TRUE)-1L)
 statistic<-function(z)c(sum(choose(colSums(z),2))/100,sum(choose(rowSums(z),2))/100,sum(choose(colSums(z),3))/10000,sum(choose(rowSums(z),3))/10000,sum((z%*%z)*z)/10,sum(diag(z%*%z%*%z))/30)
 before<-statistic(y);di<-colSums(y)
 p<-c(log(1.3),log(.8),-1,.3,-.2,.1,.2,.7,.5,.4,1.1,0,-.5,.2,-.1,.3,.2,.4,.6,.7,.2,1,.8,-.4,.9,.3,log(.5))
 adds<-which(y[1,]==0&seq_len(6)!=1);deletes<-which(y[1,]==1)
 link<-if(exponential)exp else plogis
 wa<-vapply(adds,function(v){z<-y;z[1,v]<-1;d<-statistic(z)-before;link(sum(c(1,d,log1p(di[v]),as.integer(v %in% 3:4),0)*p[3:12]))*exp(-p[20]*(100-birth[v])/365.25)},0)
 wd<-vapply(deletes,function(v){z<-y;z[1,v]<-0;d<-statistic(z)-before;link(sum(c(1,d[1:5],-log(di[v]))*p[13:19]))*exp(-p[21]*(100-birth[v])/365.25)},0)
 wp<-vapply(3:4,function(v){z<-y;z[1,2]<-0;z[1,v]<-1;d<-statistic(z)-before;link(sum(c(1,d[c(1,3,5,6)])*p[22:26])+p[10]*log1p(di[v])-p[19]*log(di[2]))*exp(-p[20]*(100-birth[v])/365.25)},0)
 ma<-exp(p[1]+if(exponential)-.3*log(3)else 0);md<-exp(p[2]+if(exponential).7*log(3)else 0)
 ra<-ma*wa[match(3,adds)]/sum(wa);rd<-md*wd[match(2,deletes)]/sum(wd);rp<-exp(p[27])*wp/sum(wp);z<-1-exp(-ma-md-exp(p[27]))
 expected<-c(1-exp(-ra-rp[1]),1-exp(-rd-sum(rp)),1-exp(-ra-rp[1])-exp(-rd-sum(rp))+exp(-ra-rd-sum(rp)))/z
 release<-c(-100,-100,0,rep(0,12));release[10]<--1000;g<-cran_graph_create(state);set.seed(1703);n<-20000L;counts<-numeric(3)
 drawp<-if(exponential)c(p[1:2],p[c(4,6,8:12)],p[c(14,16,18,19)],p[20:21],p[23:27],.2,-.1,-.3,.7)else p
 for(i in seq_len(n)){u<-cran_graph_step(g,100,'CS-pair',drawp,release,apply=FALSE);a<-any(u$add[,1]==0&u$add[,2]==2);d<-any(u$delete[,1]==0&u$delete[,2]==1);counts<-counts+c(a,d,a&&d)}
 expect_true(all(abs(counts/n-expected)<6*sqrt(expected*(1-expected)/n)))
 }
})

test_that('joint timing rejects future history and preserves complete updates', {
 state<-list(names=c('A@1','B@1'),birth=c(0,0),clock=c(0,0),edges=matrix(integer(),0,2))
 fit<-list(law='PA',parameters=c(log(1),log(.2),rep(0,6)),release_parameters=c(-100,-100,0,rep(0,12)))
 timing<-c(mu=1,alpha_short=.2,alpha_long=.1,tau_short=.1,tau_long=1,delta=0,gamma_B=0,gamma_Delta=0)
 history<-list(list(time=5,families=0L,summary=0))
 expect_error(simulate_cran_joint(state,history,fit,timing,rep(1,20),origin=4),'Future event')
 archived<-fit;archived$law<-'CS-pair'
 expect_error(simulate_cran_joint(state,list(),archived,timing,rep(1,20),origin=4),'Incompatible edge parameters')
 a<-simulate_cran_joint(state,list(),fit,timing,rep(1,20),origin=4,horizon=1,seed=19)
 b<-simulate_cran_joint(state,list(),fit,timing,rep(1,20),origin=4,horizon=1,model='M1',seed=19)
 expect_equal(a,b)
 expect_true(all(a$wait>=0&a$wait<=1))
})

test_that('native summaries preserve directed betweenness and undirected motifs', {
  state<-list(names=c('A@1','B@1','C@1'),birth=rep(0,3),clock=rep(0,3),edges=rbind(c(0L,1L),c(1L,2L)))
  g<-cran_graph_create(state);s<-cran_graph_summary(g)$metrics
  expect_equal(s[['max_betweenness']],.5)
  expect_equal(s[['mean_betweenness']],1/6)
  expect_equal(s[['weak_giant_fraction']],1)
  expect_equal(s[['triangles']],0)
  u<-list(time=1,births=character(),clocks=numeric(),sources='A@1',add=matrix(c('A@1','C@1'),ncol=2),delete=matrix(character(),0,2))
  cran_graph_replay(g,list(u));s<-cran_graph_summary(g)$metrics
  expect_equal(s[['triangles']],1)
  expect_equal(s[['transitivity']],1)
  expect_equal(s[['max_betweenness']],0)
})

test_that('calendar timing likelihood agrees with direct Hawkes integration', {
 d<-list(times=c(1,2),initial=0,tz=c(0,0),tw=c(1,1),tx=c(0,0),tu=c(0L,0L),tp=c(0L,1L,2L),cz=numeric(),cu=integer(),cp=c(0L,0L,0L),calendar=rep(1,4))
 physical<-c(2,.2,.3,.1,1.5);p<-log(physical)
 intensity<-c(2,2+sum(physical[2:3]/physical[4:5]*exp(-1/physical[4:5])))
 integral<-6+sum(vapply(c(1,2),function(t)sum(physical[2:3]*(1-exp(-(3-t)/physical[4:5]))),0))
 result<-loglik_cran_timing(d,p,'M0',start=0,end=3)
 expect_equal(result$loglik,sum(log(intensity))-integral,tolerance=1e-12)
 expect_equal(result$compensator,integral,tolerance=1e-12)
 expect_equal(loglik_cran_timing(d,c(p,0,0,0),'M2',start=0,end=3)$loglik,result$loglik,tolerance=1e-12)
})

test_that('release fitting uses supplied parameters rather than a stale cache', {
 m<-new.env(parent=emptyenv());m$n<-1L;m$year<-2020L;m$families<-1L
 m$release_data<-list(initial=matrix(c(0,0,0,0,0,0),1,6),changes=matrix(numeric(),0,6),change_ptr=c(0,0),selected=matrix(c(0,1,1),1,3),selected_ptr=c(0,1),new_count=0,time=1)
 m$release_values<-list(train=list(loglik=NaN,gradient=matrix(NA_real_,1,15)))
 result<-fit_cran_release(m,parameters=c(rep(log(.1),3),rep(0,12)),control=list(maxit=0))
 mean<-exp(result$parameters[1:3]);expected<--sum(mean[1:2])+log(-expm1(-mean[3]))-log(-expm1(-sum(mean)))
 expect_true(is.finite(result$loglik))
 expect_equal(result$loglik,expected,tolerance=1e-12)
})

test_that('proposal means distinguish new packages from new majors before edits', {
 state<-list(names=paste0(LETTERS[1:4],'@1'),birth=rep(0,4),clock=rep(0,4),edges=matrix(integer(),0,2))
 set.seed(3105)
 for(law in c('PA','CS-pair'))for(kind in c('new','major')){
  base<-if(law=='PA')8L else 20L;p<-numeric(base+4L);p[1:2]<-log(c(1.2,.7));if(law=='CS-pair')p[20]<-log(.2)
  p[base+1:4]<-c(log(2),log(3),.2,.8)
  release<-c(rep(-100,3),rep(0,12));release[if(kind=='new')1L else 2L]<-log(.3)
  residual<-variance<-numeric(1500L)
  for(i in seq_along(residual)){
   u<-cran_graph_step(cran_graph_create(state),100,law,p,release,apply=FALSE)
   b<-length(u$born);j<-4L+b-1L;r<--expm1(-1.2*(if(kind=='new')2 else 3)/j)
   residual[i]<-nrow(u$add)-b*j*r;variance[i]<-b*j*r*(1-r)
  }
  expect_true(abs(sum(residual))<6*sqrt(sum(variance)),info=paste(law,kind))
 }
})

test_that('large covariate-dependent means use the same collapsed proposal law', {
 state<-list(names=paste0('P',1:30,'@1'),birth=rep(0,30),clock=c(100,rep(0,29)),edges=cbind(0L,1:20))
 p<-c(4,4,rep(0,6),0,0,3,3);release<-c(-100,-100,0,rep(0,12));release[10]<--1000
 u<-cran_graph_step(cran_graph_create(state),100,'PA',p,release,apply=FALSE)
 expect_equal(nrow(u$add),9L);expect_equal(nrow(u$delete),20L)
})

test_that("paired additions sum correctly with overlapping migrations and singles", {
  # Three candidate additions and one candidate deletion. Two possible
  # triangle pairs overlap each other; one overlaps the migration.
  rho <- c(.11,.21,.31,.41,.51)
  r <- -expm1(-rho); mt <- .63; kappa <- 1.12
  z <- c(2,3*exp(kappa)); z <- z/sum(z)
  codes <- c(1L,2L,4L,8L,9L)
  base <- numeric(16); direct <- numeric(16)
  for(mask in 0:31) {
    take <- as.logical(intToBits(mask)[1:5])
    p <- prod(ifelse(take,r,1-r))
    s <- Reduce(bitwOr,codes[take],init=0L)
    base[s+1L] <- base[s+1L]+p
    for(j in 1:3) {
      new <- bitwOr(s,c(0L,3L,6L)[j])
      direct[new+1L] <- direct[new+1L]+p*c(exp(-mt),-expm1(-mt)*z)[j]
    }
  }
  expect_equal(sum(direct),1,tolerance=1e-14)
  for(s in 0:15) {
    a <- which(as.logical(intToBits(s)[1:3]))
    present <- bitwAnd(s,9L)==9L
    has_t <- which(bitwAnd(s,c(3L,6L))==c(3L,6L))
    m <- list(ng=1L,meanT=1L,copyT=2L,Tclass=1L,ai=if(present)match(1L,a)else integer(),
      ag=rep(0,length(a)),triangle=list(group=rep(0,length(has_t)),
        mass=matrix(c(5,3),1L),a=match(c(1L,2L)[has_t],a),b=match(c(2L,3L)[has_t],a),
        log_weight=log(c(2,3)[has_t]),copied=c(0,1)[has_t]))
    u <- if(present)log(expm1(rho[5])/(r[1]*r[4]))else numeric()
    val <- .cran_triangle_correction(m,c(log(mt),kappa),log(rho[a]),matrix(0,length(a),2),
                                     u,matrix(0,length(u),2))
    expect_equal(base[s+1L]*exp(val$loglik),direct[s+1L],tolerance=1e-13)
  }
})

test_that("native sampler draws paired additions and migrations from their joint union law", {
  state <- list(names=c("u@1","b@1","b@2","c@1","d@1"),birth=rep(0,5),clock=rep(0,5),
                edges=rbind(c(0L,1L),c(2L,3L),c(3L,4L)))
  g <- cran_graph_create(state)
  p <- numeric(28)
  p[c(1,2,20,25:27)] <- log(c(.9,.4,.5,.6,.6,.6))
  rp <- c(-1000,-1000,log(.4),rep(0,12))
  rates <- c(rep(.3,3),.4,.5); r <- -expm1(-rates)
  mass <- numeric(16)
  for(i in 0:31) {
    yes <- as.logical(intToBits(i)[1:5])
    base <- Reduce(bitwOr,c(1L,2L,4L,8L,9L)[yes],init=0L)
    pr <- prod(ifelse(yes,r,1-r))
    for(j in 1:3) {
      s <- bitwOr(base,c(0L,3L,6L)[j])
      mass[s+1L] <- mass[s+1L]+pr*c(exp(-.6),-expm1(-.6)*c(1/3,2/3))[j]
    }
  }
  mass[1] <- 0;mass <- mass/sum(mass)
  set.seed(22919)
  counts <- integer(16)
  for(i in seq_len(12000L)){
    a <- cran_graph_step(g,4,"CS-pair-triangle",p,rp,FALSE)
    if(length(a$sources)==1L&&a$sources==0L){
      mask <- sum(2^(a$add[,2]-2L))+8L*nrow(a$delete)
      counts[mask+1L] <- counts[mask+1L]+1L
    }
  }
  expect_gt(sum(counts),1500)
  expect_lt(max(abs(counts/sum(counts)-mass)),.03)
})

test_that("triangle and diamond classes compete for a single paired addition", {
  r <- -expm1(-c(.12,.23,.34))
  mt <- .4; md <- .7; active <- -expm1(-mt-md)
  # The two optional pairs share addition 2, but cannot both be proposed.
  direct <- numeric(8)
  for(s in 0:7){
    edits<-as.logical(intToBits(s)[1:3]); q<-prod(ifelse(edits,r,1-r))
    for(k in 1:3){
      union <- bitwOr(s,c(0L,3L,6L)[k])
      direct[union+1L]<-direct[union+1L]+q*c(exp(-mt-md),active*c(mt,md)/(mt+md))[k]
    }
  }
  expect_equal(sum(direct),1,tolerance=1e-14)
  for(s in 0:7){
    included <- as.logical(intToBits(s)[1:3]); ix<-which(included)
    a <- if(bitwAnd(s,3L)==3L)c(1L,2L)else integer()
    d <- if(bitwAnd(s,6L)==6L)c(2L,3L)else integer()
    m<-list(ng=1L,meanT=1L,copyT=2L,meanDiamond=3L,Tclass=1L,ai=integer(),ag=rep(0,length(ix)),
      triangle=list(group=rep(0,length(a)/2),a=match(a[1L],ix,nomatch=0L)[seq_len(length(a)/2)],
        b=match(a[2L],ix,nomatch=0L)[seq_len(length(a)/2)],log_weight=rep(log(5),length(a)/2),
        copied=rep(0,length(a)/2),mass=matrix(c(5,0),1L),
        diamond=list(group=rep(0,length(d)/2),a=match(d[1L],ix,nomatch=0L)[seq_len(length(d)/2)],
          b=match(d[2L],ix,nomatch=0L)[seq_len(length(d)/2)],log_weight=rep(log(6),length(d)/2),mass=6)))
    q0<-prod(ifelse(included,r,1-r))
    val<-.cran_triangle_correction(m,c(log(mt),0,log(md)),log(-log1p(-r[ix])),
                                   matrix(0,length(ix),3),numeric(),matrix(0,0,3))
    expect_equal(q0*exp(val$loglik),direct[s+1L],tolerance=1e-13)
  }
})

test_that("diamond sampler weights unlinked pairs by shared outgoing dependency counts", {
  # b,c,d all depend on z, and b,d also depend on z2. The shared edges
  # from b and c back to the releaser u do not count as a fourth vertex.
  state<-list(names=c("u@1","b@1","c@1","d@1","z@1","z2@1"),birth=rep(0,6),clock=rep(0,6),
              edges=rbind(c(1L,4L),c(2L,4L),c(3L,4L),c(1L,5L),c(3L,5L),c(1L,0L),c(2L,0L)))
  p<-numeric(31);p[c(1L,2L)]<--28;p[25:27]<--1e6;p[29:31]<-log(.6)
  rp<-c(-1000,-1000,log(.4),rep(0,12))
  g<-cran_graph_create(state); counts<-integer(3);set.seed(3605)
  for(i in seq_len(12000L)){
    a<-cran_graph_step(g,4,"CS-pair-diamond",p,rp,FALSE)
    if(length(a$sources)==1L&&a$sources==0L){
      # Sums of IDs distinguish b,c (3), b,d (4), c,d (5).
      stopifnot(nrow(a$add)==2L)
      j<-sum(a$add[,2])-2L
      counts[j]<-counts[j]+1L
    }
  }
  expect_gt(sum(counts),1300L)
  expect_equal(counts/sum(counts),c(1,2,1)/4,tolerance=.055)
})

test_that("the combined paired law retains the migration channel", {
  state <- list(names=c("u@1","u@2","b@1","b@2","c@1","d@1"),birth=rep(0,6),
                clock=rep(0,6),edges=rbind(c(0L,2L),c(0L,4L),c(2L,4L),c(3L,5L)))
  rp <- c(-1000,-1000,log(.4),rep(0,12))
  p <- numeric(31)
  p[c(1,2,25:27,29:31)] <- log(c(.5,.5,1,1,1,1,1,1))
  draw <- function(law,parameters) {
    graph <- cran_graph_create(state)
    set.seed(20261106)
    replicate(100,cran_graph_step(graph,4,law,parameters,rp,FALSE),simplify=FALSE)
  }
  combined <- draw("CS-pair-migration-add",p)
  expect_true(any(vapply(combined,function(x)x$pair_proposals>0L,FALSE)))
  expect_identical(combined,draw("CS-pair-diamond",p))
  baseline <- draw("CS-pair-add",p)
  expect_true(all(vapply(baseline,function(x)x$pair_proposals==0L,FALSE)))
  p[20] <- 10
  expect_identical(draw("CS-pair-add",p),baseline)
})
