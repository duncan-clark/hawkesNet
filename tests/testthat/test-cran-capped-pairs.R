test_that('two-pair mixture matches direct enumeration and its gradient', {
  edges<-rbind(c(0L,1L),c(1L,2L),c(0L,1L),c(2L,3L))
  plan<-cran_pair_union_plan(rep(0L,4),edges[,1],edges[,2],1L,4L)
  x<-c(log(c(2,3,7,1.5)),log(c(.2,.4,.12,.3,.15)))
  evaluate<-function(x){
    rates<-exp(x[5:9]);M<-sum(rates)
    dm<-matrix(c(rep(0,4),rates),1)
    cran_pair_capped2_evaluate(plan,x[1:4],cbind(diag(4),matrix(0,4,5)),x[5:8],
      cbind(matrix(0,4,4),diag(4),0),M,dm)
  }
  direct<-function(x){
    f<-exp(x[1:4]);rates<-exp(x[5:9]);M<-sum(rates)
    a<-sum(vapply(1:4,function(i)rates[i]*prod(f[edges[i,]+1]),0.))
    b<-sum(vapply(1:4,function(i)sum(vapply(1:4,function(j)
      rates[i]*rates[j]*prod(f[unique(c(edges[i,],edges[j,]))+1]),0.)),0.))
    log(exp(-M)*(1+a)+(1-exp(-M)*(1+M))/M^2*b)
  }
  expect_equal(evaluate(x)$loglik,direct(x),tolerance=1e-12)
  h<-1e-6
  numeric<-vapply(seq_along(x),function(j){d<-rep(0,length(x));d[j]<-h;(direct(x+d)-direct(x-d))/(2*h)},0.)
  expect_equal(drop(evaluate(x)$gradient),numeric,tolerance=1e-7)
  for(scale in c(1e-8,1e-4,10)){
    z<-x;z[5:9]<-z[5:9]+log(scale)
    numerical<-vapply(seq_along(z),function(j){d<-rep(0,length(z));d[j]<-h;(evaluate(z+d)$loglik-evaluate(z-d)$loglik)/(2*h)},0.)
    expect_equal(drop(evaluate(z)$gradient),numerical,tolerance=2e-7)
  }
})

test_that('capped-two pair likelihood sums to one including repeated draws', {
  n<-4L;singles<-c(.13,.27,.09,.4)
  pairs<-rbind(c(1L,2L),c(2L,3L),c(1L,2L),c(3L,4L));rates<-c(.2,.4,.12,.31);M<-sum(rates)
  actual<-numeric(16)
  for(mask in 0:15){
    selected<-which(as.logical(intToBits(mask)[1:4]))
    keep<-which(pairs[,1] %in% selected & pairs[,2] %in% selected)
    plan<-cran_pair_union_plan(rep(0L,length(keep)),pairs[keep,1]-1L,pairs[keep,2]-1L,1L,4L)
    q<-cran_pair_capped2_evaluate(plan,-log(-expm1(-singles)),matrix(0,4,0),log(rates[keep]),matrix(0,length(keep),0),M,matrix(0,1,0))
    actual[mask+1]<-exp(-sum(singles)+sum(log(expm1(singles[selected])))+q$loglik)
  }
  expect_equal(sum(actual),1,tolerance=1e-12)
  expect_equal(actual[1],exp(-sum(singles)-M),tolerance=1e-12)
})

test_that('native capped sampler permits two proposals and preserves nonempty conditioning', {
  state<-list(names=paste0(c('u','a','b','c','x'),'@1'),birth=rep(0,5),clock=rep(0,5),
    edges=rbind(c(1L,4L),c(2L,4L),c(3L,4L)))
  p<-numeric(37);p[1:2]<-log(c(.5,.2));p[20]<--1000
  p[25:27]<-log(.6);p[29:31]<-log(1.3);p[35:37]<--1000
  release<-c(-1000,-1000,log(.2),rep(0,12));g<-cran_graph_create(state)
  set.seed(20261111);counts<-integer(3);empty<-0L
  for(i in seq_len(10000L)){
    x<-cran_graph_step(g,4,'CS-pair-migration-add-consumers-capped2',p,release,FALSE)
    if(length(x$sources)!=1L||x$sources!=0L)next
    counts[x$pair_proposals+1L]<-counts[x$pair_proposals+1L]+1L
    empty<-empty+(nrow(x$add)==0L)
  }
  M<-1.9;denominator<--expm1(-M-.5)
  expected<-c(exp(-M)*-expm1(-.5),M*exp(-M),1-exp(-M)*(1+M))/denominator
  expect_equal(empty,0L);expect_gt(sum(counts),1500L)
  expect_equal(counts/sum(counts),expected,tolerance=.04)
})

test_that('capped pairs reproduce complete bundles with overlapping migrations and types', {
  state<-list(names=c('u@1','a@1','b@2','c@1','b@1','z@1'),birth=rep(0,6),clock=rep(0,6),
    edges=rbind(c(0L,4L),c(1L,2L),c(1L,3L),c(2L,3L),c(5L,1L),c(5L,2L),c(5L,3L)))
  p<-numeric(37);p[1:2]<-log(c(.7,.4));p[20]<-log(.3)
  p[25:27]<-log(1.1);p[29:31]<--1000;p[35:37]<-log(.8)
  edges<-rbind(c(1L,2L),c(1L,3L),c(2L,3L),c(1L,4L),c(2L,4L),c(3L,4L))
  weight<-c(6,8,12,2,3,4);witness<-c(1,1,2,0,0,0)
  prob<-(1.1*weight/sum(weight)+.8*weight*witness/sum(weight*witness))/1.9
  base_candidates<-c(lapply(1:5,identity),list(c(2L,5L)))
  rates<-c(rep(.7/4,4),.4,.3);expected<-numeric(32)
  mass<-c(exp(-1.9),1.9*exp(-1.9),1-exp(-1.9)*(1+1.9))
  add<-function(edits,value){index<-1L+sum(2^(unique(edits)-1L));expected[index]<<-expected[index]+value}
  for(mask in 0:63){
    active<-as.logical(intToBits(mask)[1:6]);base<-unique(unlist(base_candidates[active]))
    w<-prod(ifelse(active,-expm1(-rates),exp(-rates)))
    add(base,w*mass[1])
    for(i in 1:6){
      add(c(base,edges[i,]),w*mass[2]*prob[i])
      for(j in 1:6)add(c(base,edges[i,],edges[j,]),w*mass[3]*prob[i]*prob[j])
    }
  }
  expect_equal(sum(expected),1,tolerance=1e-12)
  expected[1]<-0;expected<-expected/sum(expected)
  release<-c(-1000,-1000,log(.2),rep(0,12));g<-cran_graph_create(state)
  set.seed(20261112);seen<-integer(32);selected<-0L
  for(i in seq_len(18000L)){
    x<-cran_graph_step(g,4,'CS-pair-migration-add-consumers-capped2',p,release,FALSE)
    if(length(x$sources)!=1L||x$sources!=0L)next
    selected<-selected+1L
    edits<-c(match(x$add[,2],c(1L,2L,3L,5L)),if(nrow(x$delete))5L)
    seen[1L+sum(2^(edits-1L))]<-seen[1L+sum(2^(edits-1L))]+1L
  }
  expect_gt(selected,2000L);expect_equal(seen[1],0L)
  expect_true(all(abs(seen/selected-expected)<6*sqrt(expected*(1-expected)/selected)+.002))
})
