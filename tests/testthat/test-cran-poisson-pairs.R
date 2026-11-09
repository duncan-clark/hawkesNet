# Independent finite enumeration of Poisson proposal-presence indicators.
pair_union_brute <- function(candidates,rates,n) {
  out<-numeric(2^n)
  for(mask in 0:(2^length(rates)-1L)){
    active<-as.logical(intToBits(mask)[seq_along(rates)])
    covered<-unique(unlist(candidates[active]))
    index<-1L+sum(2^(covered-1L))
    out[index]<-out[index]+prod(ifelse(active,-expm1(-rates),exp(-rates)))
  }
  out
}

pair_union_exact <- function(singles,edges,rates) {
  n<-length(singles);out<-numeric(2^n)
  for(mask in 0:(2^n-1L)){
    selected<-which(as.logical(intToBits(mask)[seq_len(n)]))
    eligible<-which(edges[,1] %in% selected & edges[,2] %in% selected)
    plan<-cran_pair_union_plan(rep(0L,length(eligible)),edges[eligible,1]-1L,
                              edges[eligible,2]-1L,1L,n)
    correction<-cran_pair_union_evaluate(plan,-log(-expm1(-singles)),matrix(0,n,0),
      log(rates[eligible]),matrix(0,length(eligible),0))$log_partition
    out[mask+1L]<-exp(-sum(singles)-sum(rates)+sum(log(expm1(singles[selected])))+correction)
  }
  out
}

test_that('Poisson pair union normalizes with overlaps, duplicate types and incompatible proposals', {
  # Triangle plus a fourth edit; two types propose exactly the same pair.
  edges<-rbind(c(1L,2L),c(2L,3L),c(1L,3L),c(1L,2L),c(3L,4L))
  singles<-c(.13,.27,.09,.4);rates<-c(.2,.4,.12,.07,.31)
  candidates<-c(lapply(1:4,identity),lapply(seq_len(nrow(edges)),function(i)edges[i,]))
  brute<-pair_union_brute(candidates,c(singles,rates),4L)
  exact<-pair_union_exact(singles,edges,rates)
  expect_equal(exact,brute,tolerance=1e-12)
  expect_equal(sum(exact),1,tolerance=1e-12)
  expect_equal(exact[1],exp(-sum(singles)-sum(rates)),tolerance=1e-12)
  expect_equal(sum(exact[-1]/(1-exact[1])),1,tolerance=1e-12)
})

test_that('positive union gradients agree with finite differences', {
  edges<-rbind(c(0L,1L),c(1L,2L),c(0L,2L),c(0L,1L),c(3L,4L))
  plan<-cran_pair_union_plan(rep(0L,5),edges[,1],edges[,2],1L,5L)
  x<-c(log(c(2,7,3,1.5,4)),log(c(.2,.4,.12,.07,.31)))
  f<-function(z)cran_pair_union_evaluate(plan,z[1:5],cbind(diag(5),matrix(0,5,5)),
    z[6:10],cbind(matrix(0,5,5),diag(5)))
  h<-1e-6
  numeric<-vapply(seq_along(x),function(j){d<-rep(0,length(x));d[j]<-h;
    (f(x+d)$log_partition-f(x-d)$log_partition)/(2*h)},0.)
  expect_equal(drop(f(x)$gradient),numeric,tolerance=1e-7)
  # Reordering candidate types cannot change an unordered-union law.
  permutation<-c(5,3,1,4,2)
  other<-cran_pair_union_plan(rep(0L,5),edges[permutation,1],edges[permutation,2],1L,5L)
  result<-cran_pair_union_evaluate(other,x[1:5],matrix(0,5,0),x[5+permutation],matrix(0,5,0))
  expect_equal(result$log_partition,f(x)$log_partition,tolerance=1e-12)
})

test_that('exact computation factors components, closes frontiers and refuses excess work', {
  edges<-cbind(0:98,1:99)
  plan<-cran_pair_union_plan(rep(0L,99),edges[,1],edges[,2],1L,100L)
  expect_lte(max(cran_pair_union_complexity(plan)$states),8)
  # At f=1 every pair can be present or absent freely, giving exp(sum rates).
  result<-cran_pair_union_evaluate(plan,rep(0,100),matrix(0,100,0),rep(log(.2),99),matrix(0,99,0),8)
  expect_equal(result$log_partition,19.8,tolerance=1e-11)
  clique<-t(combn(0:6,2))
  large<-cran_pair_union_plan(rep(0L,nrow(clique)),clique[,1],clique[,2],1L,7L)
  expect_error(cran_pair_union_check(large,4),'No likelihood approximation or candidate truncation')
  expect_error(cran_pair_union_evaluate(large,rep(0,7),matrix(0,7,0),rep(0,nrow(clique)),matrix(0,nrow(clique),0),4),'exceeds pair_max_states')
})

test_that('addition-pair correction accounts for an overlapping migration', {
  # Singles a,b,c,d and migration {a,d}; paired additions {a,b},{b,c}.
  # All edits are observed. Compare the integrated migration formula against
  # independent enumeration of all seven proposal-presence indicators.
  edges<-rbind(c(1L,2L),c(2L,3L));single<-c(.2,.3,.4,.25);rates<-c(.35,.15);migration<-.27
  candidates<-c(lapply(1:4,identity),list(c(1L,4L)),list(c(1L,2L)),list(c(2L,3L)))
  expected<-tail(pair_union_brute(candidates,c(single,migration,rates),4),1)
  p<-log(c(rates,single,migration))
  evaluate<-function(p){
    r<-exp(p);singles<-r[3:6];mig<-r[7];ra<--expm1(-singles)
    u<-log(expm1(mig))-log(ra[1])-log(ra[4])
    du<-matrix(0,1,7);du[1,7]<-mig/-expm1(-mig)
    du[1,3]<--singles[1]/expm1(singles[1]);du[1,6]<--singles[4]/expm1(singles[4])
    x<-list(group=0L,a=1L,b=2L,log_weight=0,copied=0,mass=matrix(c(1,0),1,2),
            diamond=list(group=0L,a=2L,b=3L,log_weight=0,mass=1))
    m<-list(triangle=x,ng=1,Tclass=1L,meanT=1L,meanDiamond=2L,meanConsumer=NULL,
      ai=1L,ag=rep(0L,3),law='CS-pair-migration-add',pair_proposals='poisson',pair_max_states=1024,
      pair_union_plan=cran_pair_union_plan(c(0L,0L),edges[,1]-1L,edges[,2]-1L,1L,3L))
    rg<-matrix(0,4,7);rg[cbind(1:4,3:6)]<-1
    correction<-hawkesNet:::.cran_triangle_correction(m,p,log(singles[1:3]),rg,u,du)
    baseline<--sum(singles)-mig+sum(log(expm1(singles)))+hawkesNet:::.cran_softplus(u)
    list(loglik=baseline+correction$loglik,correction=correction)
  }
  result<-evaluate(p)
  expect_equal(exp(result$loglik),expected,tolerance=1e-12)
  h<-1e-6
  numeric<-vapply(seq_along(p),function(j){d<-rep(0,length(p));d[j]<-h;
    (evaluate(p+d)$correction$loglik-evaluate(p-d)$correction$loglik)/(2*h)},0.)
  expect_equal(drop(result$correction$gradient),numeric,tolerance=1e-7)
})

test_that('fit metadata preserves archived one-pair simulation', {
  fit<-list(law='CS-pair-migration-add',term_specification='pair_channels')
  expect_equal(hawkesNet:::.cran_sampler_law(fit),'CS-pair-migration-add-consumers')
  fit$pair_proposals<-'poisson'
  expect_equal(hawkesNet:::.cran_sampler_law(fit),'CS-pair-migration-add-consumers-poisson')
  fit$migration_degree<-TRUE
  expect_equal(hawkesNet:::.cran_sampler_law(fit),'CS-pair-migration-add-consumers-olddegree-poisson')
})

test_that('Poisson addition pairs reproduce the exact distribution of a complete bundle', {
  # All three pairs among a,b,c are transitive and share consumer z. u also
  # has an old b dependency, allowing a migration to overlap an addition pair.
  state<-list(names=c('u@1','a@1','b@2','c@1','b@1','z@1'),
    birth=rep(0,6),clock=rep(0,6),
    edges=rbind(c(0L,4L),c(1L,2L),c(1L,3L),c(2L,3L),c(5L,1L),c(5L,2L),c(5L,3L)))
  p<-numeric(37);p[1:2]<-log(c(.7,.4));p[20]<-log(.3)
  p[25:27]<-log(1.1);p[29:31]<--1000;p[35:37]<-log(.8)
  # Target weights: a=2,b=3,c=4,z=1 (one plus in-degree).
  # Common consumers: a/b=1, a/c=1, b/c=2 (z and a).
  edges<-rbind(c(1L,2L),c(1L,3L),c(2L,3L),c(1L,4L),c(2L,4L),c(3L,4L))
  weight<-c(6,8,12,2,3,4);witness<-c(1,1,2,0,0,0)
  rho<-1.1*weight/sum(weight)+.8*weight*witness/sum(weight*witness)
  # Singles add a,b,c,z; deletion removes b@1.
  candidates<-c(lapply(1:5,identity),list(c(2L,5L)),lapply(seq_len(nrow(edges)),function(i)edges[i,]))
  rates<-c(rep(.7/4,4),.4,.3,rho)
  expected<-pair_union_brute(candidates,rates,5L);expected[1]<-0;expected<-expected/sum(expected)
  release<-c(-1000,-1000,log(.2),rep(0,12));g<-cran_graph_create(state)
  set.seed(20261109);seen<-integer(32);selected<-counts<-0L
  for(i in seq_len(18000L)){
    x<-cran_graph_step(g,4,'CS-pair-migration-add-consumers-poisson',p,release,FALSE)
    if(length(x$sources)!=1L||x$sources!=0L)next
    selected<-selected+1L
    edits<-c(match(x$add[,2],c(1L,2L,3L,5L)),if(nrow(x$delete))5L)
    seen[1L+sum(2^(edits-1L))]<-seen[1L+sum(2^(edits-1L))]+1L
    counts<-counts+x$pair_proposals
  }
  expect_gt(selected,2000L)
  expect_equal(seen[1],0L)
  tolerance<-6*sqrt(expected*(1-expected)/selected)+.002
  expect_true(all(abs(seen/selected-expected)<tolerance))
  # This count includes all pair proposals, before repeated edits collapse.
  mean_pairs<-(1.1+.8+.3)/(-expm1(-sum(rates)))
  expect_equal(counts/selected,mean_pairs,tolerance=.12)
})

test_that('transitive and open-diamond Poisson streams can both contribute several pairs', {
  state<-list(names=paste0(c('u','a','b','c','x'),'@1'),birth=rep(0,5),clock=rep(0,5),
    edges=rbind(c(1L,4L),c(2L,4L),c(3L,4L)))
  p<-numeric(37);p[1:2]<-log(c(.5,.2));p[20]<--1000
  p[25:27]<-log(.6);p[29:31]<-log(1.3);p[35:37]<--1000
  edges<-rbind(c(1L,4L),c(2L,4L),c(3L,4L),c(1L,2L),c(1L,3L),c(2L,3L))
  candidates<-c(lapply(1:4,identity),lapply(seq_len(nrow(edges)),function(i)edges[i,]))
  rates<-c(rep(.5/4,4),rep(.6/3,3),rep(1.3/3,3))
  expected<-pair_union_brute(candidates,rates,4L);expected[1]<-0;expected<-expected/sum(expected)
  release<-c(-1000,-1000,log(.2),rep(0,12));g<-cran_graph_create(state)
  set.seed(20261110);seen<-integer(16);n<-multiple<-deletions<-0L
  for(i in seq_len(14000L)){
    x<-cran_graph_step(g,4,'CS-pair-migration-add-consumers-poisson',p,release,FALSE)
    if(length(x$sources)!=1L||x$sources!=0L)next
    n<-n+1L;multiple<-multiple+(x$pair_proposals>1)
    index<-1L+sum(2^(x$add[,2]-1L));seen[index]<-seen[index]+1L
    deletions<-deletions+nrow(x$delete)
  }
  expect_equal(deletions,0L);expect_gt(n,2000L);expect_gt(multiple,500L)
  expect_true(all(abs(seen/n-expected)<6*sqrt(expected*(1-expected)/n)+.002))
})
