test_that("focal diamond scores equal the single-edit and joint changes", {
  n <- 4L
  edges <- as.matrix(expand.grid(0:3,0:3));edges <- edges[edges[,1]!=edges[,2],,drop=FALSE]
  rooted <- function(A,u=1L) {
    v <- which(A[u,]!=0)
    if(length(v)<2L)return(0)
    sum(apply(combn(v,2),2,function(pair) {
      if(A[pair[1],pair[2]]||A[pair[2],pair[1]])return(0)
      sum(A[pair[1],-u]*A[pair[2],-u])
    }))
  }
  single_error <- pair_error <- 0
  for(mask in 0:(2^nrow(edges)-1)) {
    present <- as.logical(intToBits(mask)[seq_len(nrow(edges))])
    A <- matrix(0,n,n);A[edges[present,,drop=FALSE]+1L]<-1
    state <- list(names=c("u@1","b@1","b@2","z@1"),birth=rep(0,n),clock=rep(0,n),edges=edges[present,,drop=FALSE])
    g <- cran_graph_create(state)
    score <- cran_diamond_completion(g,0L)
    before <- rooted(A)
    for(v in 2:n) {
      B<-A;B[1,v]<-1-A[1,v]
      expected<-(rooted(B)-before)*(1-2*A[1,v])
      single_error<-max(single_error,abs(score[v]-expected))
    }
    if(A[1,2]&&!A[1,3]) {
      B<-A;B[1,2]<-0;B[1,3]<-1
      cross <- if(A[2,3]||A[3,2])0 else sum(A[2,-1]*A[3,-1])
      pair_error<-max(pair_error,abs(rooted(B)-before-(score[3]-score[2]-cross)))
    }
  }
  expect_equal(single_error,0)
  expect_equal(pair_error,0)
})

test_that("single-channel simulation includes the diamond and migration preferences", {
  state <- list(names=c("u@1","b@1","b@2","w@1","z@1"),birth=rep(0,5),clock=rep(0,5),
                edges=rbind(c(0L,1L),c(1L,4L),c(2L,4L),c(3L,4L)))
  p<-numeric(19);p[c(1,2)]<-log(c(1.4,.3));p[6]<-5;p[8]<-1.2
  # Candidates are b@2, w@1 and z@1; the first two each complete one diamond.
  weights<-exp(c(.5+1.2,.5,0));rho<-1.4*weights/sum(weights)
  expected<-(-expm1(-rho))/(-expm1(-1.7))
  rp<-c(-1000,-1000,log(.2),rep(0,12))
  g<-cran_graph_create(state);seen<-integer(3);total<-0L;set.seed(60261107)
  for(i in seq_len(14000L)) {
    a<-cran_graph_step(g,4,"CS-1-patterns",p,rp,FALSE)
    if(length(a$sources)==1L&&a$sources==0L) {
      total<-total+1L;seen<-seen+as.integer(2:4 %in% a$add[,2])
    }
  }
  expect_gt(total,2000)
  expect_equal(seen/total,expected,tolerance=.035)
})

test_that("migration pairs have independent preferences in the native sampler", {
  state<-list(names=c("u@1","b@1","b@2","b@3","c@1"),birth=rep(0,5),clock=rep(0,5),
               edges=rbind(c(0L,1L),c(4L,3L)))
  p<-numeric(27);p[1:2]<-log(.001);p[20]<-log(.5);p[25]<-2
  rp<-c(-1000,-1000,log(.2),rep(0,12));g<-cran_graph_create(state)
  expected<-(-expm1(-(.5*c(1,4)/5+.001/3)))/(-expm1(-.502))
  count<-integer(2);total<-0L;set.seed(77061107)
  for(i in seq_len(12000L)) {
    a<-cran_graph_step(g,4,"CS-pair-patterns",p,rp,FALSE)
    if(length(a$sources)==1L&&a$sources==0L) {
      total<-total+1L;count<-count+as.integer(2:3 %in% a$add[,2])
    }
  }
  expect_gt(total,1800)
  expect_equal(count/total,expected,tolerance=.04)
  # A disabled migration stream must not change the single-edit mechanism.
  p[1:2]<-log(c(.8,.4));p[20]<--1e6
  draw<-function(par){set.seed(44261107);replicate(100,cran_graph_step(g,4,"CS-pair-patterns",par,rp,FALSE),simplify=FALSE)}
  before<-draw(p);p[25:27]<-c(-3,2,.4)
  expect_identical(draw(p),before)
})

test_that("single additions count shared target pairs once per existing dependency", {
  state <- list(names=c("u@1","w@1","v@1","z@1","q@1"),birth=rep(0,5),clock=rep(0,5),
                edges=rbind(c(0L,1L),c(3L,1L),c(3L,2L),c(4L,1L),c(4L,2L)))
  p <- numeric(19)
  p[1:2] <- log(c(1.4,.3))
  p[4] <- 10*log(2)
  rp <- c(-1000,-1000,log(.2),rep(0,12))
  g <- cran_graph_create(state)
  expected <- -expm1(-1.4*c(2,1,1)/4)/-expm1(-1.7)
  seen <- integer(3);total <- 0L
  set.seed(7251)
  for(i in seq_len(14000L)){
    update <- cran_graph_step(g,4,"CS-1-consumers",p,rp,FALSE)
    if(length(update$sources)==1L && update$sources==0L){
      total <- total+1L
      seen <- seen+as.integer(2:4 %in% update$add[,2])
    }
  }
  expect_gt(total,2000)
  expect_equal(seen/total,expected,tolerance=.04)
})

test_that("shared-consumer pairs are sampled in proportion to their witness mass", {
  state <- list(names=c("u@1","v@1","w@1","y@1","z@1","q@1","r@1"),
                birth=rep(0,7),clock=rep(0,7),
                edges=rbind(c(4L,1L),c(4L,2L),c(5L,1L),c(5L,3L),c(6L,1L),c(6L,2L)))
  p <- numeric(37)
  p[1:2] <- log(1e-6)
  p[c(20,25:27,29:31)] <- -1000
  p[35:37] <- log(.5)
  rp <- c(-1000,-1000,log(.2),rep(0,12))
  g <- cran_graph_create(state)
  seen <- integer(2)
  set.seed(9251)
  for(i in seq_len(14000L)){
    update <- cran_graph_step(g,4,"CS-pair-migration-add-consumers",p,rp,FALSE)
    if(length(update$sources)==1L && update$sources==0L && nrow(update$add)==2L){
      targets <- sort(update$add[,2])
      expect_true(identical(targets,c(1L,2L)) || identical(targets,c(1L,3L)))
      seen[targets[2]-1L] <- seen[targets[2]-1L]+1L
    }
  }
  expect_gt(sum(seen),1500)
  expect_equal(seen/sum(seen),c(.75,.25),tolerance=.04)
})
