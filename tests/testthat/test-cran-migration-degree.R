test_that("migration tilt matches direct full-candidate normalization", {
  x <- matrix(c(.2,1,.6,0,.8,1,1.2,0),4,2,byrow=TRUE)
  birth <- c(-1,0,.5,1);degree <- c(log(2),0,log(8),0)
  beta <- c(.4,.7);age <- .3
  base <- exp(drop(x%*%beta)+age*birth);pr<-base/sum(base)
  m <- list(ng=1,degree_history=degree[c(1,3)],XH=x[c(1,3),],
    hist=list(group=c(0,0),features=cbind(0,birth[c(1,3)])))
  A <- list(logsum=log(sum(base)),score=matrix(colSums(x*pr),1),birth=sum(birth*pr))
  for(eta in c(-3,0,.8,3)){
    result <- hawkesNet:::.cran_degree_add_sum(m,A,beta,age,eta)
    weights <- base*exp(eta*degree);probability<-weights/sum(weights)
    expect_equal(result$logsum,log(sum(weights)),tolerance=1e-10)
    expect_equal(drop(result$score),colSums(x*probability),tolerance=1e-10)
    expect_equal(result$birth,sum(birth*probability),tolerance=1e-10)
    expect_equal(result$degree,sum(degree*probability),tolerance=1e-10)
  }
})

test_that("degree effects agree with equivalent constant-degree sampler", {
  state <- list(names=c("u@1","b@1","b@2","c@1"),birth=rep(0,4),clock=rep(0,4),
                edges=rbind(c(0L,1L)))
  rp <- c(-1000,-1000,log(.2),rep(0,12));eta<-.8
  for(law in c("CS-1","CS-pair-migration-add")){
    p<-numeric(if(law=="CS-1")19 else 37)
    p[1:2]<-log(c(.8,.3))
    if(length(p)>19){p[20]<-log(.6);p[c(25:27,29:31,35:37)]<--1000}
    equivalent<-p
    if(law=="CS-1")equivalent[8]<-eta*log(2) else equivalent[20]<-p[20]+eta*log(2)
    draw<-function(par,suffix){
      g<-cran_graph_create(state);set.seed(88231)
      replicate(300,cran_graph_step(g,4,paste0(law,suffix),par,rp,FALSE),simplify=FALSE)
    }
    expect_identical(draw(c(p,eta),"-consumers-olddegree"),draw(equivalent,"-consumers"))
  }
})

test_that("single migration degree sums the eligible old targets before edits", {
  state<-list(names=c('u@1','b@1','b@2','b@3','z@1'),birth=rep(0,5),clock=rep(0,5),
              edges=rbind(c(0L,1L),c(0L,2L),c(4L,1L)))
  update<-list(time=1,births=character(),clocks=numeric(),
               add=matrix(c('u@1','b@3'),1,2),delete=matrix(c('u@1','b@1'),1,2))
  cache<-cran_migration_degree_cache(state,list(update),'u@1',0L)
  # Old targets have two and one incoming edges, respectively.
  expect_equal(cache$history,log(1+2+1))
  expect_equal(cache$observed,log(1+2+1))
})
