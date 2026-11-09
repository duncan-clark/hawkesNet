test_that('shared-person addition correction equals complete-candidate enumeration', {
  x <- cbind(c(.2,.5,.8,1.1),c(1,0,1,0))
  birth <- c(-1,0,.5,1)
  shared <- c(1,0,1,0)
  p <- c(.3,.6,.2,.7,-.4)
  m <- list(ng=1L,ia=1:2,ageA=3L,shared_people=list(group=c(0L,0L),
       base=x[c(1,3),],birth=birth[c(1,3)],design=cbind(1,x[c(1,3),2])))
  base <- exp(drop(x%*%p[1:2])+birth*p[3]);prob <- base/sum(base)
  normalizer <- list(logsum=log(sum(base)),score=matrix(colSums(x*prob),1),birth=sum(birth*prob))
  actual <- .cran_people_add_sum(m,normalizer,p,4:5)
  design <- cbind(x,shared,shared*x[,2])
  weights <- base*exp(shared*(p[4]+p[5]*x[,2]));prob <- weights/sum(weights)
  expect_equal(actual$logsum,log(sum(weights)),tolerance=1e-12)
  expect_equal(drop(actual$score),unname(colSums(design*prob)),tolerance=1e-12)
  expect_equal(actual$birth,sum(birth*prob),tolerance=1e-12)
})

test_that('shared-person TERGM correction equals independent Bernoulli enumeration', {
  addition <- rbind(c(1,.4,1),c(1,.7,0),c(1,.2,0))
  deletion <- rbind(c(1,-.3),c(1,-.6))
  p <- c(-2,.8,.5,-1,.3,.7,-.5,1.2)
  shared_add <- c(TRUE,FALSE,TRUE);shared_del <- c(FALSE,TRUE)
  xa <- addition[shared_add,,drop=FALSE];xd <- deletion[shared_del,,drop=FALSE]
  m <- list(ng=1L,shared_people=list(index=6:8,addition=xa,deletion=xd,
    addition_index=1:3,deletion_index=4:5,addition_group=c(0L,0L),deletion_group=0L))
  eta_a <- drop(addition%*%p[1:3]);eta_d <- drop(deletion%*%p[4:5])
  H <- sum(.cran_softplus(c(eta_a,eta_d)))
  dh <- matrix(c(colSums(addition*plogis(eta_a)),colSums(deletion*plogis(eta_d))),1)
  actual <- .cran_people_tergm_sum(m,H,dh,p)
  eta_a <- eta_a+shared_add*(p[6]+addition[,3]*p[8]);eta_d <- eta_d+shared_del*p[7]
  expected_gradient <- c(colSums(addition*plogis(eta_a)),colSums(deletion*plogis(eta_d)),
    sum(shared_add*plogis(eta_a)),sum(shared_del*plogis(eta_d)),sum(shared_add*addition[,3]*plogis(eta_a)))
  expect_equal(actual$H,sum(.cran_softplus(c(eta_a,eta_d))),tolerance=1e-12)
  expect_equal(drop(actual$dh),expected_gradient,tolerance=1e-12)
})

test_that('shared personnel multiply pair rates rather than only renormalizing weights', {
  state <- list(names=c('u@1','b@1','b@2','c@1'),birth=rep(0,4),clock=rep(0,4),
    edges=rbind(c(0L,1L)),people=list(u=1L,b=1L,c=integer()))
  p <- numeric(27);p[1:2]<-log(c(.8,.3));p[20]<-log(.4)
  rp <- c(-1000,-1000,log(.2),rep(0,12))
  # There is a single migration pair, u -> b@1/b@2. Its shared-person
  # rate multiplier must equal a change to the unshared reference rate.
  equivalent <- p;equivalent[20] <- p[20]+log(3)
  draw <- function(parameters,law) {
    g <- cran_graph_create(state);set.seed(9310)
    replicate(200,cran_graph_step(g,4,law,parameters,rp,FALSE),simplify=FALSE)
  }
  expect_identical(draw(c(p,0,0,log(3)),'CS-pair-consumers-people'),
                   draw(equivalent,'CS-pair-consumers'))
  for(law in c('CS-1','CS-pair','TERGM-logdegree')) {
    q <- if(law=='CS-1')c(p[1:15],rep(0,4)) else if(law=='CS-pair')p else c(-2,rep(0,9),-1,rep(0,12))
    expect_identical(draw(c(q,0,0,0),paste0(law,'-consumers-people')),
                     draw(q,paste0(law,'-consumers')))
  }
})

test_that('personnel inputs retain the latest pre-origin union and reject ambiguous rows', {
  people <- list(bin=c(1,2,3),package=c('a','b','a'),people=list(c(1,1,2),2,3))
  clean <- .cran_validate_people(people)
  expect_identical(clean$bin,1:3)
  expect_identical(clean$people[[1]],1:2)
  state <- list(names=c('a@1','a@2','b@1'))
  expect_identical(.cran_people_at(state,people,3)$people,list(a=1:2,b=2L))
  expect_error(.cran_validate_people(within(people,people[[1]]<-1.5)),'positive integer')
  expect_error(.cran_validate_people(within(people,bin<-c(1,3,2))),'sorted integer')
  expect_error(.cran_validate_people(list(bin=c(1,1),package=c('a','a'),people=list(1L,2L))),
               'one row per package/bin')
})
