test_that('a shared complete-update component is normalized and has exact scores', {
  # Enumerate all four outcomes of a small release/edge law. The shared
  # component creates dependence between the two outcomes after marginalizing.
  outcomes <- expand.grid(release=0:1,edge=0:1)
  bernoulli <- function(p,x,y){
    probability <- plogis(sum(p*x))
    list(loglik=dbinom(y,1,probability,log=TRUE),gradient=outer(y-probability,x))
  }
  local_mocked_bindings(
    .cran_single_loglik=function(m,p,per_update)bernoulli(p,c(1,.4),outcomes$edge),
     .cran_mixture_release=function(m,parts)lapply(parts$release,
      function(p)bernoulli(p,c(1,.2,-.3,rep(.1,12)),outcomes$release))
  )
  m <- list(names=c('log_m_add','log_m_delete'),lower=c(-12,-12),upper=c(4,4),
    rateA=1L,rateD=2L,law='PA',n=4L,release_data=list(),families=1L)
  layout <- .cran_mixture_layout(m)
  p <- setNames(rep(.1,length(layout$names)),layout$names)
  p[layout$shift] <- c(.8,-.2,.5,1.3,-.4,-.7)
  actual <- .cran_mixture_loglik(m,p,TRUE)
  expect_equal(sum(exp(actual$loglik)),1,tolerance=1e-13)
  parts <- .cran_mixture_components(m,p)
  component <- function(k)exp(bernoulli(parts$edge[[k]],c(1,.4),outcomes$edge)$loglik+
    bernoulli(parts$release[[k]],c(1,.2,-.3,rep(.1,12)),outcomes$release)$loglik)
  expected <- (1-plogis(parts$logit))*component(1)+plogis(parts$logit)*component(2)
  expect_equal(exp(actual$loglik),expected,tolerance=1e-13)
  h <- 1e-6
  numeric <- sapply(seq_along(p),function(j){
    delta <- numeric(length(p));delta[j] <- h
    (.cran_mixture_loglik(m,p+delta,TRUE)$loglik-.cran_mixture_loglik(m,p-delta,TRUE)$loglik)/(2*h)
  })
  expect_equal(unname(actual$gradient),numeric,tolerance=2e-8)
  expect_equal(actual$edge+actual$release,actual$loglik,tolerance=1e-13)
  collapsed <- p;collapsed[layout$shift] <- 0
  single <- bernoulli(parts$edge[[1]],c(1,.4),outcomes$edge)$loglik+
    bernoulli(parts$release[[1]],c(1,.2,-.3,rep(.1,12)),outcomes$release)$loglik
  expect_equal(.cran_mixture_loglik(m,collapsed,TRUE)$loglik,single,tolerance=1e-13)
})

test_that('component shifts cover all rate classes and preserve target preferences', {
  m <- list(names=c('log_m_add','log_m_delete','copy','log_m_migration',paste0('pair',1:9)),
    lower=rep(-12,13),upper=rep(4,13),rateA=1L,rateD=2L,meanP=4L,
    meanT=5:7,meanDiamond=8:10,meanConsumer=11:13,law='CS-pair-migration-add')
  layout <- .cran_mixture_layout(m)
  p <- setNames(seq_along(layout$names)/100,layout$names)
  parts <- .cran_mixture_components(m,p)
  expect_equal(parts$edge[[2]][3],parts$edge[[1]][3])
  expect_equal(parts$edge[[2]][5:13]-parts$edge[[1]][5:13],rep(p['mixture_pair_shift'],9),ignore_attr=TRUE)
  expect_equal(parts$release[[2]][4:15],parts$release[[1]][4:15])
  expect_equal(sum(layout$lower<layout$upper),length(m$names)+15+8)
  expect_equal(layout$lower[layout$shift['mixture_new_shift']],0)
})

test_that('named mixture vectors are matched by coefficient rather than position', {
  m <- list(names=c('log_m_add','log_m_delete'),lower=c(-12,-12),upper=c(4,4),
            rateA=1L,rateD=2L,law='PA')
  layout <- .cran_mixture_layout(m)
  p <- setNames(seq_along(layout$names)/10,layout$names)
  expect_equal(.cran_mixture_components(m,p),.cran_mixture_components(m,rev(p)))
  bad <- p;names(bad)[2] <- names(bad)[1]
  expect_error(.cran_mixture_components(m,bad),'Names of joint mixture parameters')
})
