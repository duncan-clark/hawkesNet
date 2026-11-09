test_that('kernel-specific amplitudes match direct likelihood and nest shared weights', {
  # One historical batch has all three release types. A later update is
  # scored at the resulting intensity. The entire window also has survival.
  data <- list(times=c(.2,.8), initial=rep(0,3), tz=rep(0,4),
    tw=c(rep(1/3,3),1), tx=c(rep(log(3),3),0), tu=c(0L,1L,2L,0L),
    tp=c(0L,3L,4L), cz=numeric(), cu=integer(), cp=rep(0L,3),
    calendar=rep(1,2), release_class=c(0L,1L,2L,0L))
  base <- c(log(2), log(.4), log(.2), log(.1), log(2), .5)
  shared <- log(c(3,2))
  split <- log(c(4,.5,.1,2.5))
  evaluate <- function(p, separate=TRUE)
    loglik_cran_timing(data,p,'batch',start=.2,end=1,details=TRUE,
                      release_types=TRUE,type_by_kernel=separate)
  actual <- evaluate(c(base,split))
  mass <- exp(base[2:3])
  delay <- exp(base[4:5])
  a <- 3^(-.5)*c(1+4+.5,1+.1+2.5)
  intensity <- 2 + sum(mass/delay*a*exp(-.6/delay))
  comp <- 2*.8 + sum(mass*(a*(-expm1(-.8/delay)) + (-expm1(-.2/delay))))
  expect_equal(actual$intensity[2],intensity,tolerance=1e-12)
  expect_equal(actual$loglik,log(2)+log(intensity)-comp,tolerance=1e-12)
  p <- c(base,split)
  for(i in 7:10){
    h <- numeric(10); h[i] <- 1e-5
    expect_equal(unname(actual$gradient[i]),
      (evaluate(p+h)$loglik-evaluate(p-h)$loglik)/2e-5,tolerance=1e-8)
  }
  expect_equal(evaluate(c(base,shared),FALSE)$loglik,
               evaluate(c(base,rep(shared,2)))$loglik,tolerance=1e-12)
  expect_error(loglik_cran_timing(within(data,rm(release_class)),p,'batch',
    start=.2,end=1,release_types=TRUE,type_by_kernel=TRUE),'Missing timing release classes')
})

test_that('joint forecasts initialize separate kernels from recorded line types', {
  state <- list(names=c('a@1','a@2','b@1'),birth=c(0,.5,.5),clock=c(0,.5,.5),
                edges=matrix(integer(),0,2))
  history <- list(list(time=.5,families=c(0L,0L,1L),summary=rep(0,3),
                       release_class=c(0L,2L,1L)))
  fit <- list(law='PA',parameters=rep(0,8),release_parameters=rep(0,15))
  timing <- c(mu=2,alpha_short=.4,alpha_long=.2,tau_short=.1,tau_long=2,
     delta=.5,gamma_B=0,gamma_Delta=0,
     setNames(log(c(4,.5,.1,2.5)),.cran_kernel_type_amplitudes))
  a <- 3^(-.5)*c(1+4+.5,1+.1+2.5)
  tau <- c(.1,2)
  z <- simulate_cran_joint(state,history,fit,timing,c(1,1),origin=1,
                          horizon=1e-9,model='M1',directed_waits=TRUE)
  expect_equal(z$initial_intensity,2+sum(c(.4,.2)/tau*exp(-.5/tau)*a),
               tolerance=1e-12)
  expect_named(z$directed_wait,c('Migration','Transitive','Diamond','Consumer'))
  history[[1]]$time <- .5
  z <- simulate_cran_joint(state,history,fit,timing,c(1,1),origin=.5,
                          horizon=1e-9,model='M1',
                          migration_cohort=list(target='a@2',exposed='b'))
  expect_equal(z$initial_intensity,2+sum(c(.4,.2)/tau*a),tolerance=1e-12)
  expect_length(z$migration_cohort$time,0)
})
