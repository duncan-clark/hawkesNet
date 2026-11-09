cran_audit_timing_data <- function() {
  list(times=c(.2,1.2,2.2),initial=0,tz=rep(0,3),tw=rep(1,3),tx=rep(0,3),
       tu=rep(0L,3),tp=0:3,cz=numeric(),cu=integer(),cp=rep(0L,4),
       calendar=c(2,3,4),calendar_phase=0:2,release_class=c(0L,1L,2L))
}

test_that('publication calendar controls match direct Poisson likelihood and scores', {
  data <- cran_audit_timing_data()
  p <- c(log_mu=log(2),log_break=log(.1),log_reopening=log(5))
  expected_rates <- 2*c(2,.3,5)
  result <- loglik_cran_timing(data,p,'Poisson',start=0,end=3,
                              calendar_effects=TRUE,details=TRUE)
  expect_equal(result$loglik,sum(log(expected_rates))-sum(expected_rates),tolerance=1e-12)
  expect_equal(result$intensity,expected_rates,tolerance=1e-12)
  expect_equal(unname(result$gradient),c(3-sum(expected_rates),1-.6,1-10),tolerance=1e-12)
  expect_true(is.na(result$end_intensity))
  expect_equal(loglik_cran_timing(data,rev(p),'Poisson',start=0,end=3,
                                  calendar_effects=TRUE)$loglik,result$loglik)
  fitted <- fit_cran_timing(data,model='Poisson',start=0,end=3,calendar_effects=TRUE)
  expect_equal(fitted$loglik,-3,tolerance=1e-6)
  expect_lt(fitted$projected_gradient,1e-5)
})

test_that('all manuscript timing models have matching calendar and amplitude scores', {
  data <- cran_audit_timing_data()
  for(model in c('one_exponential','M0','M1a','M1b','M1c','M2')) {
    typed <- model=='M2'
    layout <- .cran_timing_layout(model,typed,typed,TRUE)
    p <- setNames(rep(.1,length(layout$names)),layout$names)
    p[layout$tau] <- log(c(.3,2)[seq_len(layout$k)])
    evaluate <- function(p) loglik_cran_timing(data,p,model,start=0,end=3,
      release_types=typed,type_by_kernel=typed,calendar_effects=TRUE)
    result <- evaluate(p)
    numerical <- vapply(seq_along(p),function(j) {
      d <- rep(0,length(p));d[j]<-1e-5
      (evaluate(p+d)$loglik-evaluate(p-d)$loglik)/2e-5
    },0.)
    expect_equal(unname(result$gradient),numerical,tolerance=2e-6,info=model)
  }
  expect_length(.cran_timing_names('M1b',calendar_effects=TRUE),12)
  expect_length(.cran_timing_names('M2',TRUE,TRUE,TRUE),14)
})

test_that('invalid timing data fail before native array access', {
  data <- cran_audit_timing_data()
  run <- function(d) loglik_cran_timing(d,c(0,-1,-1,-2,1),'M0',start=0,end=3)
  expect_error(run(within(data,tp<-c(0L,2L,1L,3L))),'pointers')
  expect_error(run(within(data,tu[1]<-9L)),'indices')
  expect_error(run(within(data,times[2]<-times[1])),'strictly increasing')
  expect_error(run(within(data,calendar[1]<- -1)),'calendar')
  expect_error(loglik_cran_timing(data,c(0,-1,-1,-2,1),'M0',start=2,end=1),'boundaries')
})

test_that('graph input validation rejects malformed states and native sampler vectors', {
  state <- list(names=c('u@1','v@1'),birth=c(0,0),clock=c(0,0),edges=matrix(integer(),0,2))
  expect_error(cran_graph_create(within(state,birth<-0)),'vectors must match')
  expect_error(cran_graph_create(within(state,names[2]<-names[1])),'unique')
  expect_error(cran_graph_create(within(state,names[2]<-'v@-1')),'nonnegative-integer')
  expect_error(cran_graph_create(within(state,edges<-matrix(c(0L,2L),1,2))),'absent vertex')
  expect_error(cran_graph_create(within(state,edges<-matrix(c(0L,1L,0L,1L),2,2,byrow=TRUE))),'unique')
  g <- cran_graph_create(state)
  expect_error(cran_graph_step(g,1,'PA',numeric(),rep(0,15)),'parameter vector')
  expect_error(cran_graph_step(g,1,'PA',rep(0,8),c(rep(-1000,3),rep(0,12))),'total release')
})

test_that('motif-free metadata preserves migration preferences and native layouts', {
  m <- list2env(list(names=c('addition_transitive','pair_diamond','addition_shared_consumer',
                'migration','popularity_add'),lower=rep(-12,5),upper=rep(12,5)))
  .cran_fix_motif_terms(m)
  expect_equal(m$lower[1:3],rep(0,3))
  expect_equal(m$upper[1:3],rep(0,3))
  expect_equal(m$lower[4:5],rep(-12,2))
  expect_identical(.cran_sampler_law(list(law='PA',term_specification='motif_free')),'PA')
  expect_identical(.cran_sampler_law(list(law='CS-pair',term_specification='motif_free',
                  shared_people=TRUE,mixture='two_component')),'CS-pair-consumers-people-mixture2')
})

test_that('timing fits supply physical parameters that initialize joint forecasts exactly', {
  data <- cran_audit_timing_data()
  state <- list(names='a@1',birth=0,clock=0,edges=matrix(integer(),0,2))
  history <- lapply(seq_along(data$times),function(i)list(time=data$times[i],
    families=0L,summary=0,release_class=data$release_class[i]))
  mark <- list(law='PA',parameters=rep(0,8),release_parameters=rep(0,15))
  for(model in c('Poisson','one_exponential','M1b','M1c','M2')) {
    typed <- model=='M2'
    fitted <- fit_cran_timing(data,model=model,start=0,end=2.5,
      release_types=typed,type_by_kernel=typed,calendar_effects=TRUE,control=list(maxit=0))
    actual <- simulate_cran_joint(state,history,mark,fitted$timing,fitted$calendar,
      origin=2.5,horizon=0,model=fitted$simulation_model)$initial_intensity
    expected <- loglik_cran_timing(data,fitted$parameters,model,start=0,end=2.5,
      release_types=typed,type_by_kernel=typed,calendar_effects=TRUE)$end_intensity
    expect_equal(actual,expected,tolerance=1e-12,info=model)
  }
})

test_that('custom M2 clipping is preserved from timing fit to forecast', {
  data <- cran_audit_timing_data()
  data$tz <- c(0,log(5),log(5));data$cz <- log(5);data$cu <- 0L;data$cp <- c(0L,0L,1L,1L)
  parameters <- c(log_mu=0,log_alpha_short=log(.2),log_alpha_long=log(.2),
    log_tau_short=log(.1),log_tau_long=log(2),gamma_B=.1,gamma_Delta=.8,delta=0)
  fitted <- fit_cran_timing(data,parameters,'M2',start=0,end=2.5,clip=.2,control=list(maxit=0))
  state <- list(names=paste0(letters[1:5],'@1'),birth=rep(0,5),clock=rep(0,5),
                edges=cbind(1:4,0L))
  history <- lapply(seq_along(data$times),function(i)list(time=data$times[i],
    families=0L,summary=data$tz[i]))
  mark <- list(law='PA',parameters=rep(0,8),release_parameters=rep(0,15))
  forecast <- function(clip=NULL)simulate_cran_joint(state,history,mark,fitted$timing,
    fitted$calendar,origin=2.5,horizon=0,model=fitted$simulation_model,clip=clip)$initial_intensity
  expected <- loglik_cran_timing(data,fitted$parameters,'M2',start=0,end=2.5,clip=.2)$end_intensity
  expect_equal(forecast(),expected,tolerance=1e-12)
  expect_gt(abs(forecast(log(4))-expected),1e-5)
  expect_error(forecast(-1),'clip must be finite and nonnegative')
})
