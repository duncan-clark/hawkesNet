test_that('a shared mixture component preserves release and edit conditioning', {
  state <- list(names=c('a@1','b@1'),birth=c(0,0),clock=c(0,0),
                edges=matrix(integer(),0,2))
  weight <- .3
  release_means <- c(.3,2.5)
  addition_means <- c(.2,1.2)
  fit <- list(law='PA',mixture='two_component',
    parameters=c(log(addition_means[1]),log(.1),rep(0,6),
                 log(addition_means[2]/addition_means[1]),0,0,0),
    release_parameters=c(log(release_means[1]),-1000,-1000,rep(0,12),
                         log(release_means[2]/release_means[1]),0,0,qlogis(weight)))
  expect_silent(.cran_check_fit(fit))
  expect_identical(.cran_sampler_law(fit),'PA-mixture2')
  set.seed(20261118)
  draws <- replicate(6000L,{
    update <- cran_graph_step(cran_graph_create(state),1,.cran_sampler_law(fit),
                              fit$parameters,fit$release_parameters,FALSE)
    c(component=update$mixture_component,births=length(update$born),
      zero_additions=nrow(update$add)==0)
  })
  expect_true(all(draws['births',]>=1))
  expect_equal(mean(draws['component',]==2),weight,tolerance=.025)
  for(component in 1:2){
    selected <- draws[,draws['component',]==component,drop=FALSE]
    lambda <- release_means[component]
    expected_mean <- lambda/-expm1(-lambda)
    expected_zero <- (exp(lambda*expm1(-addition_means[component]))-exp(-lambda))/-expm1(-lambda)
    expect_equal(mean(selected['births',]),expected_mean,tolerance=.08)
    expect_lt(abs(mean(selected['zero_additions',])-expected_zero),
              6*sqrt(expected_zero*(1-expected_zero)/ncol(selected))+.002)
  }
  expected <- sum(c(1-weight,weight)*
    (exp(release_means*expm1(-addition_means))-exp(-release_means))/-expm1(-release_means))
  expect_lt(abs(mean(draws['zero_additions',])-expected),
            6*sqrt(expected*(1-expected)/ncol(draws))+.002)
})

test_that('mixture rate shifts route through every standard native sampler', {
  state <- list(names=c('u@1','a@1','b@2','c@1','b@1','z@1'),birth=rep(0,6),clock=rep(0,6),
                edges=rbind(c(0L,4L),c(1L,2L),c(1L,3L),c(2L,3L),c(5L,1L),c(5L,2L),c(5L,3L)))
  release <- c(log(c(.2,.1,.7)),rep(0,12))
  release_shift <- c(.7,.4,-.2)
  shifts <- c(.6,-.4,.8,1.1)
  for(law in c('PA','CS-1-consumers','CS-pair-migration-add-consumers-capped2','TERGM-consumers')){
    parameters <- switch(law,
      PA=c(log(c(.7,.4)),rep(0,10)),
      'CS-1-consumers'=c(log(c(.7,.4)),rep(0,17)),
      'CS-pair-migration-add-consumers-capped2'=rep(0,37),
      'TERGM-consumers'=c(-3,rep(0,9),-.2,rep(0,8)))
    transformed <- parameters
    transformed[c(1,if(grepl('^TERGM',law))11 else 2)] <-
      transformed[c(1,if(grepl('^TERGM',law))11 else 2)]+shifts[1:2]
    if(grepl('^CS-pair',law)){
      transformed[20] <- transformed[20]+shifts[3]
      transformed[c(25:27,29:31,35:37)] <- transformed[c(25:27,29:31,35:37)]+shifts[4]
    }
    for(seed in 1:12){
      set.seed(seed)
      mixture <- cran_graph_step(cran_graph_create(state),100,paste0(law,'-mixture2'),
        c(parameters,shifts),c(release,release_shift,1000),FALSE)
      set.seed(seed);runif(1)
      direct <- cran_graph_step(cran_graph_create(state),100,law,transformed,
        c(release[1:3]+release_shift,release[-(1:3)]),FALSE)
      expect_equal(mixture$mixture_component,2L)
      mixture$mixture_component <- direct$mixture_component
      expect_equal(mixture,direct,info=law)
    }
  }
})

test_that('mixture validation supports legacy fits and catches incomplete vectors', {
  legacy <- list(law='PA',parameters=c(log(c(1,.2)),rep(0,6)),
                 release_parameters=c(log(c(.2,.2,.3)),rep(0,12)))
  expect_silent(.cran_check_fit(legacy))
  expect_identical(.cran_sampler_law(legacy),'PA')
  mixture <- legacy;mixture$mixture <- 'two_component'
  expect_error(.cran_check_fit(mixture),'Incompatible edge parameters')
  mixture$parameters <- c(mixture$parameters,rep(0,4))
  expect_error(.cran_check_fit(mixture),'19 finite release parameters')
  mixture$release_parameters <- c(mixture$release_parameters,rep(0,4))
  expect_silent(.cran_check_fit(mixture))
  mixture$release_parameters[19] <- NA_real_
  expect_error(.cran_check_fit(mixture),'19 finite release parameters')
})

test_that('complete mixture bundles agree with an independently enumerated law', {
  state <- list(names=paste0(c('a','b','c'),'@1'),birth=rep(0,3),clock=rep(0,3),
                edges=matrix(integer(),0,2))
  releases <- c(.4,1.8)
  additions <- c(.35,1.4)
  weight <- .3
  dyads <- rbind(c(0L,1L),c(0L,2L),c(1L,0L),c(1L,2L),c(2L,0L),c(2L,1L))
  expected <- numeric(64)
  for(mask in 1:63){
    selected <- as.logical(intToBits(mask)[1:6])
    count <- rowSums(matrix(selected,3,2,byrow=TRUE))
    releasing <- count>0
    component <- vapply(1:2,function(k){
      selected_probability <- -expm1(-releases[k]/3)
      release_probability <- prod(ifelse(releasing,selected_probability,1-selected_probability))/-expm1(-releases[k])
      inclusion <- -expm1(-additions[k]/2)
      edit_probability <- prod(inclusion^count[releasing]*(1-inclusion)^(2-count[releasing])/-expm1(-additions[k]))
      release_probability*edit_probability
    },0.)
    expected[mask+1] <- sum(c(1-weight,weight)*component)
  }
  expect_equal(sum(expected),1,tolerance=1e-13)
  parameters <- c(log(additions[1]),log(.1),rep(0,6),log(additions[2]/additions[1]),0,0,0)
  release_parameters <- c(-1000,-1000,log(releases[1]),rep(0,12),0,0,log(releases[2]/releases[1]),qlogis(weight))
  g <- cran_graph_create(state)
  counts <- integer(64)
  set.seed(20261120)
  for(i in 1:12000){
    update <- cran_graph_step(g,1,'PA-mixture2',parameters,release_parameters,FALSE)
    keys <- paste(update$add[,1],update$add[,2])
    selected <- paste(dyads[,1],dyads[,2]) %in% keys
    counts[1L+sum(2^(which(selected)-1L))] <- counts[1L+sum(2^(which(selected)-1L))]+1L
  }
  expect_equal(counts[1],0L)
  expect_true(all(abs(counts/sum(counts)-expected)<6*sqrt(expected*(1-expected)/sum(counts))+.001))
})
