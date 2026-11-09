test_that('release selection cache preserves both component probabilities and all scores', {
  data <- list(
    time=c(2,3,4),new_count=c(0,1,0),
    initial=rbind(c(0,0,.5,1.2,.8,.3),c(1,.2,1,.4,.6,.9)),
    changes=rbind(c(2,2,.3,.5,.6,.7),c(0,0,.8,1.5,.9,.4)),
    change_ptr=c(0,1,2,2),
    selected=rbind(c(0,0,1),c(1,1,2),c(2,0,1),c(0,1,1),c(1,0,1),c(2,1,1)),
    selected_ptr=c(0,2,4,6))
  first <- c(log(c(.2,.4,.8)),c(-.1,-.02,.2,.3,-.1,.1,-.2,-.01,.3,.1,.2,-.1))
  second <- first;second[1:3] <- second[1:3]+c(.8,-.3,.5)
  cache <- cran_release_mixture_cache(data,first,3L)
  for(rates in list(c(-25,-26,-24),first[1:3],c(1,2,.7))){
    first[1:3] <- rates
    actual <- cran_release_mixture_cached(cache,first,second)
    expected <- lapply(list(first,second),function(p)cran_release_evaluate(data,p,3L))
    expect_equal(actual,expected,tolerance=1e-12)
  }
  prefix <- cran_release_mixture_cache(data,first,3L,2L)
  actual <- cran_release_mixture_cached(prefix,first,second)
  expected <- lapply(list(first,second),function(p){
    result <- cran_release_evaluate(data,p,3L)
    result$loglik[3] <- 0;result$gradient[3,] <- 0;result
  })
  expect_equal(actual,expected,tolerance=1e-12)
  changed <- first;changed[5] <- changed[5]+.1
  expect_error(cran_release_mixture_cached(cache,changed,second),'Stale release selection cache')
})
