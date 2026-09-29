construction_history <- function(varied_times = FALSE) {
  g <- network::network.initialize(4L, directed = FALSE)
  network::set.vertex.attribute(g, "time", if (varied_times) c(.1,.2,.3,.4) else rep(0,4))
  network::add.edges(g, 1L, 2L, names.eval = "time", vals.eval = .5)
  g
}

construction_params <- function() list(mu=3, K=.2, beta_overall=1.2,
  m=1.3, beta_edges=.4, node_lambda=0, CS_params=c(0,.8,.15,-.03))

construction_sets <- function(g) {
  pairs <- t(utils::combn(seq_len(network::network.size(g)), 2L))
  present <- vapply(seq_len(nrow(pairs)), function(j)
    length(network::get.edgeIDs(g,pairs[j,1],pairs[j,2])) > 0L, logical(1))
  pairs <- pairs[!present,,drop=FALSE]
  states <- t(vapply(0:(2^nrow(pairs)-1), function(mask)
    as.logical(intToBits(mask)[seq_len(nrow(pairs))]), logical(nrow(pairs))))
  marks <- lapply(seq_len(nrow(states)), function(j) {
    h <- network::network.copy(g)
    on <- which(states[j,])
    if (length(on)) network::add.edges(h,pairs[on,1],pairs[on,2],
                                      names.eval="time",vals.eval=rep(1,length(on)))
    h
  })
  list(pairs=pairs,states=states,marks=marks,size=rowSums(states))
}

construction_distribution <- function(mode, params, g, support,
                                      nonempty=FALSE, rhs="edges + triangles + star(c(2,3))") {
  vapply(support$marks, function(h) PMF_mark_CS(1,params,g,mark=h,
    formula_RHS=rhs,truncation=4,cs_mode=mode,
    condition_nonempty=nonempty)$mark_density, numeric(1))
}

test_that("both new CS constructions normalize over complete unordered supports", {
  g <- construction_history(); support <- construction_sets(g)
  params <- construction_params()
  for (mode in c("independent","size_conditional")) {
    p <- construction_distribution(mode,params,g,support)
    expect_equal(sum(p),1,tolerance=1e-11,info=mode)
    expect_true(all(p > 0),info=mode)
    p_nonempty <- construction_distribution(mode,params,g,support,TRUE)
    expect_equal(p_nonempty[1],0,info=mode)
    expect_equal(sum(p_nonempty),1,tolerance=1e-11,info=mode)
    expect_equal(p_nonempty[-1],p[-1]/(1-p[1]),tolerance=1e-11,info=mode)
  }
})

test_that("CS-1 freezes single-edge statistics and uses the collapsed product law", {
  g <- construction_history(); support <- construction_sets(g)
  params <- construction_params()
  rhs <- "edges + triangles + star(c(2,3))"
  before <- .cs_change_model(g,rhs)$statistics()
  changes <- t(vapply(seq_len(nrow(support$pairs)), function(j) {
    h <- network::network.copy(g)
    network::add.edges(h,support$pairs[j,1],support$pairs[j,2])
    unname(.cs_change_model(h,rhs)$statistics()-before)
  },numeric(4)))
  # All endpoint ages are equal here, so the age factor cancels exactly.
  weights <- stats::plogis(as.vector(changes %*% params$CS_params))
  r <- -expm1(-params$m*weights/sum(weights))
  expected <- apply(support$states,1,function(on)prod(ifelse(on,r,1-r)))
  observed <- construction_distribution("independent",params,g,support)
  expect_equal(observed,expected,tolerance=1e-11)
  expect_equal(observed[1],exp(-params$m),tolerance=1e-11)
  # No single absent edge closes a triangle, although a two-edge update can.
  params$CS_params <- c(0,0,0,0)
  at_zero <- construction_distribution("independent",params,g,support)
  params$CS_params[2] <- 2
  expect_equal(construction_distribution("independent",params,g,support),
               at_zero,tolerance=1e-11)
})

test_that("CS-2 separates the truncated-Poisson count from structural selection", {
  g <- construction_history(); support <- construction_sets(g)
  params <- construction_params(); d <- nrow(support$pairs)
  expected_count <- stats::dpois(0:d,params$m)/stats::ppois(d,params$m)
  p <- construction_distribution("size_conditional",params,g,support)
  actual_count <- vapply(0:d,function(k)sum(p[support$size==k]),numeric(1))
  expect_equal(actual_count,expected_count,tolerance=1e-11)
  changed <- params; changed$CS_params <- c(3,-1,.8,-.1)
  changed$beta_edges <- 2
  q <- construction_distribution("size_conditional",changed,g,support)
  expect_equal(vapply(0:d,function(k)sum(q[support$size==k]),numeric(1)),
               expected_count,tolerance=1e-11)
  edge_only <- params; edge_only$CS_params[1] <- 9
  expect_equal(construction_distribution("size_conditional",edge_only,g,support),
               p,tolerance=1e-11)
})

test_that("CS-2 weights interactions created by the whole simultaneous update", {
  g <- construction_history(); support <- construction_sets(g)
  params <- construction_params(); params$CS_params <- c(0,log(3),0,0)
  p <- construction_distribution("size_conditional",params,g,support)
  selected <- function(edges) {
    keys <- paste(support$pairs[,1],support$pairs[,2],sep="-")
    desired <- keys %in% edges
    which(apply(support$states,1,function(on)identical(unname(on),desired)))
  }
  triangle <- selected(c("1-3","2-3"))
  no_triangle <- selected(c("1-3","1-4"))
  expect_length(triangle,1L); expect_length(no_triangle,1L)
  expect_equal(p[triangle]/p[no_triangle],3,tolerance=1e-11)
})

test_that("CS-1 and CS-2 simulation, direct evaluation and parameter caches agree", {
  g <- construction_history(TRUE); params <- construction_params()
  params$node_lambda <- .7
  rhs <- "edges + triangles + star(c(2,3))"
  set.seed(20260922)
  for (mode in c("independent","size_conditional")) {
    for (j in seq_len(8)) {
      draw <- PMF_mark_CS(1,params,g,generate_mark=TRUE,formula_RHS=rhs,
                          truncation=4,cs_mode=mode)
      direct <- PMF_mark_CS(1,params,g,mark=draw$mark_sample,formula_RHS=rhs,
                            truncation=4,cs_mode=mode)
      expect_equal(draw$log_mark_sample_density,direct$log_mark_density,
                   tolerance=1e-10,info=mode)
      expect_equal(direct$log_density_func(params),direct$log_mark_density,
                   tolerance=1e-10,info=mode)
      changed <- params; changed$m <- .8; changed$node_lambda <- 1.1
      changed$beta_edges <- .9; changed$CS_params <- c(.2,-.4,.3,-.07)
      fresh <- PMF_mark_CS(1,changed,g,mark=draw$mark_sample,formula_RHS=rhs,
                           truncation=4,cs_mode=mode)
      expect_equal(direct$log_density_func(changed),fresh$log_mark_density,
                   tolerance=1e-10,info=mode)
      expect_true(network::network.size(draw$mark_sample)>4 ||
                    network::network.edgecount(draw$mark_sample)>1,info=mode)
      before <- .cs_change_model(g,rhs)$statistics()
      after <- .cs_change_model(draw$mark_sample,rhs)$statistics()
      expect_equal(unname(draw$mark_change_stats),unname(after-before),
                   tolerance=1e-10,info=mode)
    }
  }
  expect_equal(network::network.edgecount(g),1)
})

test_that("new CS modes normalize the joint birth-and-edge mark, not each birth separately", {
  old <- network::network.initialize(2,directed=FALSE)
  network::set.vertex.attribute(old,"time",c(0,0))
  params <- construction_params(); params$node_lambda <- .4
  rhs <- "edges + triangles + star(c(2,3))"
  for (mode in c("independent","size_conditional")) {
    mass <- conditional_mass <- numeric(5)
    for (b in 0:4) {
      born <- network::network.copy(old)
      if (b) network::add.vertices(born,b)
      network::set.vertex.attribute(born,"time",c(0,0,rep(1,b)))
      cands <- get_truncated_candidates(born,2+b,2,3,"node_entrance",FALSE)
      d <- length(cands$tails)
      for (mask in 0:(2^d-1)) {
        h <- network::network.copy(born)
        selected <- which(as.logical(intToBits(mask)[seq_len(d)]))
        if (length(selected)) network::add.edges(h,cands$tails[selected],
          cands$heads[selected],names.eval="time",vals.eval=rep(1,length(selected)))
        arguments <- list(time=1,params=params,mark_filtration=old,mark=h,
          formula_RHS=rhs,truncation=3,cs_mode=mode)
        q <- do.call(PMF_mark_CS,c(arguments,list(condition_nonempty=FALSE)))
        qc <- do.call(PMF_mark_CS,c(arguments,list(condition_nonempty=TRUE)))
        mass[b+1] <- mass[b+1]+q$mark_density
        conditional_mass[b+1] <- conditional_mass[b+1]+qc$mark_density
        if (b==0 && mask==0) excluded_empty <- q$mark_density
      }
    }
    expect_equal(mass,stats::dpois(0:4,params$node_lambda),tolerance=1e-11,info=mode)
    expected <- mass; expected[1] <- expected[1]-excluded_empty
    expected <- expected/(1-excluded_empty)
    expect_equal(conditional_mass,expected,tolerance=1e-11,info=mode)
    # Account analytically for omitted B>4 rather than silently truncating the
    # birth law when testing a finite portion of its infinite support.
    expect_equal(sum(conditional_mass)+
      stats::ppois(4,params$node_lambda,lower.tail=FALSE)/(1-excluded_empty),
      1,tolerance=1e-11,info=mode)
  }
})

test_that("explicit nonlegacy CS modes require their count or attempt parameter", {
  g <- construction_history(); params <- construction_params(); params$m <- NULL
  rhs <- "edges + triangles + star(c(2,3))"
  for (mode in c("independent","size_conditional","joint")) {
    required <- expected_params_PMF_mark_CS(g,rhs,cs_mode=mode)$required
    expect_true("m" %in% required,info=mode)
    expect_error(validate_params_for_PMF(params,PMF_mark_CS,g,
      formula_RHS=rhs,cs_mode=mode),"requires the following parameters: m")
    expect_error(PMF_mark_CS(1,params,g,generate_mark=TRUE,formula_RHS=rhs,
      cs_mode=mode),"requires params")
  }
  expect_false("m" %in% expected_params_PMF_mark_CS(g,rhs,cs_mode="legacy")$required)
  # R's `$m` partially matches `mu` when m is absent. Compatibility dispatch
  # must test the exact parameter name rather than interpreting mu as m.
  params$node_lambda <- .7
  mark <- construction_sets(g)$marks[[2]]
  automatic <- PMF_mark_CS(1,params,g,mark=mark,formula_RHS=rhs,truncation=4)
  legacy <- PMF_mark_CS(1,params,g,mark=mark,formula_RHS=rhs,truncation=4,cs_mode="legacy")
  expect_equal(automatic$log_mark_density,legacy$log_mark_density)
  changed <- params; changed$mu <- params$mu * 3
  changed_density <- PMF_mark_CS(1,changed,g,mark=mark,formula_RHS=rhs,
                                 truncation=4,cs_mode="legacy")
  expect_equal(changed_density$log_mark_density,legacy$log_mark_density)
  expect_equal(legacy$log_density_func(changed),legacy$log_mark_density)
})

test_that("the public CS-2 fitter refuses an unidentified free edge-count coefficient", {
  g <- construction_history(); params <- construction_params()
  params$CS_params <- c(0,.8)
  fit_args <- list(params_init=params,time_window=c(0,1),mark_filtration=g,
    PMF_mark=PMF_mark_CS,cs_mode="size_conditional",formula_RHS="edges + triangles",
    truncation=4,maxit=5,method="L-BFGS-B",get_hessian=FALSE,verbose=FALSE)
  expect_error(suppressMessages(do.call(fit_hawkesNet,fit_args)),
               "CS-2 conditions on edge count: fix CS_params1")
  fit_args$formula_RHS <- "triangles + edges"
  expect_error(suppressMessages(do.call(fit_hawkesNet,fit_args)),
               "CS-2 conditions on edge count: fix CS_params2")
  fit_args$formula_RHS <- "edges + triangles"
  fit_args$params_init$CS_params <- c(edge=0,triangle=.8)
  expect_error(suppressMessages(do.call(fit_hawkesNet,fit_args)),
               "CS-2 conditions on edge count: fix CS_params.edge")
  fit_args$formula_RHS <- "edges"; fit_args$params_init$CS_params <- 0
  expect_error(suppressMessages(do.call(fit_hawkesNet,fit_args)),
               "CS-2 conditions on edge count: fix CS_params")
})

test_that("fixing CS-2's cancelling coefficient permits a production cached fit", {
  params <- construction_params(); params$CS_params <- c(0,.6); params$node_lambda <- .8
  options <- list(cs_mode="size_conditional",formula_RHS="edges + triangles",truncation=4)
  set.seed(20260922)
  sim <- do.call(sim_hawkesNet,c(list(params=params,time_window=c(0,3),
    PMF_mark=PMF_mark_CS,cond_intensity=cond_intensity,verbose=FALSE),options))
  fit <- suppressMessages(do.call(fit_hawkesNet,c(list(params_init=params,
    time_window=c(0,3),mark_filtration=sim$net,PMF_mark=PMF_mark_CS,
    fixed_params=c("mu","K","beta_overall","beta_edges","node_lambda",
                   "CS_params1","CS_params2"),method="L-BFGS-B",maxit=80,
    get_hessian=FALSE,verbose=FALSE),options)))
  expect_equal(fit$fit$convergence,0)
  expect_identical(fit$params$CS_params,params$CS_params)
})

test_that("new CS marks preserve exact categorical support and tiny probabilities", {
  old <- construction_history()
  network::set.vertex.attribute(old,"group",rep("a",4))
  born_a <- network::network.copy(old)
  network::add.vertices(born_a,1)
  network::set.vertex.attribute(born_a,"time",c(rep(0,4),1))
  network::set.vertex.attribute(born_a,"group",rep("a",5))
  born_b <- network::network.copy(born_a)
  network::set.vertex.attribute(born_b,"group",c(rep("a",4),"b"))
  params <- construction_params(); params$node_lambda <- .7
  params$vertex_categorical <- list(group=c(a=0))
  params$vertex_categorical_levels <- list(group=c("a","b"))
  rhs <- "edges + triangles + star(c(2,3))"
  for (mode in c("independent","size_conditional")) {
    evaluate <- function(p,mark) PMF_mark_CS(1,p,old,mark=mark,formula_RHS=rhs,
                                           truncation=3,cs_mode=mode)
    no_attributes <- params
    no_attributes$vertex_categorical <- NULL
    no_attributes$vertex_categorical_levels <- NULL
    base <- evaluate(no_attributes,born_a)
    for (pa in c(0,1)) {
      p <- params; p$vertex_categorical$group <- c(a=pa)
      a <- evaluate(p,born_a); b <- evaluate(p,born_b)
      forbidden <- if(pa==0) a else b
      permitted <- if(pa==0) b else a
      expect_identical(forbidden$log_mark_density,-Inf,info=mode)
      expect_equal(forbidden$mark_density,0,info=mode)
      expect_equal(permitted$log_mark_density,base$log_mark_density,
                   tolerance=1e-11,info=mode)
      expect_equal(a$mark_density+b$mark_density,base$mark_density,
                   tolerance=1e-11,info=mode)
      other <- p; other$vertex_categorical$group <- c(a=1-pa)
      expect_equal(a$log_density_func(other),evaluate(other,born_a)$log_mark_density,
                   tolerance=1e-11,info=mode)
      expect_equal(b$log_density_func(other),evaluate(other,born_b)$log_mark_density,
                   tolerance=1e-11,info=mode)
    }
    tiny <- params; tiny$vertex_categorical$group <- c(a=1e-200)
    a <- evaluate(tiny,born_a); b <- evaluate(tiny,born_b)
    expect_equal(a$log_mark_density-base$log_mark_density,log(1e-200),
                 tolerance=1e-11,info=mode)
    expect_equal(a$mark_density+b$mark_density,base$mark_density,
                 tolerance=1e-11,info=mode)
    expect_equal(a$log_density_func(tiny),a$log_mark_density,
                 tolerance=1e-11,info=mode)
  }
  expect_equal(expand_vertex_categorical_probs(c(b=.3,a=.1),c("a","b","c")),
               c(a=.1,b=.3,c=.6),tolerance=1e-12)
})
