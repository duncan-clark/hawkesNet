#' Conditional Intensity Function for hawkesGrowthNet Model
#'
#' Computes the conditional intensity function for the hawkesGrowthNet model,
#' defined as \deqn{\lambda(t,m | \mathcal{H}_t) = q(m | t,\mathcal{H}_t) \cdot \lambda(t | \mathcal{H}_t)}
#' where time decay is governed by a kernel and the mark PMF is supplied via \link{PMF_mark}.
#'
#' @param new_net Network. The network representing the mark to evaluate.
#' @param t Numeric. The time at which the mark occurs.
#' @param mark_filtration Network or NULL. Filtration/history of the network up to time \code{t}. Default: \code{NULL}.
#' @param PMF_mark Function. Function to compute the mark PMF (see \link{PMF_mark}).
#' @param params List. Model parameters, e.g. \code{mu}, \code{K}, \code{beta_overall}, \code{beta_edges}.
#' @param new_edge_hash Hash or NULL. Optional hash for fast edge lookup. Default: \code{NULL}.
#' @param ... Additional arguments passed to \code{PMF_mark}.
#'
#' @return List with elements:
#'   \item{result}{Numeric. The conditional intensity value.}
#'   \item{lambda}{Numeric. The background rate.}
#'   \item{kernel_sum}{Numeric. The sum of kernel decays.}
#'   \item{decays}{Numeric vector. Individual kernel decay values.}
#'   \item{diffs}{Numeric vector. Time differences.}
#'   \item{last_edge_probs}{Numeric vector. Edge probabilities from mark PMF.}
#'   \item{mark_density}{Numeric. Mark PMF value.}
#'   \item{mark_grad, decay_grad}{Gradients if computed.}
#'
#' @details
#' Computes the conditional intensity function for the \code{hawkesGrowthNet} model
#' defined as
#' \deqn{
#'   \lambda(t,m \mid \mathcal{H}_t) = 
#'     q(m \mid t,\mathcal{H}_{t}) \left( \lambda(t \mid \mathcal{H}_{t}) \right)
#' }
#' where the time decay kernel function, \eqn{g(t_1,t_2 \mid \mathcal{H}_t)},
#' is given by \eqn{g(t_1,t_2 \mid \mathcal{H}_t) = \exp(-\beta (t_1 - t_2))}.
#' The mark PMF, \eqn{q(m \mid t,\mathcal{H}_{t})}, is supplied via \link{PMF_mark}.
#' @seealso
#'  \code{\link{PMF_mark}}
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname loglik_hawkesGrowthNet
#' @export
cond_intensity <- function(new_net,
                           t,
                           mark_filtration = NULL,
                           PMF_mark,
                           params,
                           new_edge_hash = NULL,
                           ...){
    tmp <- PMF_mark(time = t,
                    params = params,
                    mark_filtration = mark_filtration,
                    mark = new_net,
                    generate_mark = FALSE,
                    new_edge_hash = new_edge_hash,
                    ...
                    )
    times <- get_times(mark_filtration)
    times <- times$times
    times <- times[times<t]
    diffs  <- t - times

    decays <- exp(-params$beta_overall*diffs)
    log_result <- tmp$log_mark_density + log(params$mu + params$K*sum(decays))
    result <- exp(log_result)

    return(list(result = result,
                lambda = params$mu,
                kernel_sum = params$K*sum(decays),
                decays = decays,
                diffs = diffs,
                last_edge_probs = tmp$edge_probs,
                mark_density = tmp$mark_density,
                mark_grad = tmp$mark_grad,
                decay_grad = tmp$decay_grad
                ))
}

#' Simulate Network Growth Under Hawkes Process
#'
#' Simulates network and event arrival under a Hawkes process, using a thinning approach, for the hawkesGrowthNet model.
#'
#' @param params List. Model parameters (e.g. \code{mu}, \code{K}, \code{beta_overall}, \code{beta_edges}, etc).
#' @param time_window Numeric vector of length 2. Time interval to simulate over.
#' @param PMF_mark Function. Computes/generates the mark PMF (see \link{PMF_mark}).
#' @param cond_intensity Function. Computes the conditional intensity (see \link{cond_intensity}).
#' @param hashed_edges Logical. If \code{TRUE}, use hashed edges for fast lookup. Default: \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, print simulation details. Default: \code{FALSE}.
#' @param mu_multiplier Numeric. Multiplier for background rate. Default: \code{10}.
#' @param joint_accept Logical. If \code{TRUE}, use joint acceptance criterion. Default: \code{FALSE}.
#' @param n_mark_sample Integer or NULL. If set, use importance sampling for marks. Default: \code{NULL}.
#' @param ... Additional arguments passed to \code{PMF_mark}.
#'
#' @return List with elements:
#'   \item{events}{List of event times and densities.}
#'   \item{net}{Final network at end of simulation.}
#'   \item{accept_probs}{Acceptance probabilities for each event.}
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[network]{as.edgelist}}
#'  \code{\link[hash]{hash}}
#' \code{\link{cond_intensity}}
#' \code{\link{PMF_mark}}
#' @rdname sim_hawkesGrowthNet
#' @export
sim_hawkesGrowthNet <- function(params,
                                time_window,
                                PMF_mark,
                                cond_intensity, 
                                hashed_edges = FALSE,
                                verbose = FALSE,
                                mu_multiplier = 10,
                                joint_accept = FALSE,
                                n_mark_sample = NULL,
                                ... ){
    t1 <- proc.time()
    ## simulate the background points (can only simulate their times right now)
    mu <- params$mu
    ## poisson in time lambda
    ## theta <- params$theta
    ## beta <- params$beta
    ## K <- params$K

    lambda <- mu_multiplier*mu

    ## Initialize the output list of events
    events = list()
    events$n = 0
    events$t = c()
    events$mark_density <- c()

    accept_probs <- c()

    ## propose points to be thinned:
    n_bg = rpois(1, lambda * (time_window[2] - time_window[1]))
    event_queue <- data.table(time = sort(runif(n_bg, min=0, max=time_window[2])))
    ## maintain order so no need to sort
    setkey(event_queue, time)
    ## Initialize the list to store new events
    new_events_list <- list()
    list_index <- 1
    tot <- 0
    tot_attempt <- 0
    current_net <- network::network(matrix(1),directed = FALSE)
    delete.vertices(current_net,1)
    while (nrow(event_queue) > 0) {
        t <- proc.time()
        current_event <- event_queue[1, ,drop = FALSE]
        event_queue <- event_queue[-1, ,drop = FALSE]  # Remove the processed event

        if(is.null(current_net %v% 'n')){
            mark_sample <- PMF_mark(time = current_event$time,
                                    params = params,
                                    mark_filtration = current_net,
                                    mark = NULL,
                                    generate_mark = TRUE,
                                    new_edge_hash = NULL,
                                    ...)
            net <- mark_sample$mark_sample
            accept <- 1
        }
        else{
                                        # get the mark samples
                                        # debug(PMF_mark)
            mark_sample <- PMF_mark(time = current_event$time,
                                    params = params,
                                    mark_filtration = current_net,
                                    mark = NULL,
                                    generate_mark = TRUE,
                                    generate_density = FALSE,
                                    new_edge_hash = NULL,
                                    ...)
            net <- mark_sample$mark_sample
            new_nodes <- (net %n% 'n') - (current_net %n% 'n')
            if(hashed_edges && length(net$mel)!=0){
                                        # hash the network edge list for fast lookup:
                edges <- network::as.edgelist(net)
                keys_vec <- paste(edges[,1], edges[,2], sep = "-")
                edge_hash <- hash::hash(keys = keys_vec, values = rep(TRUE, length(keys_vec)))
            }else{
                edge_hash <- NULL
            }
            if(joint_accept){
                intensity <- cond_intensity(new_net = net,
                                            t = current_event$time,
                                            mark_filtration = current_net,
                                            PMF_mark = PMF_mark,
                                            params = params,
                                            new_edge_hash = edge_hash,
                                            ...
                                            )$result
            }else{
                if(!is.null(n_mark_sample)){
                    imp_sample <- sapply(1:n_mark_sample,function(i){
                        mark_sample <- PMF_mark(time = current_event$time,
                                                params = params,
                                                mark_filtration = current_net,
                                                mark = NULL,
                                                generate_mark = TRUE,
                                                new_edge_hash = TRUE,
                                                ...
                                                )
                        net <- mark_sample$mark_sample
                        new_nodes <- (net %n% 'n') - (current_net %n% 'n')
                        if(hashed_edges){
                                        # hash the network edge list for fast lookup:
                            edges <- network::as.edgelist(net)
                            keys_vec <- paste(edges[,1], edges[,2], sep = "-")
                            edge_hash <- hash::hash(keys = keys_vec, values = rep(TRUE, length(keys_vec)))
                        }else{
                            edge_hash <- NULL
                        }
                        intensity <- cond_intensity(new_net = net,
                                                    t = current_event$time,
                                                    mark_filtration = current_net,
                                                    PMF_mark = PMF_mark,
                                                    params = params,
                                                    new_edge_hash = edge_hash,
                                                    ...
                                                    )
                        return(intensity$result/mark_sample$mark_sample_density)
                    })
                    intensity <- mean(imp_sample)
                }else{
                    tmp <- cond_intensity(new_net = net,
                                          t = current_event$time,
                                          mark_filtration = current_net,
                                          PMF_mark = PMF_mark,
                                          params = params,
                                          new_edge_hash = edge_hash,
                                          ...
                                          )
                    intensity <- tmp$lambda + tmp$kernel_sum
                }

                mark_sample <- PMF_mark(time = current_event$time,
                                        params = params,
                                        mark_filtration = current_net,
                                        mark = NULL,
                                        generate_mark = TRUE,
                                        new_edge_hash = NULL,
                                        ...
                                        )
                net <- mark_sample$mark_sample
            }

            accept <- intensity/lambda
        }
        accept_probs <- c(accept_probs,accept)
                                        # if we accept the point add it in
        if(verbose){print(paste0("accept prob is: ",accept))}
        if(runif(1) < accept){
            if(verbose) print('accepted!')
            current_net <- net
            events$t[length(events$t)+1] <- current_event$time
            events$n <- events$n + 1 
            if(length(events$t) >2){
                events$mark_density <- c(events$mark_density,mark_sample$mark_density)
            }
        }else{
            ## do nothing since we rejected the point
        }
        if(verbose){
            print(paste0("time is ",current_event$t, " size of net is ",
                         current_net %n% 'n',' number of edges is ',length(current_net$mel)))
            print(paste0("time is ",current_event$t, " this iteration of while loop took ",
                         round((proc.time()-t)[3],2)," seconds"))
        }
                                        # Concatenate new events to event_queue only if we have only one event left to go
        old_n <- dim(event_queue)[1]
        event_queue <- rbindlist(list(event_queue, rbindlist(new_events_list, use.names = TRUE)))
        new_events_list <- vector("list", length(new_events_list))  # Reset the list
        list_index <- 1
    }
    t1 <- proc.time() - t1
    print(paste0("simulation took ",round(t1[3],2)," seconds"))
    return(list(events = events,
                net = current_net,
                accept_probs = accept_probs))
}

#' Log-Likelihood for hawkesGrowthNet Model
#'
#' Computes the log-likelihood and optionally the gradients for a hawkesGrowthNet model, given observed network growth.
#'
#' @param params List. Model parameters (see \link{sim_hawkesGrowthNet}).
#' @param time_window Numeric vector. Time interval for likelihood computation.
#' @param mark_filtration Network. Filtration/history of the network over time.
#' @param PMF_mark Function. Computes the mark PMF.
#' @param edge_hash_list List or NULL. Optional list of edge hashes for fast lookup. Default: \code{NULL}.
#' @param verbose Logical. Print details. Default: \code{FALSE}.
#' @param do_grad Logical. If \code{TRUE}, also compute gradients. Default: \code{FALSE}.
#' @param ... Additional arguments to pass to \code{PMF_mark}.
#'
#' @return List with elements:
#'   \item{loglik}{Numeric. Log-likelihood value.}
#'   \item{grads}{Named vector/list of parameter gradients (if \code{do_grad = TRUE}).}
#' @examples
#' \dontrun{
#' if(interactive()){
#' data(net, package = "hawkesGrowthNet")
#' params_ba <-  list(mu = 10,
#' beta_overall = 0.1,
#' K = 0.1, beta_edges = 0.1)
#' mark_filtration <-  network::get.inducedSubgraph(net, v = 1:10)
#' loglik_ba <- loglik_hawkesGrowthNet(params = params_ba,
#'                                        time_window = c(0,max(get_times(mark_filtration)$times)),
#'                                        mark_filtration =  mark_filtration,
#'                                        PMF_mark = PMF_mark,  
#'                                        verbose = FALSE)
#' ## CS
#'require(ernm)
#' params_cs <-  list(mu = 10,
#'                      beta_overall = 0.1,
#'                       K = 0.1,
#'                       beta_edges = 0.1,
#'                       node_lambda = 1,
#'                       CS_params =  c(-10,0,0,0))
#' loglik_cs <- loglik_hawkesGrowthNet(params = params_cs,
#'                                        time_window = c(0,max(get_times(mark_filtration)$times)),
#'                                        mark_filtration =  mark_filtration,
#'                                        PMF_mark = PMF_mark, type = "CS",
#'                                        truncation = 1,
#'                                        formula_RHS = "edges + triangles + star(c(2,3))",
#'                                        max_node_time = 1,
#'                                        verbose = FALSE)
#'  }
#' }
#' @seealso
#'  \code{\link[network]{as.edgelist}}
#'  \code{\link[hash]{hash}}
#'  \code{\link{sim_hawkesGrowthNet}}
#'  \code{\link{cond_intensity}}
#' @rdname loglik_hawkesGrowthNet
#' @export
loglik_hawkesGrowthNet = function(params,
                                  time_window,
                                  mark_filtration,
                                  PMF_mark,
                                  edge_hash_list = NULL,
                                  verbose = FALSE,
                                  do_grad = FALSE,
                                  ...){
    t<-proc.time()
                                        # don't allow negative parameters in first 2
    if(any(sapply(params[1:min(length(params),2)],function(x){x<0}))){
        return(list(loglik = -(10**(100)),
                    grads = rep(0,length(params)))
               )
    }
    times <- get_times(mark_filtration)
    times <- times$times

    tval <- time_window[2]-time_window[1]
    max_t <- max(times)
    if(tval < max(times) - min(times)){
        stop("realization has points outside time window")
    }

                                        # do the sum of the intensities:
    intens_sum <- 0
    intens_list <- list()
    
                                        # if parallelize do that here with PSOCK for simplicity:
    intens_func <- function(i){
                                        # need to do this on the fly otherwise too storage intensive
        current_net <- filtration_to_net(mark_filtration,times[i],equals = TRUE)
        if("formula_RHS" %in% names(list(...))){
            model = createCppModel(as.formula(paste("current_net ~ ",list(...)$formula_RHS)))
        }else{
            model <- NULL
        }
        intensity <- cond_intensity(new_net = current_net,
                                    t = times[i],
                                    mark_filtration = current_net,
                                    PMF_mark = PMF_mark,
                                    params = params,
                                    model = model,
                                    ...)
        return(intensity$result)
                                        #return(intensity)
    }
    
    if("cores" %in% names(list(...))){
        if(verbose){
            print(paste0("using ",list(...)$cores," cores on ",length(times), " objects"))
        }
        cores <- list(...)$cores
                                        #cl <- makeForkCluster(cores)
        t <- proc.time()
                                        # intens_list <- parLapply(cl,X=1:length(times),fun = function(x){intens_func(x)})
        intens_list <- parallel::mclapply(seq_along(times), intens_func, mc.cores = cores,mc.preschedule=FALSE)
        print(paste0("intens list ", round((proc.time()-t)[3],2)," seconds"))
                                        #stopCluster(cl)
    }else{
        intens_list <- lapply(1:length(times),intens_func)
    }
                                        # tmp <- sapply(intens_list,function(x){x$result})
    tmp <- unlist(intens_list)
    if (exists("verbose") && isTRUE(verbose)) print(tmp)
    
    if(any(is.na(tmp))){
        tmp[is.na(tmp)] <- min(tmp[!is.na(tmp)])/2
    }
    if(sum(tmp==0)!=0){
        warning("some of the intens lists have zero - something is probably wrong")
        tmp[tmp==0] <- min(tmp[tmp>0])/2
        # print(summary(tmp))
    }
    intens_sum <- sum(log(tmp))

                                        # Integral due to kernel being density:


    max_t <- max(times)
    pieces <- sapply(times,function(x){
        (1-exp(-params$beta_overall*(tval-x)))
    })
    integral <- params$mu * tval + (1/params$beta_overall)*params$K*sum(pieces)
    loglik <- intens_sum - integral
                                        # print(paste0("integral is ",integral))
                                        # print(paste0("intens_sum is ",intens_sum))
                                        # print(paste0("trigger part of integral is  ",(1/params$beta_overall)*params$K*sum(pieces)))
                                        # print(paste0("result is :",loglik))
                                        # print(paste0("this iteration of loglik took ", round((proc.time()-t)[3],2)," seconds"))

                                        
    
    grads <- list()
    if(do_grad){
        # calcualte the gradients:

        kernel_sum <- sapply(intens_list,function(x){x$kernel_sum})
        decays <- lapply(intens_list,function(x){x$decays})
        diffs <- lapply(intens_list,function(x){x$diffs})
                                        # hawkes_mu
        grads$mu <- sum(1/(params$mu + params$K*kernel_sum)) - tval
                                        # hawkes K
        grads$K <-  sum(kernel_sum/(params$mu + params$K*kernel_sum)) -  (1/params$beta_overall)*sum(pieces)
                                        # hawkes beta
        grads$beta_overall <- sum(sapply(1:length(decays),function(i){
            (params$K * sum(-diffs[[i]]*decays[[i]])) / (params$mu + kernel_sum[i])
        })) +
            (-params$K)/(params$beta_overall**2) * ( length(kernel_sum) + sum((tval - times - 1/(params$beta_overall^2))*exp(-params$beta_overall*(tval-times))))
        
                                        # hawkes edge decay from mark generator:
        decay_grads <- sapply(intens_list,function(x){x$decay_grad})
        grads$beta_edges <- sum(decay_grads)
        
                                        # hawkes theta_params (from mark generator)
                                        # get the mark PMF grads:
        mark_grads <- lapply(intens_list,function(x){x$mark_grad})
        grads$CS_params <- colSums(do.call(rbind,mark_grads))
    }


                                        # TODO fix the mark grad issues with simplification

    t<-proc.time() - t
    if(verbose){
        print(paste0("this iteration of loglik took ", round(t[3],2)," seconds"))
    }

    return(list(loglik = loglik,
                grads = grads))
}

#' Fit hawkesGrowthNet Model via Maximum Likelihood
#'
#' Fits (estimates parameters for) a hawkesGrowthNet model by maximizing the log-likelihood.
#'
#' @param params_init List. Initial parameter values.
#' @param time_window Numeric vector. Time interval for model fitting.
#' @param mark_filtration Network. Network history/filtration.
#' @param PMF_mark Function. Computes the mark PMF.
#' @param trace Integer. Optimizer trace level. Default: 0.
#' @param REPORT Integer. Optimizer report frequency. Default: 10.
#' @param reltol Numeric. Optimization relative tolerance. Default: 1e-8.
#' @param maxit Integer. Maximum number of iterations.
#' @param get_hessian Logical. Compute Hessian for SEs? Default: \code{FALSE}.
#' @param fixed_params Character vector or NULL. Parameter names to hold fixed. Default: \code{NULL}.
#' @param ... Additional arguments passed to log-likelihood function.
#'
#' @return List with component:
#'   \item{fit}{Object returned by \code{optim} (MLE fit).}
#'
#' @seealso \code{\link{loglik_hawkesGrowthNet}}, \code{\link{sim_hawkesGrowthNet}}
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname fit_hawkesGrowthNet
#' @export
fit_hawkesGrowthNet <- function(params_init,
                                time_window,
                                mark_filtration,
                                PMF_mark,
                                trace = 0,
                                REPORT = 10,
                                reltol = 1e-8,
                                maxit,
                                get_hessian = FALSE,
                                fixed_params = NULL,
                                ...){
    params_init_old <- params_init
    if(!is.null(fixed_params)){
        for(k in fixed_params){
            params_init[[k]] <- NULL
        }
    }
    
    optim_func <- function(params,...){
        param_vec <- params
        params <- relist(params,skeleton = params_init)
        params[fixed_params] <- params_init_old[fixed_params]
        
                                        # Don't think K always needs to be less than 1 ? 
                                        # if(params$K>1){
                                        #   return(list(value = -10**(20),
                                        #               grad = NULL))
                                        # }

        result <- loglik_hawkesGrowthNet(params = params,
                                         time_window = time_window,
                                         mark_filtration = mark_filtration,
                                         PMF_mark = PMF_mark,
                                         ...)

                                        # =================
                                        # USE numDERIV TO DEBUG !
                                        # =================
        wrapper <- function(x,...){
            x <- relist(x, skeleton = params_init)
            return(loglik_hawkesGrowthNet(params = x,
                                          ...)$loglik)
        }

                                        # numgrad <- grad(func = wrapper,
                                        #                 x    = param_vec,
                                        #                 time_window = time_window,  # or whatever your data is
                                        #                 events      = events,
                                        #                 PMF_mark    = PMF_mark,
                                        #                 ...)
                                        #
                                        # print("Numeric gradient vs analytic gradient:")
                                        # print(numgrad)
                                        #print("analytic gradient")
                                        #print(unlist(result$grads))

                                        # print("params are:")
                                        # print(params)
                                        # print("Loglik is :")
                                        # print(result$loglik)
                                        # print("Grads are:")
                                        # print(result$grads)
        return(list(value = result$loglik,
                                        #grad = numgrad
                    grad = unlist(result$grads)
                    ))
    }

    fn_wrapper <- function(par,...) {
        res <- optim_func(par, ...)
        return(res$value)
    }
    gr_wrapper <- function(par, ...) {
        res <- optim_func(par, ...)

        return(res$grad)
    }
    t<-proc.time()
    fit <- optim(par = unlist(params_init),
                 fn = fn_wrapper,
                                        # gr = gr_wrapper,
                 gr = NULL,
                                        # method = 'BFGS',
                                        # method = 'CG',
                 method = "Nelder-Mead",
                 control = list(fnscale = -1,
                                trace=trace,
                                maxit=maxit,
                                reltol = reltol,
                                abstol = NULL),
                 hessian = get_hessian,
                 ...)
    print(paste0("fitting took ",round((proc.time()-t)[3],2)," seconds"))
    
                                        # Numerically estimate Hessian at optimal params
                                        # THIS IS GONNA TKE FOREVER - NEED TO CODE UP GRADIENT!
                                        # hessian_estimate <- numDeriv::hessian(
                                        #   func = function(p) {
                                        #     cat(sprintf("Parameters: %s\n",paste(round(p, 4), collapse = ", ")))
                                        #     p <- relist(p, skeleton = params_init)
                                        #     result <- -loglik_hawkesGrowthNet(params = p,
                                        #                                       time_window = time_window,
                                        #                                       mark_filtration = mark_filtration,
                                        #                                       PMF_mark = PMF_mark,
                                        #                                       ...)$loglik
                                        #     cat(sprintf("  --> loglik: %f\n", -result))
                                        #     
                                        #     
                                        #     return(result)
                                        #   },
                                        #   x = fit$par
                                        # )
                                        # 
                                        # # Fisher information approximation is the negative Hessian at the optimum
                                        # fisher_info <- hessian_estimate
                                        # vcov_matrix <- solve(fisher_info)
    

    
    return(list(fit=fit))
}


#' Compensator Process for hawkesGrowthNet
#'
#' Computes compensator (integrated intensity) values for a hawkesGrowthNet realization.
#'
#' @param params List. Model parameters.
#' @param time_window Numeric vector. Time interval.
#' @param mark_filtration Network. Network history/filtration.
#'
#' @return Numeric vector. Compensator values at observed event times.
#'
#' @seealso \code{\link{loglik_hawkesGrowthNet}}
#' @rdname ks_test_pval_hawkesGrowthNet
#' @export
compensators_hawkesGrowthNet <- function(params,
                                         time_window,
                                         mark_filtration){
    times <- get_times(mark_filtration)
    times <- times$times

    tval <- time_window[2]-time_window[1]
    max_t <- max(times)
    if(tval < max(times) - min(times)){
        stop("realization has points outside time window")
    }
                                        # Integral due to kernel being density:
    pieces <- sapply(times,function(x){
        (1-exp(-params$beta_overall*(tval-x)))
    })
    incremental <- sapply(1:length(times),function(i){
        integral <- params$mu * times[i] + (1/params$beta_overall)*params$K*sum(pieces[1:i])
        integral <- params$mu * times[i] + params$K*sum(pieces[1:i])
        return(integral)
    })
    return(incremental)
}

#' Kolmogorov-Smirnov Test for hawkesGrowthNet Compensators
#'
#' Performs a KS test on transformed compensator increments to assess model fit.
#'
#' @param params List. Model parameters.
#' @param time_window Numeric vector. Time interval.
#' @param mark_filtration Network. Network history/filtration.
#'
#' @return Numeric. p-value from KS test.
#'
#' @seealso \code{\link{compensators_hawkesGrowthNet}}
#' @rdname ks_test_pval_hawkesGrowthNet
#' @export
ks_test_pval_hawkesGrowthNet <- function(params,
                                         time_window,
                                         mark_filtration){
    compensators <- compensators_hawkesGrowthNet(params = unlist(params),
                                                 time_window = time_window,
                                                 mark_filtration = mark_filtration
                                                 )
    compensator_incs <- diff(compensators)
    test_dist <- 1 - exp(-compensator_incs)
    test <- ks.test(test_dist,"punif")
    return(test$p.value)
}


