#' Prepare and simulate CRAN graphs
#'
#' The state is a persistent graph of published package-major lines. Vertices
#' are born isolated before simultaneous edge edits. The release and edge
#' parameters have exactly the same interpretation as in [cran_updates()].
#' Node identities in the returned edge matrices are zero based.
#' A two-component fit draws one latent component for each complete update.
#' Its release and edit proposal rates apply to every releasing line in that
#' update, with nonempty conditioning performed within the selected component.
#' @param path CRAN application cache directory.
#' @param date UTC date whose preceding graph is required.
#' @param people Optional dated personnel data as in [cran_update_data()].
#'   The graph stores package personnel as of the forecast origin. Simulations
#'   hold these sets fixed, inherit them for new majors of existing families,
#'   and assign empty sets to newly simulated package families.
#' @param state A portable list containing names, birth days, family clocks and
#'   a zero-based directed edge matrix, or a native graph pointer.
#' @param fit A fitted CRAN update model.
#' @param times Update times in days since 1 January 2019.
#' @param checkpoints Times at which graph summaries are returned.
#' @param seed R random-number seed.
#' @return A portable state or simulated summaries and update sizes.
#' @name cran_simulation
NULL
.cran_check_fit<-function(fit){
 if(!is.character(fit$law)||length(fit$law)!=1L)stop('Invalid CRAN update-law fit')
 if(!is.null(fit$mixture)&&!fit$mixture %in% c('single','two_component'))stop('Invalid CRAN mixture specification')
 mixture<-identical(fit$mixture,'two_component')
 if(!is.null(fit$pair_proposals)&&!fit$pair_proposals %in% c('single','poisson','capped2'))stop('Invalid pair proposal law')
 if(isTRUE(fit$pair_proposals %in% c('poisson','capped2'))&&!fit$law %in% c('CS-pair-triangle','CS-pair-diamond','CS-pair-add','CS-pair-migration-add'))stop('Poisson addition pairs require a paired-addition law')
 expected<-switch(fit$law,PA=8L,`CS-1`=21L,`CS-1U`=16L,`CS-pair`=27L,TERGM=19L,`TERGM-logdegree`=23L,NA_integer_)
 if(fit$law %in% c('CS-1','CS-pair')&&identical(fit$weights,'exponential'))expected<-if(fit$law=='CS-1')15L else 20L
 if(identical(fit$proposal_means,'covariate'))expected<-expected+4L
 if(fit$law=='CS-pair-triangle'&&identical(fit$weights,'exponential')&&identical(fit$proposal_means,'covariate'))expected<-28L
 if(fit$law %in% c('CS-pair-diamond','CS-pair-add','CS-pair-migration-add')&&identical(fit$weights,'exponential')&&identical(fit$proposal_means,'covariate'))expected<-31L
 if(isTRUE(fit$term_specification %in% c("pattern_preferences","shared_consumer","pair_channels","motif_free"))&&fit$law %in% c("CS-pair","CS-pair-migration-add"))expected<-expected+3L
 if(isTRUE(fit$term_specification %in% c("shared_consumer","pair_channels","motif_free"))&&fit$law=="CS-pair-migration-add")expected<-expected+3L
 if(isTRUE(fit$migration_degree))expected<-expected+1L
 if(isTRUE(fit$shared_people)){
  if(!fit$law %in% c("CS-1","CS-pair","TERGM-logdegree") || !mixture ||
     !isTRUE(fit$term_specification %in% c("shared_consumer","motif_free")) ||
     (fit$law!="TERGM-logdegree" && (!identical(fit$weights,"exponential") ||
                                    !identical(fit$proposal_means,"covariate"))))
    stop("Shared-person simulations require a motif-free two-component CS or TERGM fit")
  expected<-expected+3L
 }
 if(isTRUE(fit$shared_people)||identical(fit$term_specification,"motif_free")){
  motif_terms<-switch(fit$law,`CS-1`=c(3:6,10:12),
    `CS-pair`=c(3:6,8L,10:12,16:19),`TERGM-logdegree`=c(2:7,12:16),integer())
  if(anyNA(fit$parameters[motif_terms])||any(fit$parameters[motif_terms]!=0))
    stop("Motif-free simulations require fixed-zero motif terms")
 }
 if(mixture)expected<-expected+4L
 if(is.na(expected)||length(fit$parameters)!=expected||any(!is.finite(fit$parameters)))stop('Incompatible edge parameters: refit with the current model specification')
 release_count<-if(mixture)19L else 15L
 if(length(fit$release_parameters)!=release_count||any(!is.finite(fit$release_parameters)))stop('Supply ',release_count,' finite release parameters')
 invisible(TRUE)
}
.cran_sampler_law <- function(fit){
 law <- if(fit$law == "PA") "PA" else if(isTRUE(fit$migration_degree))paste0(fit$law,"-consumers-olddegree") else
  if(isTRUE(fit$term_specification %in% c("shared_consumer","pair_channels","motif_free")))paste0(fit$law,"-consumers") else
  if(identical(fit$term_specification,"pattern_preferences"))paste0(fit$law,"-patterns") else fit$law
 # A missing flag identifies an archived one-pair fit, never a new Poisson fit.
 law<-if(identical(fit$pair_proposals,'poisson'))paste0(law,'-poisson')else if(identical(fit$pair_proposals,'capped2'))paste0(law,'-capped2')else law
 if(isTRUE(fit$shared_people))law<-paste0(law,'-people')
 if(identical(fit$mixture,'two_component'))paste0(law,'-mixture2')else law
}
.cran_people_at <- function(state,people,bin){
 people<-.cran_validate_people(people)
 families<-unique(sub("@[^@]+$","",state$names))
 known<-which(people$bin<bin&people$package%in%families)
 known<-known[order(people$bin[known])]
 known<-known[!duplicated(people$package[known],fromLast=TRUE)]
 state$people<-setNames(people$people[known],people$package[known])
 state
}
.cran_publications<-function(path){
 db<-normalizePath(file.path(path,'../cran_history.sqlite'),mustWork=TRUE)
 if(!nzchar(Sys.which('sqlite3')))stop('sqlite3 is needed to read the existing publication warehouse')
 dest<-tempfile(fileext='.csv');on.exit(unlink(dest))
 query<-"SELECT package,version,epoch FROM releases WHERE date_source='Date/Publication' AND epoch>=852076800 AND epoch<1767225600 AND version IS NOT NULL ORDER BY epoch,package,release_id"
 status<-system2('sqlite3',c('-readonly','-header','-csv',shQuote(db),shQuote(query)),stdout=dest)
 if(status!=0)stop('Cannot read publication warehouse')
 x<-utils::read.csv(dest,stringsAsFactors=FALSE);valid<-grepl('^[0-9]+([.-][0-9]+)*$',x$version);x<-x[valid,]
 x$name<-paste0(x$package,'@',as.integer(sub('^([0-9]+).*','\\1',x$version)));x[!duplicated(x$name),c('name','epoch')]
}
.cran_records<-function(path,publications=NULL){
 con<-file(file.path(path,'marks.jsonl'),'r');on.exit(close(con));records<-list();i<-0L
 repeat{lines<-readLines(con,n=1000L,warn=FALSE);if(!length(lines))break
  for(line in lines){r<-jsonlite::fromJSON(line,simplifyVector=FALSE);adds<-deletes<-matrix(character(),0,2);births<-as.character(unlist(r$births,use.names=FALSE));sources<-unique(c(births,vapply(r$updates,`[[`,'','package')))
   for(u in r$updates){before<-vapply(u$before,function(x)x[[2]],'');after<-vapply(u$after,function(x)x[[2]],'');a<-setdiff(after,before);d<-setdiff(before,after);if(length(a))adds<-rbind(adds,cbind(u$package,a));if(length(d))deletes<-rbind(deletes,cbind(u$package,d))}
   clocks<-if(is.null(publications))rep(r$time,length(births))else (publications$epoch[match(births,publications$name)]-1546300800)/86400
   if(anyNA(clocks))stop('Missing first-publication clock in replay')
   i<-i+1L;records[[i]]<-list(time=r$time,births=births,clocks=clocks,sources=sources,add=adds,delete=deletes)
  }
 };records
}
#' @rdname cran_simulation
#' @export
cran_graph_data<-function(path,date='2025-01-01',people=NULL){
 path<-normalizePath(path,mustWork=TRUE);when<-as.POSIXct(date,tz='UTC');cutoff<-as.numeric(when);epoch0<-1546300800
 if(length(cutoff)!=1L||!is.finite(cutoff))stop('date must identify one valid UTC forecast origin')
 publications<-.cran_publications(path)
 # Snapshot at the end of the preceding quarter; replay any subsequent updates.
 snapshot_dates<-sub('_nodes.csv.gz$','',list.files(file.path(path,'summaries'),pattern='_nodes.csv.gz$'))
 possible<-as.POSIXct(paste0(snapshot_dates,' 23:59:59'),tz='UTC')
 eligible<-which(possible<when);if(!length(eligible))stop('No preceding snapshot')
 snap<-snapshot_dates[eligible[which.max(possible[eligible])]]
 nodes<-utils::read.csv(gzfile(file.path(path,'summaries',paste0(snap,'_nodes.csv.gz'))),stringsAsFactors=FALSE)$node
 edges<-utils::read.csv(gzfile(file.path(path,'summaries',paste0(snap,'_edges.csv.gz'))),stringsAsFactors=FALSE)
 first<-publications$epoch[match(nodes,publications$name)];if(anyNA(first))stop('Missing publication time in snapshot')
 state<-list(names=nodes,birth=(floor(first/600)*600-epoch0)/86400,clock=(first-epoch0)/86400,edges=cbind(match(edges$source,nodes)-1L,match(edges$target,nodes)-1L))
 if(anyNA(state$edges))stop('Snapshot edges reference absent nodes')
 records<-.cran_records(path,publications);time<-(cutoff-epoch0)/86400;after<-(as.numeric(as.POSIXct(paste0(snap,' 23:59:59'),tz='UTC'))-epoch0)/86400
 g<-cran_graph_create(state);take<-vapply(records,function(r)r$time>after&&r$time<time,FALSE)
 if(any(take))cran_graph_replay(g,records[take])
 state<-cran_graph_snapshot(g)
 if(!is.null(people))state<-.cran_people_at(state,people,floor(cutoff/600))
 list(state=state,summary=cran_graph_summary(g),records=records,origin=time,date=date)
}
#' @rdname cran_simulation
#' @export
simulate_cran_updates<-function(state,fit,times,checkpoints=if(length(times)) max(times) else numeric(),seed=1L){
 .cran_check_fit(fit)
 if(!is.numeric(times)||any(!is.finite(times))||any(times<0)||anyDuplicated(times))stop('Update times must be finite, nonnegative and distinct')
 if(!is.numeric(checkpoints)||any(!is.finite(checkpoints))||any(checkpoints<0))stop('Checkpoints must be finite and nonnegative')
 set.seed(seed);g<-if(typeof(state)=='externalptr')cran_graph_create(cran_graph_snapshot(state))else cran_graph_create(state)
 times<-sort(times);checkpoints<-sort(checkpoints);snapshots<-vector('list',length(checkpoints));sizes<-matrix(0,length(times),5,dimnames=list(NULL,c('births','additions','deletions','releasing_lines','paired_proposals')));ci<-1L
 for(i in seq_along(times)){
  while(ci<=length(checkpoints)&&times[i]>=checkpoints[ci]){snapshots[[ci]]<-cran_graph_summary(g);ci<-ci+1L}
  u<-cran_graph_step(g,times[i],.cran_sampler_law(fit),fit$parameters,fit$release_parameters)
  sizes[i,]<-c(length(u$born),nrow(u$add),nrow(u$delete),length(u$sources),u$pair_proposals)
 }
 while(ci<=length(checkpoints)){snapshots[[ci]]<-cran_graph_summary(g);ci<-ci+1L}
 list(law=fit$law,seed=seed,times=times,checkpoints=checkpoints,snapshots=snapshots,sizes=sizes)
}
#' Simulate joint CRAN update timing and content
#'
#' Uses exact thinning for the two-exponential calendar Hawkes specification.
#' In M2 every past release retains its own summary at occurrence and subsequent
#' changes use the same clipping bound as the fitted timing model. Stops once
#' the requested first-creation
#' waits are known, or at the forecast horizon. Cohort forecasts always run
#' through the horizon. The event guard raises an error.
#' @param state Portable graph state at the forecast origin.
#' @param history Past update records, each with `time`, zero-based `families`,
#'   and their post-update log-degree `summary`. Supply history strictly before
#'   the origin, except for the release update of a migration cohort.
#'   For type-specific timing include integer `release_class` per line
#'   (0 existing line, 1 new package, 2 new major).
#' @param fit A CRAN update-law fit.
#' @param timing Named vector `mu`, `alpha_short`, `alpha_long`, `tau_short`,
#'   `tau_long`, `delta`, `gamma_B`, `gamma_Delta`, in physical units.
#'   Optionally include both `log_amp_new_package` and `log_amp_new_major`
#'   for M1 or M2. For separate effects on the two kernels supply all four
#'   coefficients using the same names with `_short` and `_long` suffixes.
#' @param calendar Daily baseline profile from 1 January 2019 onward.
#' @param origin Forecast origin in days since 1 January 2019.
#' @param horizon Forecast horizon in days. Zero returns the initialized intensity
#'   without generating updates.
#' @param model Timing amplitude specification, M0, M1 or M2.
#' @param seed R seed.
#' @param directed_waits Also record waits for new transitive triples, open
#'   diamonds and shared-consumer rectangles in the directed graph.
#' @param clip Nonnegative bound for later log-degree changes. NULL uses the
#'   clipping attribute attached by [fit_cran_timing()], or log(4) for a plain
#'   timing vector. Set this explicitly when supplying other fitted parameters.
#' @param migration_cohort Optional list with the new major-line name `target`
#'   and character vector of exposed consumer-family names `exposed`. Run
#'   through the horizon and record each exposed family's first migration.
#'   The target's release update at the origin is part of the history.
#' @export
simulate_cran_joint<-function(state,history,fit,timing,calendar,origin,horizon=14,model=c('M0','M1','M2'),seed=1L,directed_waits=FALSE,migration_cohort=NULL,clip=NULL){
 .cran_check_fit(fit)
 if(is.null(clip))clip<-attr(timing,'clip',exact=TRUE)
 if(is.null(clip))clip<-log(4)
 if(length(clip)!=1L||!is.finite(clip)||clip<0)stop('clip must be finite and nonnegative')
 if(length(origin)!=1L||!is.finite(origin)||origin<0||length(horizon)!=1L||!is.finite(horizon)||horizon<0)stop('Use finite nonnegative origin and horizon')
 if(!is.numeric(calendar)||any(!is.finite(calendar))||any(calendar<0)||length(calendar)<max(floor(origin)+1L,ceiling(origin+horizon)))stop('The nonnegative daily calendar must cover the forecast horizon')
 model<-match.arg(model);required<-c('mu','alpha_short','alpha_long','tau_short','tau_long','delta','gamma_B','gamma_Delta')
 shared<-any(.cran_type_amplitudes %in% names(timing))
 separate<-any(.cran_kernel_type_amplitudes %in% names(timing))
 if(shared&&separate)stop('Choose shared or kernel-specific release amplitudes')
 if(shared)required<-c(required,.cran_type_amplitudes)
 if(separate)required<-c(required,.cran_kernel_type_amplitudes)
 if(!all(required %in% names(timing)))stop('Missing timing parameters')
 if(any(!is.finite(timing[required])))stop('Invalid timing parameters')
 if(any(timing[c('mu','alpha_short','alpha_long','tau_short','tau_long')]<0)||any(timing[c('tau_short','tau_long')]<=0))stop('Invalid positive timing parameters')
 set.seed(seed);cran_joint_path(state,history,origin,horizon,.cran_sampler_law(fit),fit$parameters,fit$release_parameters,unname(timing[required]),calendar,model,directed_waits=directed_waits,migration_cohort=migration_cohort,clip=clip)
}
