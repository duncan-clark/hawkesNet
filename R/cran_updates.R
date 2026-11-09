#' CRAN complete-update models
#'
#' Exact complete-update likelihoods for the CRAN application. The release law
#' selects births and releasing major lines. Each existing releasing line must
#' change at least one edge; a newly born line may have no edges. All changes in
#' a recorded bundle are simultaneous. Addition and deletion coefficients,
#' including age decay, are separate. The manuscript uses `PA`, `CS-1`,
#' `CS-pair` and `TERGM-logdegree`. The default specification uses motif-free
#' weights and a two-component mixture. CS-pair proposes migrations jointly;
#' transitive triples, open diamonds and shared-consumer rectangles are
#' reserved for goodness-of-fit checks. The CS laws and TERGM also accept
#' dated shared-personnel covariates through `people`.
#'
#' For fitting the paper models, select `law` and supply `people` for the
#' CS and TERGM laws. See the installed `cran/README.md` guide for the cache
#' layout and the workflow from fit to forecast.
#'
#' @param path Directory containing the parameter-independent CRAN numerical
#'   caches (`marks`, `actor_selection`, `cs_popularity_20260924`). The importer
#'   reads the legacy NPY/NPZ format directly using R and C++; Python is not used.
#' @param law `PA`, `CS-1`, `CS-pair` or `TERGM-logdegree` for the paper.
#'   Earlier CRAN diagnostics also use `CS-1U`, `CS-pair-triangle`,
#'   `CS-pair-diamond`, `CS-pair-add`, `CS-pair-migration-add` and `TERGM`
#'   with explicit exploratory term specifications.
#' @param weights CS proposal weights, exponential by default. Logistic weights
#'   are retained for reproducing archived fits. Both TERGM variants use direct
#'   log odds. TERGM-logdegree replaces out-stars with the CS sender covariates.
#' @param proposal_means Covariate-dependent proposal means (default), or
#'   constant means for an ablation. Addition log means include indicators for
#'   new packages and new majors and log one plus source out-degree. Deletion
#'   log means include log one plus source out-degree. Existing lines are the
#'   reference class. The paired-migration proposal mean remains constant.
#' @param term_specification The default `motif_free` fixes motif coefficients to zero
#'   and is the manuscript specification. It supports PA, exponential CS-1,
#'   CS-pair and TERGM-logdegree. `legacy` reproduces archived fits. `reduced`
#'   fixes cyclic-triple coefficients to zero and, for laws with paired
#'   migrations, also fixes the single-addition migration coefficient to zero.
#'   Single edits and paired migrations remain admissible. In-star terms are
#'   retained. `pattern_preferences` removes all in-stars, replaces the cyclic
#'   addition term with a focal open-diamond completion count, and retains
#'   migration preferences in the single-only CS and TERGM laws. In paired CS
#'   laws, the single-addition migration coefficient is zero. Migration-pair
#'   transitivity, diamond, popularity and age coefficients are independent of
#'   the single-edit coefficients. Use `TERGM-logdegree` without out-stars.
#'   `shared_consumer` extends these patterns with a focal count of
#'   dependency pairs already used together by another package. It adds
#'   separate single-addition and migration-pair coefficients and, in
#'   `CS-pair-migration-add`, three release-class rates for a pair
#'   proposed through a common consumer.
#'   `pair_channels` retains CS-1 and one CS-pair-migration-add law; in the
#'   latter, single-edit migration, transitive, diamond and consumer coefficients
#'   are fixed at zero, with popularity, age, copying and rate regressions retained.
#' @param data An object returned by `cran_update_data`.
#' @param migration_degree Test an old-target in-degree interaction. In CS-1,
#'   add migration times log one plus old-target degree to addition log weights.
#'   If several older targets are present their in-degrees are summed. In the
#'   paired law, multiply each normalized migration-proposal intensity by
#'   `(1 + old_target_indegree)^eta`; this changes the total migration rate as
#'   well as pair selection. It is not another normalized popularity coefficient.
#' @param parameters Named or positional edge-law parameters. For mixture
#'   likelihoods, supply the joint vector in `fit$mixture_parameters` or the
#'   fitted object itself. For mixture fitting, an existing single-component
#'   or mixture fit can be supplied as a warm start.
#' @param pair_proposals Addition-pair count law: `capped2` for
#'   min(Poisson(M), 2), unrestricted `poisson`, or the archived
#'   `single` draw. NULL defaults to `capped2` for the paired `pair_channels`
#'   specification and `single` otherwise. Existing rate parameters are reused;
#'   they must be refitted when changing the count law. All pair candidates and
#'   weights use the pre-edit graph, and repeated edits are collapsed.
#' @param mixture `NULL` selects `two_component` with `motif_free` terms and
#'   `single` for exploratory term specifications. `two_component`
#'   selects one component for a complete update, with component-specific
#'   release and edit rates and shared target preferences. Both component laws
#'   are normalized separately, and their exact probabilities are mixed.
#'   Release and edge parameters are fitted jointly for this specification.
#' @param pair_max_states Computational budget for exact overlap enumeration
#'   with unrestricted `poisson` proposals (unused for `capped2`),
#'   default 262144 frontier states per connected component. Exceeding it raises
#'   an error before fitting; it never drops candidates, bundles or proposals.
#'   Increasing this budget can require substantial memory and computation.
#' @param people Optional dated personnel data with columns bin, package
#'   and people (integer person IDs per package). Available for the
#'   motif-free two-component exponential CS-1, CS-pair and TERGM-logdegree fits. Bin numbers
#'   count UTC ten-minute intervals from 1 January 1970.
#' @param window Estimation window: 2020--2023 (`train`) or 2020--2025 (`full`).
#' @param starts Number of deterministic optimisation starts, at least two.
#' @param control Additional `optim` controls. `polish_tol` optionally sets
#'   the stopping tolerance for the per-update projected score in refinement.
#' @param per_update Return individual update contributions and gradients.
#' @param checkpoint Optional RDS path for retaining optimisation iterates.
#' @return Data, likelihood contributions, or a fitted model with predictive
#'   likelihoods, boundary flags and optimisation diagnostics.
#' @name cran_updates
#' @useDynLib hawkesNet, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

.cran_npz <- function(path) {
  directory <- tempfile("cran-cache-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  files <- utils::unzip(path, exdir = directory)
  ans <- lapply(files, cran_read_npy)
  names(ans) <- sub("\\.npy$", "", basename(files)); ans
}
.cran_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
.cran_parameters <- function(parameters, expected, label="parameters") {
  if(!is.numeric(parameters) || length(parameters)!=length(expected) || any(!is.finite(parameters)))
    stop("Invalid ",label,": expected ",length(expected)," finite values")
  if(!is.null(names(parameters))) {
    if(anyDuplicated(names(parameters)) || !setequal(names(parameters),expected))
      stop("Names of ",label," do not match this model")
    parameters <- parameters[expected]
  }
  setNames(as.numeric(parameters),expected)
}

.cran_sum <- function(group, x, n) {
  if (is.null(dim(x))) return(drop(cran_group_sum(group, matrix(x, ncol=1), n)))
  cran_group_sum(group, x, n)
}
.cran_norm <- function(group, x, n) {
  # Group IDs are zero based, as in the stored candidate arrays.
  mx <- rep(-Inf,n); v <- tapply(x, group, max); mx[as.integer(names(v))+1L] <- v
  w <- exp(x - mx[group+1L]); s <- .cran_sum(group,w,n)
  list(logsum=ifelse(s>0, mx+log(s), 0), probability=w/s[group+1L])
}
.cran_features <- function(d,o,tr=0,cy=0,pop=d,reverse=NULL,flags=NULL,delete=FALSE,undir=FALSE,pa=FALSE,exponential=FALSE) {
  d<-drop(d);o<-drop(o);pop<-drop(pop);n<-length(d)
  if(pa) return(if(delete) matrix(log1p(d),ncol=1) else cbind(log1p(d),if(is.null(flags)) matrix(0,n,2) else flags))
  v<-d-as.integer(delete);u<-o-as.integer(delete);sgn<-if(delete) -1 else 1
  motif <- if(undir) cbind((v+u)/100,(v*(v-1)+u*(u-1))/20000,rep_len(tr,n)/10) else
    cbind(v/100,u/100,v*(v-1)/20000,u*(u-1)/20000,rep_len(tr,n)/10,rep_len(cy,n)/10)
  motif<-sgn*motif
  if(!is.null(reverse)) motif[as.logical(reverse),]<-0
  if(delete && !undir) {if(any(motif[,6]!=0)) stop("Cyclic deletion is not identically zero");motif<-motif[,-6,drop=FALSE]}
  x<-cbind(1,motif,if(delete) -log(pmax(pop,1)) else log1p(pop))
  x<-if(delete)x else cbind(x,if(is.null(flags)) matrix(0,n,2) else flags)
  if(exponential&&!undir)x<-x[,if(delete)c(2,4,6,7)else c(2,4,6:10),drop=FALSE]
  x
}
.cran_directed_pairs <- function(m) {
  file<-file.path(m$path,'hawkesNet_cache','directed_pairs.rds')
  if(file.exists(file))return(readRDS(file))
  d<-cran_graph_data(m$path,'2020-01-01')
  updates<-d$records[vapply(d$records,function(r)r$time>=d$origin,FALSE)]
  groups<-.cran_json(file.path(m$path,'cs2_paired_20260924/groups.json'))
  sources<-vapply(groups,`[[`,'','source');batches<-vapply(groups,`[[`,0L,'batch')
  stopifnot(length(updates)==m$n,identical(as.integer(m$gb),batches))
  x<-cran_directed_pair_cache(d$state,updates,sources,batches)
  old<-.cran_npz(file.path(m$path,'cs2_paired_20260924/pairs.npz'))$values
  # Candidate identities and the old undirected scores are independently
  # checked before adopting the new directed scores.
  stopifnot(nrow(x$values)==nrow(old),max(abs(x$values[,c(1,11:13,8:10)]-old))<1e-9)
  obs<-x$observed
  key<-function(group,features,birth)apply(cbind(group,features,birth),1,function(z)paste(sprintf('%.12g',z),collapse='/'))
  ak<-key(m$ag,m$XA,m$a$added_birth)
  dk<-key(m$sg,m$XD[m$sel,,drop=FALSE],m$a$existing_birth[m$sel])
  xa<-.cran_features(obs[,3],obs[,4],obs[,5],obs[,6],flags=cbind(rep(1,nrow(obs)),rep(0,nrow(obs))),exponential=m$exponential)
  xd<-.cran_features(obs[,8],obs[,9],obs[,10],obs[,11],delete=TRUE,exponential=m$exponential)
  match_all<-function(want,available){pool<-split(seq_along(available),available);vapply(want,function(k){v<-pool[[k]];if(!length(v))stop('Migration endpoint does not match the single-edit cache');pool[[k]]<<-v[-1L];v[1L]},0L)}
  x$add_index<-match_all(key(obs[,2],xa,obs[,7]),ak)
  x$delete_index<-length(m$ag)+match_all(key(obs[,2],xd,obs[,12]),dk)
  stopifnot(!anyDuplicated(c(x$add_index,x$delete_index)))
  x$schema<-'directed-migration-v1'
  dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE);saveRDS(x,file);x
}
.cran_pattern_preferences <- function(m, sources) {
  folder <- file.path(m$path,"hawkesNet_cache")
  consumers <- m$term_specification %in% c("shared_consumer", "pair_channels", "motif_free")
  stem <- if(consumers)"shared_consumer_v1" else "single_pattern_preferences_v2"
  file <- file.path(folder,paste0(stem,".rds"))
  raw <- file.path(folder,paste0(stem,".bin"))
  if(file.exists(file) && file.exists(raw)) x <- readRDS(file) else {
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    d <- cran_graph_data(m$path,"2020-01-01")
    records <- d$records[vapply(d$records,function(r)r$time>=d$origin,FALSE)]
    stopifnot(length(records)==m$n)
    temporary <- tempfile(stem,tmpdir=folder)
    on.exit(unlink(c(temporary,paste0(temporary,".rds"))))
    x <- cran_pattern_preference_cache(d$state,records,sources,as.integer(m$gb),temporary,consumers)
    saveRDS(x,paste0(temporary,".rds"))
    stopifnot(file.rename(temporary,raw),file.rename(paste0(temporary,".rds"),file))
  }
  stopifnot(identical(x$schema,if(consumers)"shared-consumer-v1" else "single-pattern-preferences-v2"),
            file.info(raw)$size == x$candidate_records*8)
  a <- m$a; e <- x$observed
  key <- function(g,d,o,tr,cy,b,mig,copy)paste(g,d,o,tr,cy,round(b*365.25*144),mig,copy,sep="/")
  oldkey <- key(m$ag,drop(a$added_degree),drop(a$added_outdegree),drop(a$added_triangle),
                drop(a$added_cycle),a$added_birth,m$hist$selected_flags[,1],m$hist$selected_flags[,2])
  newkey <- key(e[,1],e[,2],e[,3],e[,4],e[,5],e[,6],e[,7],e[,8])
  if(!identical(sort(oldkey),sort(newkey)))stop("Observed diamond features disagree with graph replay")
  # Use replayed identities for the observed additions. This also fixes the
  # mapping of overlapping migration/addition pairs when old feature rows tie.
  m$ag <- a$added_group <- e[,1]
  a$added_degree <- matrix(e[,2],ncol=1)
  a$added_outdegree <- matrix(e[,3],ncol=1)
  a$added_triangle <- matrix(e[,4],ncol=1)
  a$added_cycle <- matrix(e[,9],ncol=1)
  a$added_birth <- e[,6]
  m$a <- a
  hist <- x$history; hist$selected_flags <- e[,7:8,drop=FALSE];m$hist <- hist
  if(!is.null(m$oi)) {
    ix <- match(m$oi-1L,x$pair_observed[,1])
    if(anyNA(ix))stop("Observed migration identities disagree with replay")
    m$ai <- x$pair_observed[ix,2]+1L
  }
  if(!is.null(m$triangle)) {
    pack <- function(v)list(group=v[,1],a=v[,2]+1L,b=v[,3]+1L,log_weight=log(v[,4]))
    triangle <- c(pack(x$triangles),list(copied=x$triangles[,5],mass=m$triangle$mass,
                                        normaliser_checks=m$triangle$normaliser_checks))
    if(!is.null(m$triangle$diamond))triangle$diamond <- c(pack(x$diamonds),list(mass=m$triangle$diamond$mass))
    if(consumers)triangle$consumer <- c(pack(x$consumer_pairs),list(mass=x$consumer_mass))
    m$triangle <- triangle
  }
  m$XA <- .cran_features(a$added_degree,a$added_outdegree,a$added_triangle,a$added_cycle,
                          flags=hist$selected_flags,exponential=m$exponential)
  h <- hist$features
  sender_degree <- if(consumers)drop(a$group_outdegree)[hist$group+1L] else h[,3]
  m$XH0 <- .cran_features(h[,1],sender_degree,h[,4],h[,5],exponential=m$exponential)
  if(consumers){
    # Reuse the retired in-three-star slot in the native parameter layout.
    slot <- if(m$exponential)2L else 4L
    m$XA[,slot] <- e[,11]/10
    m$XH0[,slot] <- h[,3]/10
    m$XE[,slot] <- m$XS[,slot] <- 0
    m$names[m$names %in% c("addition_in3star","addition_in_3star")] <- "addition_shared_consumer"
  }
  m$XH <- m$XH0;m$XH[,(ncol(m$XH)-1):ncol(m$XH)] <- hist$flags
  m$cache <- cran_native_cache(a,x$index,x$codes,raw,hist,m$time,FALSE)
  m$exp_cache <- NULL
  m$names[m$names=="addition_cyclic"] <- "addition_diamond"
  if(!is.null(m$ip)) {
    stopifnot(length(x$pair_diamond)==nrow(m$PX))
    m$PX[,ncol(m$PX)] <- x$pair_diamond
    j <- match("pair_cyclic",m$names);m$names[j] <- "pair_diamond"
    m$lower[j] <- -12;m$upper[j] <- 12
    if(consumers){
      stopifnot(length(x$pair_consumer)==nrow(m$PX))
      j <- match("pair_in3star",m$names)
      m$PX[,match(j,m$ip)] <- x$pair_consumer
      m$names[j] <- "pair_shared_consumer"
      m$lower[j] <- -12;m$upper[j] <- 12
    }
    m$pair_preferences <- length(m$names)+seq_len(3L)
    m$names <- c(m$names,"pair_popularity_add","pair_popularity_delete","pair_age_decay_add")
    m$lower <- c(m$lower,-12,-12,0);m$upper <- c(m$upper,12,12,2)
    # The single-only model retains migration preference. In paired models
    # the migration stream supplies it, as requested for this comparison.
    j <- match("migration",m$names);m$lower[j] <- m$upper[j] <- 0
  }
  if(m$law %in% c("TERGM","TERGM-logdegree"))
    m$obs[,7] <- .cran_sum(m$ag,a$added_cycle/10,m$ng)
  if(consumers && m$law=="TERGM-logdegree")
    m$obs[,4] <- .cran_sum(m$ag,e[,11]/10,m$ng)
  if(consumers && !is.null(m$triangle)){
    m$meanConsumer <- length(m$names)+1:3
    m$names <- c(m$names,paste0("log_m_consumer_",c("existing","new_package","new_major")))
    m$lower <- c(m$lower,rep(log(1e-8),3));m$upper <- c(m$upper,rep(log(10),3))
  }
  removed <- grep("(in2star|in3star|in_[23]star)$",m$names)
  m$lower[removed] <- m$upper[removed] <- 0
  m$pattern_cache <- list(schema=x$schema,candidate_records=x$candidate_records,max_diamond=x$max_diamond)
  if(consumers)m$pattern_cache$consumer_checks <- x$consumer_checks
  invisible(m)
}

.cran_triangle_pairs <- function(m,sources) {
  diamond<-m$law %in% c("CS-pair-diamond","CS-pair-add","CS-pair-migration-add")
  file<-file.path(m$path,"hawkesNet_cache",if(diamond)"triangle_diamond_pairs_v2.rds"else"triangle_pairs_v1.rds")
  if(file.exists(file)) x<-readRDS(file) else {
    d<-cran_graph_data(m$path,"2020-01-01")
    take<-vapply(d$records,function(r)r$time>=d$origin,FALSE)
    stopifnot(sum(take)==m$n)
    x<-cran_triangle_pair_cache(d$state,d$records[take],sources,as.integer(m$gb),diamond)
    dir.create(dirname(file),showWarnings=FALSE,recursive=TRUE);saveRDS(x,file)
  }
  e<-x$edits; a<-m$a
  stopifnot(nrow(e)==length(m$ag),nrow(x$mass)==m$ng)
  # Match rates using exact sufficient statistics; legacy NPZ stores no target
  # names. Pin compatible migration endpoints in the same order as their cache.
  key<-function(g,d,o,tr,cy,b,mig,copy)paste(g,d,o,tr,cy,round(b*365.25*144),mig,copy,sep="/")
  k<-key(e[,1],e[,2],e[,3],e[,4],e[,5],e[,6],e[,7],e[,8])
  ka<-key(m$ag,drop(a$added_degree),drop(a$added_outdegree),drop(a$added_triangle),
          drop(a$added_cycle),a$added_birth,m$hist$selected_flags[,1],m$hist$selected_flags[,2])
  pin<-which(e[,9]>=0);ix<-integer(nrow(e))
  stopifnot(identical(sort(as.integer(e[pin,9])),seq_along(m$ai)-1L))
  ix[pin]<-m$ai[e[pin,9]+1L]; stopifnot(all(k[pin]==ka[ix[pin]]))
  left<-setdiff(seq_along(k),pin); right<-setdiff(seq_along(ka),ix[pin]);
  u<-left[order(k[left],method="radix")];v<-right[order(ka[right],method="radix")]
  if(!identical(k[u],ka[v]))stop("Paired-addition endpoints disagree with observed single-edit features")
  ix[u]<-v
  pack<-function(obs){
    a<-ix[obs[,2]+1L];b<-ix[obs[,3]+1L]
    stopifnot(all(m$ag[a]==obs[,1]),all(m$ag[b]==obs[,1]))
    list(group=obs[,1],a=a,b=b,log_weight=log(obs[,4]))
  }
  result<-c(pack(x$observed),list(copied=x$observed[,5],mass=x$mass,normaliser_checks=x$normaliser_checks))
  if(diamond)result$diamond<-c(pack(x$observed_diamonds),list(mass=x$diamond_mass))
  result
}
#' @rdname cran_updates
#' @export
cran_update_data <- function(path,
                             law = c("PA", "CS-1", "CS-pair", "TERGM-logdegree",
                                     "CS-1U", "CS-pair-triangle", "CS-pair-diamond",
                                     "CS-pair-add", "CS-pair-migration-add", "TERGM"),
                             weights = c("exponential", "logistic"),
                             proposal_means = c("covariate", "constant"),
                             term_specification = c("motif_free", "legacy", "reduced",
                                                    "pattern_preferences", "shared_consumer",
                                                    "pair_channels"),
                             migration_degree = FALSE, pair_proposals = NULL,
                             pair_max_states = 262144, mixture = NULL, people = NULL) {
  law<-if(identical(law,"CS-pair-migration"))"CS-pair"else match.arg(law);path<-normalizePath(path,mustWork=TRUE)
  m<-new.env(parent=emptyenv());m$path<-path;m$law<-law;m$pa<-law=="PA";m$undir<-law=="CS-1U"
  m$term_specification <- match.arg(term_specification)
  if (is.null(mixture)) {
    mixture <- if (m$term_specification == "motif_free") "two_component" else "single"
  }
  m$mixture <- match.arg(mixture, c("single", "two_component"))
  tergm <- law %in% c("TERGM", "TERGM-logdegree")
  m$weights<-if(m$undir&&missing(weights))"logistic"else match.arg(weights)
  m$exponential<-!m$pa&&!m$undir&&!tergm&&m$weights=="exponential"
  m$proposal_means<-if(tergm)"direct"else match.arg(proposal_means)
  if(m$term_specification == "motif_free" && !law %in% c("PA","CS-1","CS-pair","TERGM-logdegree"))
    stop("The motif-free specification supports PA, CS-1, CS-pair and TERGM-logdegree")
  if(!is.null(people)) {
    if(!law %in% c("CS-1","CS-pair","TERGM-logdegree") ||
       !m$term_specification %in% c("shared_consumer","motif_free") ||
       m$mixture!="two_component" || (!tergm && (!m$exponential || m$proposal_means!="covariate")))
      stop("Personnel terms require CS-1, CS-pair or TERGM-logdegree with motif_free terms, two_component mixture and covariate exponential CS weights")
    people <- .cran_validate_people(people)
  }
  paired_additions<-law %in% c("CS-pair-triangle","CS-pair-diamond","CS-pair-add","CS-pair-migration-add")
  if(is.null(pair_proposals))pair_proposals<-if(paired_additions&&m$term_specification=='pair_channels')'capped2'else'single'
  m$pair_proposals <- match.arg(pair_proposals,c('capped2','poisson','single'))
  if(m$pair_proposals!='single'&&!paired_additions)stop('Poisson addition pairs require a paired-addition law')
  if(length(pair_max_states)!=1L||!is.finite(pair_max_states)||pair_max_states<4)stop('pair_max_states must be finite and at least four')
  m$pair_max_states <- pair_max_states
  if(paired_additions&&(!m$exponential||m$proposal_means!="covariate"))stop("Paired additions use the current exponential CS-pair specification")
  if(m$undir&&m$weights=="exponential")stop("Archived CS-1U requires weights='logistic'; use directed CS-1 for exponential weights")
  folder<-if(m$undir) file.path(path,"cs_popularity_20260924/undirected") else file.path(path,"marks")
  a<-m$a<-.cran_npz(file.path(folder,"broad_arrays.npz"))
  hist<-m$hist<-.cran_npz(file.path(path,"cs_popularity_20260924",if(m$undir) "undirected_features.npz" else "directed_features.npz"))
  rec<-.cran_json(file.path(path,"marks/batches.json"));m$year<-vapply(rec,`[[`,0L,"year")
  m$time<-(vapply(rec,`[[`,0,"epoch")-1577836800)/(365.25*86400)
  m$n<-length(rec);m$gb<-a$group_batch;m$ng<-length(m$gb);m$old<-abs(a$group_source_birth-m$time[m$gb+1L])>1e-12
  m$eg<-a$existing_group;m$ag<-a$added_group;m$sel<-which(a$existing_outcome==1);m$sg<-m$eg[m$sel]
  if(any(!a$existing_outcome %in% 0:1))stop("Invalid edge outcomes")
  motif<-if(m$undir) c("two_stars","three_stars","triangles") else c("in2star","out2star","in3star","out3star","transitive","cyclic")
  m$namesA<-if(m$pa)c("pa_add","migration","copy") else c(if(!m$exponential)"addition_intercept",paste0("addition_",if(m$exponential)setdiff(motif,c("out2star","out3star"))else motif),"popularity_add","migration","copy")
  m$namesD<-if(m$pa)"pa_delete" else c(if(!m$exponential)"deletion_intercept",paste0("deletion_",setdiff(motif,if(m$exponential)c("cyclic","out2star","out3star")else"cyclic")),"popularity_delete")
  m$na<-length(m$namesA);m$nd<-length(m$namesD);m$ia<-2+seq_len(m$na);m$id<-2+m$na+seq_len(m$nd)
  m$ageA<-3+m$na+m$nd;m$ageD<-m$ageA+1L
  m$names<-c("log_m_add","log_m_delete",m$namesA,m$namesD,"age_decay_add","age_decay_delete")
  m$lower<-c(rep(-12,2+m$na+m$nd),0,0);m$upper<-c(4,4,rep(12,m$na+m$nd),2,2)
  f<-function(...) .cran_features(...,undir=m$undir,pa=m$pa,exponential=m$exponential)
  m$XA<-f(a$added_degree,a$added_outdegree,a$added_triangle,a$added_cycle,pop=if(is.null(a$added_popdegree))a$added_degree else a$added_popdegree,reverse=a$added_reverse,flags=hist$selected_flags)
  m$XD<-f(a$existing_degree,a$existing_outdegree,a$existing_triangle,a$existing_cycle,pop=if(is.null(a$existing_popdegree))a$existing_degree else a$existing_popdegree,reverse=a$existing_reverse,delete=TRUE)
  m$XE<-f(a$existing_degree,a$existing_outdegree,pop=if(is.null(a$existing_popdegree))a$existing_degree else a$existing_popdegree)
  m$XS<-f(a$group_source_degree,a$group_outdegree,pop=if(is.null(a$group_source_popdegree))a$group_source_degree else a$group_source_popdegree)
  h<-hist$features;m$XH0<-f(h[,1],h[,3],h[,4],h[,5],pop=if(ncol(h)==7)h[,7] else h[,1],reverse=if(ncol(h)==7)h[,6] else NULL)
  m$XH<-m$XH0;m$XH[,(ncol(m$XH)-1):ncol(m$XH)]<-hist$flags
  m$exposure<-cbind(a$group_has_addition,a$group_has_existing)
  if(!m$pa){
    grouped<-file.path(path,"cs_popularity_20260924",if(m$undir)"undirected_grouped"else"directed_grouped")
    index<-if(file.exists(file.path(grouped,"index.npz"))).cran_npz(file.path(grouped,"index.npz"))else NULL
    m$cache<-cran_native_cache(a,.cran_npz(file.path(folder,"motif_index.npz")),cran_read_npy(file.path(folder,"motif_codes.npy")),file.path(folder,"motif_candidates.bin"),hist,m$time,m$undir,index,file.path(grouped,"birth_ids.bin"))
  }
  if(law=="CS-pair"||paired_additions) {
    pairs<-.cran_directed_pairs(m);v<-pairs$values
    m$pg<-v[,1];m$PX<-if(m$exponential)v[,c(2,4,6,7),drop=FALSE]else cbind(1,v[,c(2,4,6,7)]);m$Pb<-v[,8];m$Ppop<-v[,9:10]
    m$hasP<-.cran_sum(m$pg,rep(1,nrow(v)),m$ng)>0
    m$ip<-length(m$names)+seq_len(ncol(m$PX));m$meanP<-length(m$names)+ncol(m$PX)+1L
    m$names<-c(m$names,if(!m$exponential)"pair_intercept","pair_in2star","pair_in3star","pair_transitive","pair_cyclic","log_m_migration")
    m$lower<-c(m$lower,if(!m$exponential)-24,rep(-12,4),log(1e-8));m$upper<-c(m$upper,if(!m$exponential)24,rep(12,4),4)
    zero<-which(colSums(abs(m$PX))==0);m$lower[m$ip[zero]]<-m$upper[m$ip[zero]]<-0
    m$oi<-pairs$observed[,1]+1L;m$ai<-pairs$add_index;m$di<-pairs$delete_index;m$og<-m$pg[m$oi]
    if(anyDuplicated(c(m$ai,m$di)))stop("Overlapping observed pairs require full coverage enumeration")
  }
  m$rateA<-1L;m$rateD<-2L;m$ZA<-m$ZD<-matrix(1,m$ng,1)
  if(m$proposal_means=="covariate"||law=="TERGM-logdegree"){
    groups<-.cran_json(file.path(path,"cs2_paired_20260924/groups.json"));sources<-vapply(groups,`[[`,"","source")
    stopifnot(length(sources)==m$ng,identical(as.integer(m$gb),vapply(groups,`[[`,0L,"batch")))
    pub<-.cran_publications(path);family<-sub("@[^@]+$","",pub$name);first<-tapply(pub$epoch,family,min)
    source_family<-sub("@[^@]+$","",sources);born<-floor(pub$epoch[match(sources,pub$name)]/600);family_born<-floor(first[source_family]/600)
    if(anyNA(born)||anyNA(family_born))stop("Missing release-class publication evidence")
    newpackage<-!m$old&born==family_born;newmajor<-!m$old&!newpackage
    m$ZA<-cbind(1,newpackage,newmajor,log1p(drop(a$group_outdegree)));m$ZD<-cbind(1,log1p(drop(a$group_outdegree)))
    if(!tergm){
      ix<-length(m$names)+seq_len(4);m$rateA<-c(1L,ix[1:3]);m$rateD<-c(2L,ix[4])
      m$names<-c(m$names,"mean_add_new_package","mean_add_new_major","mean_add_log_outdegree","mean_delete_log_outdegree")
      m$lower<-c(m$lower,-6,-6,-3,-3);m$upper<-c(m$upper,6,6,3,3)
    }
    m$release_class<-ifelse(m$old,"existing",ifelse(newpackage,"new_package","new_major"))
  }
  if(paired_additions){
    m$triangle<-.cran_triangle_pairs(m,sources)
    m$Tclass<-match(m$release_class,c("existing","new_package","new_major"))
    m$meanT<-length(m$names)+1:3;m$copyT<-length(m$names)+4L
    m$names<-c(m$names,"log_m_triangle_existing","log_m_triangle_new_package","log_m_triangle_new_major","triangle_copy_both")
    m$lower<-c(m$lower,rep(log(1e-8),3),-12);m$upper<-c(m$upper,rep(log(10),3),12)
    if(law %in% c("CS-pair-diamond","CS-pair-add","CS-pair-migration-add")){
      m$meanDiamond<-length(m$names)+1:3
      m$names<-c(m$names,paste0("log_m_diamond_",c("existing","new_package","new_major")))
      m$lower<-c(m$lower,rep(log(1e-8),3));m$upper<-c(m$upper,rep(log(10),3))
    }
    # Keep inactive entries fixed in the shared simulator parameter layout.
    if(law %in% c("CS-pair-add","CS-pair-migration-add"))m$lower[m$copyT]<-m$upper[m$copyT]<-0
    if(law=="CS-pair-add"){
      inactive<-c(m$ip,m$meanP)
      m$lower[inactive]<-m$upper[inactive]<-0
    }
  }
  if(tergm) {
    m$names<-c(paste0("addition_",c("intercept","in_2star","out_2star","in_3star","out_3star","transitive","cyclic","popularity","migration","copy")),paste0("deletion_",c("intercept","in_2star","out_2star","in_3star","out_3star","transitive","popularity")),"addition_age_decay","deletion_age_decay")
    m$lower<-c(-20,rep(-12,9),-20,rep(-12,6),0,0);m$upper<-c(5,rep(12,9),5,rep(12,6),2,2)
    m$obs<-cbind(.cran_sum(m$ag,m$XA,m$ng),.cran_sum(m$sg,m$XD[m$sel,,drop=FALSE],m$ng),-.cran_sum(m$ag,m$time[m$gb[m$ag+1L]+1L]-a$added_birth,m$ng),-.cran_sum(m$sg,m$time[m$gb[m$sg+1L]+1L]-a$existing_birth[m$sel],m$ng))
    if(law=="TERGM-logdegree"){
      out_stars <- grep("_out_[23]star$", m$names)
      m$lower[out_stars] <- m$upper[out_stars] <- 0
      m$names <- c(m$names, "addition_new_package", "addition_new_major",
                   "addition_log_outdegree", "deletion_log_outdegree")
      m$lower <- c(m$lower, -6, -6, -3, -3)
      m$upper <- c(m$upper, 6, 6, 3, 3)
      add_counts <- .cran_sum(m$ag, rep(1, length(m$ag)), m$ng)
      delete_counts <- .cran_sum(m$sg, rep(1, length(m$sg)), m$ng)
      m$obs <- cbind(m$obs, m$ZA[, -1, drop = FALSE] * add_counts, m$ZD[, 2] * delete_counts)
    }
  }
  if(m$term_specification == "reduced") {
    # Retain zero slots so archived parameter vectors and the exact native
    # likelihood/sampler keep their shared layout; these are not fitted terms.
    removed <- grep("(^|_)cyclic$", m$names)
    if(law %in% c("CS-pair", "CS-pair-triangle", "CS-pair-diamond", "CS-pair-migration-add"))
      removed <- c(removed, match("migration", m$names))
    m$lower[removed] <- m$upper[removed] <- 0
  }
  if(m$term_specification %in% c("pattern_preferences","shared_consumer","pair_channels","motif_free") && !m$pa) {
    if(!law %in% c("CS-1","CS-pair","CS-pair-migration-add","TERGM-logdegree"))
      stop("Pattern preferences support CS-1, both retained paired laws and TERGM-logdegree")
    if(!tergm && (!m$exponential || m$proposal_means!="covariate"))
      stop("Pattern preferences require exponential weights and covariate proposal means")
    .cran_pattern_preferences(m,sources)
    if(m$term_specification == "pair_channels") {
      if(!law %in% c("CS-1","CS-pair-migration-add"))
        stop("The two-law comparison supports CS-1 and CS-pair-migration-add")
      if(law == "CS-pair-migration-add") {
        background <- match(c("migration", "addition_transitive", "addition_diamond",
                              "addition_shared_consumer", "deletion_transitive"), m$names)
        stopifnot(!anyNA(background))
        m$lower[background] <- m$upper[background] <- 0
      }
    }
  }
  if(m$term_specification == "motif_free") {
    .cran_fix_motif_terms(m)
  }
  m$migration_degree <- isTRUE(migration_degree)
  if(m$migration_degree) {
    if(m$term_specification!="pair_channels")stop("Degree extension requires pair_channels specification")
    m$degree_index <- length(m$names)+1L
    m$names <- c(m$names,"migration_old_logdegree")
    m$lower <- c(m$lower,-3);m$upper <- c(m$upper,3)
    if(law=="CS-1") {
      file <- file.path(path,"hawkesNet_cache/migration_old_degree_v1.rds")
      if(file.exists(file))degree <- readRDS(file) else {
        d <- cran_graph_data(path,"2020-01-01")
        records <- d$records[vapply(d$records,function(r)r$time>=d$origin,FALSE)]
        degree <- cran_migration_degree_cache(d$state,records,sources,as.integer(m$gb))
        saveRDS(degree,file)
      }
      stopifnot(length(degree$history)==nrow(m$XH),length(degree$observed)==nrow(m$XA),
                identical(degree$history>0,as.logical(m$hist$flags[,1])),
                identical(degree$observed>0,as.logical(m$hist$selected_flags[,1])))
      m$degree_history <- degree$history;m$degree_observed <- degree$observed
    } else m$degree_pair <- log1p(exp(-m$Ppop[,2]))
  }
  if(m$pair_proposals %in% c('poisson','capped2')){
    x<-m$triangle;parts<-list(x,x$diamond,x$consumer)
    group<-unlist(lapply(parts,`[[`,'group'));a<-unlist(lapply(parts,`[[`,'a'));b<-unlist(lapply(parts,`[[`,'b'))
    m$pair_union_plan<-cran_pair_union_plan(as.integer(group),as.integer(a-1L),as.integer(b-1L),m$ng,length(m$ag))
    m$pair_union_complexity<-cran_pair_union_complexity(m$pair_union_plan)
  }
  m$release_data<-.cran_npz(file.path(path,"actor_selection/risk_sets.npz"))
  m$families<-length(.cran_json(file.path(path,"actor_selection/families.json")))
  if(!is.null(people)) .cran_attach_people(m,people,sources)
  class(m)<-"cran_update_data";m
}
.cran_weights<-function(x,b,beta,eta,pa){z<-drop(x %*% beta);if(pa)list(logw=z+eta*b,score=x)else list(logw=stats::plogis(z,log.p=TRUE)+eta*b,score=x*stats::plogis(-z))}
.cran_old_fit<-function(data,window){p<-file.path(data$path,"current_fits",if(window=="train")"chronological/JOINT_FITS.json" else "JOINT_FITS.json");kind<-switch(data$law,"CS-1"="CS1","CS-1U"="CS1U","CS-pair"="CS1","CS-pair-triangle"="CS1","CS-pair-diamond"="CS1","CS-pair-add"="CS1","CS-pair-migration-add"="CS1","TERGM"="CS1","TERGM-logdegree"="CS1",data$law);if(file.exists(p)).cran_json(p)[[paste0("all_two_M0__",kind)]]else NULL}
.cran_old_start<-function(m,window){old<-.cran_old_fit(m,window);p<-setNames(unlist(old$parameters),unlist(old$parameter_names));if(m$law %in% c("TERGM","TERGM-logdegree"))return(c(-10,rep(0,6),1,2,7,-1,rep(0,5),.1,.2,.1, if(m$law=="TERGM-logdegree")rep(0,4), rep(0,length(m$names)-if(m$law=="TERGM-logdegree")23L else 19L)))
 out<-numeric(length(m$names))
 if(is.null(old)) {
   out[1:2]<-log(c(1,.2))
   out[c(m$ageA,m$ageD)]<-.1
   if(!is.null(m$meanP))out[m$meanP]<-log(.1)
   if(!is.null(m$meanT))out[m$meanT]<-log(.1)
   if(!is.null(m$meanDiamond))out[m$meanDiamond]<-log(.1)
   if(!is.null(m$meanConsumer))out[m$meanConsumer]<-log(.1)
   return(pmin(m$upper,pmax(m$lower,out)))
 }
 out[1:2]<-p[c("log_m_add","log_m_delete")]
 if(m$pa){out[m$ia]<-p[c("pa_add","migration","copy")];out[m$id]<-p["pa_delete"]}else{terms<-if(m$undir)c("star2","star3","triangle") else if(m$exponential)c("in2star","in3star","transitive","cycle")else c("in2star","out2star","in3star","out3star","transitive","cycle");out[m$ia]<-p[c(if(!m$exponential)"alpha",terms,"popularity_add","migration","copy")];out[m$id]<-p[c(if(!m$exponential)"alpha",setdiff(terms,"cycle"),"popularity_delete")]}
 out[c(m$ageA,m$ageD)]<-p["age_decay"]
 if(!is.null(m$meanP)&&m$law!="CS-pair-add"){out[m$ip]<-c(if(!m$exponential)p["alpha"]+p["migration"],p[c("in2star","in3star","transitive","cycle")]);out[m$meanP]<-log(.2)}
 if(!is.null(m$triangle))out[m$meanT]<-log(.2)
 if(!is.null(m$meanDiamond))out[m$meanDiamond]<-log(.1)
 pmin(m$upper,pmax(m$lower,out))}
.cran_release_parameters<-function(m,window){if(!is.null(m$release_fits[[window]]))return(m$release_fits[[window]]$parameters);old<-.cran_old_fit(m,window);p<-setNames(unlist(old$parameters),unlist(old$parameter_names));n<-c("log_new_package_proposals","log_major_proposals","log_update_proposals",paste0(rep(c("major","update"),each=6),"_",rep(c("age","age_squared","in1","out1","in2","out2"),2)));setNames(if(is.null(old))c(log(c(.3,.2,.7)),rep(0,12))else p[n],n)}
#' @rdname cran_updates
#' @export
fit_cran_release<-function(data,window=c('train','full'),parameters=NULL,control=list()){
 m<-data;window<-match.arg(window)
 if(identical(m$mixture,"two_component"))stop("Fit complete-update mixture release and edge parameters jointly with fit_cran_updates()")
 take<-if(window=='train')m$year<=2023 else rep(TRUE,m$n);n<-sum(take)
 initial<-if(is.null(parameters)).cran_release_parameters(m,window)else parameters
 if(length(initial)!=15L||any(!is.finite(initial)))stop('Supply 15 finite release parameters')
 last<-initial;raw<-if(is.null(parameters)&&!is.null(m$release_values[[window]]))m$release_values[[window]]else cran_release_evaluate(m$release_data,initial,m$families)
 evaluate<-function(p){if(!identical(as.numeric(p),as.numeric(last))){last<<-p;raw<<-cran_release_evaluate(m$release_data,p,m$families)};list(value=-sum(raw$loglik[take])/n,gradient=-colSums(raw$gradient[take,,drop=FALSE])/n)}
 fit<-stats::optim(initial,function(p)evaluate(p)$value,function(p)evaluate(p)$gradient,method='L-BFGS-B',lower=c(rep(-15,3),rep(-Inf,12)),upper=c(rep(5,3),rep(Inf,12)),control=utils::modifyList(list(maxit=500L,pgtol=2e-6,factr=100),control))
 v<-evaluate(fit$par);names(fit$par)<-names(initial)
 list(parameters=fit$par,window=window,loglik=-n*v$value,gradient=-n*v$gradient,projected_gradient=max(abs(v$gradient)),convergence=fit$convergence,values=raw)
}
.cran_add_sum <- function(m,beta,age) {
 if(!identical(m$mixture,"two_component"))return(.cran_add_sum_uncached(m,beta,age))
 key <- c(as.numeric(beta),age)
 cached <- m$mixture_add_cache
 if(!is.null(cached)&&identical(cached$key,key))return(cached$value)
 value <- .cran_add_sum_uncached(m,beta,age)
 m$mixture_add_cache <- list(key=key,value=value)
 value
}
.cran_add_sum_uncached<-function(m,beta,age){a<-m$a;n<-m$n;ng<-m$ng
 if(m$pa){order<-order(a$update_batch);order<-order[a$update_batch[order]<n];ub<-a$update_batch[order];uo<-a$update_old_degree[order];un<-a$update_new_degree[order];bt<-a$update_birth[order]
  val<-function(d,b)exp(beta[1]*log1p(pmax(d,0))+age*b);wi<-val(a$initial_degree,a$initial_birth);wo<-val(uo,bt)*(uo>=0);wn<-val(un,bt)
  prefix<-function(fi,fo,fn)(sum(wi*fi)+cumsum(.cran_sum(ub,wn*fn-wo*fo,n)))[m$gb+1L]
  S<-prefix(1,1,1);G<-matrix(0,ng,m$na);G[,1]<-prefix(log1p(a$initial_degree),log1p(pmax(uo,0)),log1p(un));B<-prefix(a$initial_birth,bt,bt)
 }else{rp<-if(m$undir)c(beta[1],beta[2],beta[2],beta[3],beta[3],beta[4],0,beta[5],age)else if(m$exponential)c(0,beta[1],0,beta[2],0,beta[3:5],age)else c(beta[1:8],age)
  raw<-cran_candidate_sum(m$cache,rp,m$undir,exponential=m$exponential);S<-raw[,1];G<-matrix(0,ng,m$na);B<-raw[,10]
  if(m$exponential){G[,1:5]<-raw[,c(3,5,7:9)];scale<-raw[,11]}else if(m$undir)G[,1:5]<-cbind(raw[,2],raw[,3]+raw[,4],raw[,5]+raw[,6],raw[,7],raw[,9])else G[,1:8]<-raw[,2:9]
 }
 for(which in c("E","S")){x<-m[[paste0("X",which)]];b<-if(which=="E")a$existing_birth else a$group_source_birth;group<-if(which=="E")m$eg else 0:(ng-1);w<-.cran_weights(x,b,beta,age,m$pa||m$exponential);v<-exp(w$logw);S<-S-.cran_sum(group,v,ng);G<-G-.cran_sum(group,w$score*v,ng);B<-B-.cran_sum(group,v*b,ng)}
 b<-m$hist$features[,2];group<-m$hist$group;w0<-.cran_weights(m$XH0,b,beta,age,m$pa||m$exponential);w1<-.cran_weights(m$XH,b,beta,age,m$pa||m$exponential);v0<-exp(w0$logw);v1<-exp(w1$logw)
 S<-S+.cran_sum(group,v1-v0,ng);G<-G+.cran_sum(group,v1*w1$score-v0*w0$score,ng);B<-B+.cran_sum(group,(v1-v0)*b,ng)
 if(m$exponential&&any(!is.finite(S)|!is.finite(rowSums(G))|S<1e-9*scale)){
  bad<-which(!is.finite(S)|!is.finite(rowSums(G))|S<1e-9*scale)
  if(is.null(m$exp_cache))m$exp_cache<-cran_exponential_cache(m$cache)
  raw<-cran_exponential_sum(m$exp_cache,c(beta,age),bad-1L);S[bad]<-1;G[bad,]<-0;B[bad]<-0
  ans<-list(logsum=log(S),score=G/S,birth=B/S);ans$logsum[bad]<-raw[bad,1];ans$score[bad,]<-raw[bad,2:8,drop=FALSE];ans$birth[bad]<-raw[bad,9];return(ans)
 }
 if(any(S[m$exposure[,1]>0]<=0))stop("Nonpositive full-candidate weight sum")
 S[m$exposure[,1]==0]<-1;list(logsum=log(S),score=G/S,birth=B/S)
}
.cran_logexpm1<-function(l){h<-exp(l);ifelse(l< -30,l,h+log(-expm1(-h)))}
.cran_bernoulli_score<-function(l){h<-exp(l);ifelse(h<1e-4,1-h/2+h*h/12,h*exp(-h)/(-expm1(-h)))}
.cran_logbern<-function(l)ifelse(l< -30,l,log(-expm1(-exp(l))))
.cran_softplus<-function(x)pmax(x,0)+log1p(exp(-abs(x)))
.cran_triangle_correction<-function(m,p,la,rg,migration_log_odds,migration_gradient){
  x<-m$triangle; d<-x$diamond; shared<-x$consumer; n<-m$ng
  rate_col<-m$meanT[m$Tclass]; log_mean<-p[rate_col]
  mt<-exp(log_mean)*(x$mass[,1]>0)
  rate_d<-m$meanDiamond[m$Tclass];md<-if(is.null(d))numeric(n)else exp(p[rate_d])*(d$mass>0)
  rate_c<-m$meanConsumer[m$Tclass];mc<-if(is.null(shared))numeric(n)else exp(p[rate_c])*(shared$mass>0)
  mean<-mt+md+mc
  # Removing the requirement on an observed addition gives the inverse of its
  # Bernoulli probability, except when it can share a migration with a deletion.
  log_ra<-.cran_logbern(la); log_f<--log_ra
  df<--rg[seq_along(la),,drop=FALSE]*.cran_bernoulli_score(la)
  additions_only<-identical(m$law,"CS-pair-add")
  if(!additions_only){
    a<-m$ai; u<-migration_log_odds
    log_f[a]<-log_f[a]+.cran_softplus(u+log_ra[a])-.cran_softplus(u)
    df[a,]<-df[a,,drop=FALSE]*(1-stats::plogis(u+log_ra[a]))+
      migration_gradient*(stats::plogis(u+log_ra[a])-stats::plogis(u))
  }

  log_mix<--mean; gradient<-matrix(0,n,length(p)); pos_t<-cbind(seq_len(n),rate_col)
  gradient[pos_t]<--mt
  if(!is.null(d)){pos_d<-cbind(seq_len(n),rate_d);gradient[pos_d]<--md}
  if(!is.null(shared)){pos_c<-cbind(seq_len(n),rate_c);gradient[pos_c]<--mc}
  group<-c(x$group,d$group,shared$group);i<-group+1L
  if(length(i)){
    pair_copy<-!additions_only&&!identical(m$law,"CS-pair-migration-add")
    copy_effect<-if(pair_copy)p[m$copyT] else 0
    wc<-x$mass[,2]*exp(copy_effect);Z<-x$mass[,1]-x$mass[,2]+wc;it<-x$group+1L
    lt<-log_mean[it]+x$log_weight+copy_effect*x$copied-log(Z[it])
    id<-d$group+1L;ld<-if(is.null(d))numeric()else p[rate_d[id]]+d$log_weight-log(d$mass[id])
    ic<-shared$group+1L;lc<-if(is.null(shared))numeric()else p[rate_c[ic]]+shared$log_weight-log(shared$mass[ic])
    if(isTRUE(m$pair_proposals %in% c('poisson','capped2'))){
      # Multiple types can propose the same pair. The native evaluator sums
      # their Poisson rates, then enumerates the positive union polynomial.
      lp<-c(lt,ld,lc);dp<-matrix(0,length(lp),length(p))
      row<-seq_along(lp);column<-c(rate_col[it],rate_d[id],rate_c[ic])
      dp[cbind(row,column)]<-1
      if(pair_copy)dp[seq_along(it),m$copyT]<-x$copied-wc[it]/Z[it]
      position<-pos_t;exposure<-mt
      if(!is.null(d)){position<-rbind(position,pos_d);exposure<-c(exposure,md)}
      if(!is.null(shared)){position<-rbind(position,pos_c);exposure<-c(exposure,mc)}
      if(identical(m$pair_proposals,'capped2')){
        union<-cran_pair_capped2_evaluate(m$pair_union_plan,log_f,df,lp,dp,mean,-gradient)
        return(list(loglik=union$loglik,gradient=union$gradient,mean=mean,position=position,exposure=exposure))
      }
      union<-cran_pair_union_evaluate(m$pair_union_plan,log_f,df,lp,dp,m$pair_max_states)
      return(list(loglik=log_mix+union$log_partition,gradient=gradient+union$gradient,
                  mean=mean,position=position,exposure=exposure))
    }
    # Archived at-most-one paired addition, with one contribution per type.
    positive<-mean>0; common<-numeric(n);common[positive]<-log(-expm1(-mean[positive])/mean[positive])
    a<-c(x$a,d$a,shared$a);b<-c(x$b,d$b,shared$b);lw<-common[i]+c(lt,ld,lc)+log_f[a]+log_f[b]
    norm<-.cran_norm(group,lw,n); groups<-unique(i)
    odds<-norm$logsum+mean
    log_mix[groups]<-log_mix[groups]+.cran_softplus(odds[groups])
    any_prob<-numeric(n);any_prob[groups]<-stats::plogis(odds[groups])
    resp<-any_prob[i]*norm$probability
    count<-.cran_sum(c(a,b)-1L,c(resp,resp),length(la))
    gradient<-gradient+.cran_sum(m$ag,df*count,n)
    dcommon<-numeric(n);dcommon[positive]<-(.cran_bernoulli_score(log(mean[positive]))-1)/mean[positive]
    rt<-seq_along(it); rd<-length(it)+seq_along(id)
    gradient[pos_t]<-gradient[pos_t]+any_prob*mt*(1+dcommon)+.cran_sum(x$group,resp[rt],n)
    if(pair_copy)gradient[,m$copyT]<-.cran_sum(x$group,resp[rt]*(x$copied-wc[it]/Z[it]),n)
    if(!is.null(d))gradient[pos_d]<-gradient[pos_d]+any_prob*md*(1+dcommon)+.cran_sum(d$group,resp[rd],n)
    if(!is.null(shared)){
      rc<-length(it)+length(id)+seq_along(ic)
      gradient[pos_c]<-gradient[pos_c]+any_prob*mc*(1+dcommon)+.cran_sum(shared$group,resp[rc],n)
    }
  }
  position<-pos_t;exposure<-mt
  if(!is.null(d)){position<-rbind(position,pos_d);exposure<-c(exposure,md)}
  if(!is.null(shared)){position<-rbind(position,pos_c);exposure<-c(exposure,mc)}
  list(loglik=log_mix,gradient=gradient,mean=mean,position=position,exposure=exposure)
}
#' @rdname cran_updates
#' @export
loglik_cran_updates <- function(data, parameters, per_update=FALSE) {
 if(identical(data$mixture,"two_component"))
  return(.cran_mixture_loglik(data,parameters,per_update))
 if(inherits(parameters,"cran_update_fit"))parameters<-parameters$parameters
 .cran_single_loglik(data,parameters,per_update)
}
.cran_single_loglik<-function(data,parameters,per_update=FALSE){m<-data;p<-unname(.cran_parameters(parameters,m$names,"edge parameters"))
 if(identical(m$pair_proposals,'poisson'))cran_pair_union_check(m$pair_union_plan,m$pair_max_states)
 if(m$law %in% c("TERGM", "TERGM-logdegree")){
  logdegree <- m$law == "TERGM-logdegree"
  shifts <- NULL
  if(logdegree){
   add_covariates <- m$ZA[, -1, drop = FALSE]
   delete_covariate <- m$ZD[, 2]
   shifts <- cbind(drop(add_covariates %*% p[20:22]), delete_covariate * p[23])
  }
  raw <- cran_candidate_sum(m$cache, p[1:19], tergm = TRUE, sender_offsets = shifts,
                            motif_free = !is.null(m$shared_people))
  H <- raw[, 1]
  dh <- raw[, -1]
  if(logdegree){
   # Each sender covariate shifts every candidate logit for its operation.
   dh <- cbind(dh, add_covariates * dh[, 1], delete_covariate * dh[, 11])
  }
  if(!is.null(m$shared_people)){
   adjusted <- .cran_people_tergm_sum(m,H,dh,p)
   H <- adjusted$H
   dh <- adjusted$dh
  }
  ll <- drop(m$obs %*% p) - H
  ll[m$old] <- ll[m$old] - log(-expm1(-H[m$old]))
  w <- rep(1, m$ng)
  w[m$old] <- 1 / (-expm1(-H[m$old]))
  G <- m$obs - dh * w
 }else{
  a<-m$a;ng<-m$ng;dim<-length(p);A<-.cran_add_sum(m,p[m$ia],p[m$ageA]);wa<-.cran_weights(m$XA,a$added_birth,p[m$ia],p[m$ageA],m$pa||m$exponential);wd<-.cran_weights(m$XD,a$existing_birth,p[m$id],p[m$ageD],m$pa||m$exponential)
  add_index <- m$ia
  delete_index <- m$id
  if(!is.null(m$shared_people)){
    people <- m$shared_people
    add_people <- people$index[c(1L,if(m$law=="CS-1")3L)]
    A <- .cran_people_add_sum(m,A,p,add_people)
    wa$logw <- wa$logw+drop(people$observed_add%*%p[add_people])
    wa$score <- cbind(wa$score,people$observed_add)
    wd$logw <- wd$logw+p[people$index[2L]]*people$deletion
    wd$score <- cbind(wd$score,people$deletion)
    add_index <- c(add_index,add_people)
    delete_index <- c(delete_index,people$index[2L])
  }
  if(isTRUE(m$migration_degree)&&m$law=="CS-1"){
    A <- .cran_degree_add_sum(m,A,p[m$ia],p[m$ageA],p[m$degree_index])
    wa$logw <- wa$logw+p[m$degree_index]*m$degree_observed
  }
  dn<-.cran_norm(m$eg,wd$logw,ng);D<-.cran_sum(m$eg,wd$score*dn$probability,ng);Db<-.cran_sum(m$eg,a$existing_birth*dn$probability,ng)
  meanA<-exp(drop(m$ZA%*%p[m$rateA]));meanD<-exp(drop(m$ZD%*%p[m$rateD]))
  la<-log(meanA[m$ag+1L])+wa$logw-A$logsum[m$ag+1L];ld<-log(meanD[m$sg+1L])+wd$logw[m$sel]-dn$logsum[m$sg+1L];lr<-c(la,ld);nn<-length(la)
  rg<-matrix(0,length(lr),dim);rg[seq_len(nn),m$rateA]<-m$ZA[m$ag+1L,,drop=FALSE];rg[nn+seq_along(ld),m$rateD]<-m$ZD[m$sg+1L,,drop=FALSE]
  rg[seq_len(nn),add_index]<-wa$score-A$score[m$ag+1L,,drop=FALSE];rg[nn+seq_along(ld),delete_index]<-wd$score[m$sel,,drop=FALSE]-D[m$sg+1L,,drop=FALSE]
  rg[seq_len(nn),m$ageA]<-a$added_birth-A$birth[m$ag+1L];rg[nn+seq_along(ld),m$ageD]<-a$existing_birth[m$sel]-Db[m$sg+1L]
  if(isTRUE(m$migration_degree)&&m$law=="CS-1")
    rg[seq_len(nn),m$degree_index] <- m$degree_observed-A$degree[m$ag+1L]
  groups<-c(m$ag,m$sg);HA<-m$exposure[,1]*meanA;HD<-m$exposure[,2]*meanD;H<-HA+HD;ll<--H+.cran_sum(groups,.cran_logexpm1(lr),ng)
  mult<-exp(lr)+.cran_bernoulli_score(lr);G<-.cran_sum(groups,rg*mult,ng);dH<-matrix(0,ng,dim);dH[,m$rateA]<-m$ZA*HA;dH[,m$rateD]<-m$ZD*HD;G<-G-dH
  u<-du<-NULL
  if(!is.null(m$meanP)&&m$law!="CS-pair-add"){
   pp <- if(is.null(m$pair_preferences))c(m$ia[m$na-2],tail(m$id,1),m$ageA)else m$pair_preferences
   z<-drop(m$PX%*%p[m$ip])+p[pp[1]]*m$Ppop[,1]+p[pp[2]]*m$Ppop[,2]
   lw<-(if(m$exponential)z else stats::plogis(z,log.p=TRUE))+p[pp[3]]*m$Pb;pn<-.cran_norm(m$pg,lw,ng);lp<-p[m$meanP]+lw-pn$logsum[m$pg+1L]
   dl<-matrix(0,length(lp),dim);s<-if(m$exponential)1 else stats::plogis(-z);dl[,m$ip]<-m$PX*s;dl[,pp[1]]<-s*m$Ppop[,1];dl[,pp[2]]<-s*m$Ppop[,2];dl[,pp[3]]<-m$Pb
   dl<-dl-.cran_sum(m$pg,dl*pn$probability,ng)[m$pg+1L,,drop=FALSE];dl[,m$meanP]<-dl[,m$meanP]+1
   if(isTRUE(m$migration_degree)||!is.null(m$shared_people)){
    tilt_index <- if(isTRUE(m$migration_degree))m$degree_index else m$shared_people$index[3L]
    tilt <- if(isTRUE(m$migration_degree))m$degree_pair else m$shared_people$pair
    lp <- lp+p[tilt_index]*tilt
    dl[,tilt_index] <- tilt
    HP <- .cran_sum(m$pg,exp(lp),ng)
    pair_dH <- .cran_sum(m$pg,exp(lp)*dl,ng)
   } else {
    HP <- exp(p[m$meanP])*m$hasP
    pair_dH <- matrix(0,ng,dim);pair_dH[,m$meanP] <- HP
   }
   H<-H+HP;ll<-ll-HP;G<-G-pair_dH;dH<-dH+pair_dH
   u<-.cran_logexpm1(lp[m$oi])-.cran_logbern(lr[m$ai])-.cran_logbern(lr[m$di]);ll<-ll+.cran_sum(m$og,.cran_softplus(u),ng)
   du<-(exp(lp[m$oi])+.cran_bernoulli_score(lp[m$oi]))*dl[m$oi,,drop=FALSE]-.cran_bernoulli_score(lr[m$ai])*rg[m$ai,,drop=FALSE]-.cran_bernoulli_score(lr[m$di])*rg[m$di,,drop=FALSE]
   G<-G+.cran_sum(m$og,stats::plogis(u)*du,ng)
  }
  if(!is.null(m$triangle)){
   tri<-.cran_triangle_correction(m,p,la,rg,u,du)
   ll<-ll+tri$loglik;G<-G+tri$gradient;H<-H+tri$mean;dH[tri$position]<-tri$exposure
  }
  ll[m$old]<-ll[m$old]-log(-expm1(-H[m$old]));G[m$old,]<-G[m$old,,drop=FALSE]-dH[m$old,,drop=FALSE]/expm1(H[m$old])
 }
 ll<-.cran_sum(m$gb,ll,m$n);G<-.cran_sum(m$gb,G,m$n);if(per_update)list(loglik=ll,gradient=G)else list(loglik=sum(ll),gradient=colSums(G))
}
.cran_projected <- function(p,g,lo,hi){g[(p<=lo+1e-7 & g>0)|(p>=hi-1e-7 & g<0)]<-0;g}
#' @rdname cran_updates
#' @export
fit_cran_updates <- function(data,window=c("train","full"),starts=2L,control=list(),parameters=NULL,checkpoint=NULL){
 m<-data;window<-match.arg(window);if(length(starts)!=1L||!is.finite(starts)||starts<2L||starts!=floor(starts))stop("Use at least two integer starts")
 if(identical(m$mixture,"two_component"))
  return(.cran_fit_mixture(m,window,starts,control,parameters,checkpoint))
 if(identical(m$pair_proposals,'poisson'))cran_pair_union_check(m$pair_union_plan,m$pair_max_states)
 polish_tol<-if(is.null(control$polish_tol))2e-7 else control$polish_tol
 if(length(polish_tol)!=1L||!is.finite(polish_tol)||polish_tol<=0||polish_tol>2e-6)stop("Invalid polishing tolerance")
 control$polish_tol<-NULL
 take<-if(window=="train")m$year<=2023 else rep(TRUE,m$n);n<-sum(take)
 if(inherits(parameters,"cran_update_fit"))parameters<-parameters$parameters
 supplied<-if(is.list(parameters))lapply(parameters,.cran_parameters,expected=m$names,label="starting parameters")else NULL
 initial<-if(is.null(parameters)).cran_old_start(m,window)else if(!is.null(supplied))supplied[[1]]else .cran_parameters(parameters,m$names,"starting parameters")
 if(length(initial)!=length(m$names))stop("Wrong start length")
 # Apply fixed-term restrictions to supplied and archived warm starts alike.
 initial <- pmin(m$upper, pmax(m$lower, initial))
 if(!is.null(supplied)) supplied <- lapply(supplied, function(p) {
   if(length(p)!=length(m$names))stop("Wrong start length")
   pmin(m$upper, pmax(m$lower, p))
 })
 last<-NULL;value<-NULL;calls<-0L
 evaluate<-function(p){if(!identical(p,last)){v<-loglik_cran_updates(m,p,TRUE);last<<-p;value<<-list(value=-sum(v$loglik[take])/n,gradient=-colSums(v$gradient[take,,drop=FALSE])/n);calls<<-calls+1L;if(calls%%50L==0L)message(m$law," ",window," evaluation ",calls,": ",format(-value$value*n,digits=12))};value}
 raw<-loglik_cran_updates(m,initial,TRUE)$gradient[take,,drop=FALSE]
 scales<-pmin(100,pmax(.001,1/sqrt(colMeans(raw^2)+1e-10)))
 objective<-function(p)evaluate(p)$value;gradient<-function(p)evaluate(p)$gradient
 polish<-function(p){
  for(iter in seq_len(12L)){
   g<-gradient(p);pg<-.cran_projected(p,g,m$lower,m$upper);if(max(abs(pg))<polish_tol)break
   free<-which(!((p<=m$lower+1e-7 & g>0)|(p>=m$upper-1e-7 & g<0)))
   H<-matrix(0,length(free),length(free));for(k in seq_along(free)){j<-free[k];h<-1e-7*max(1,abs(p[j]));d<-rep(0,length(p));d[j]<-h;H[,k]<-(gradient(p+d)-gradient(p-d))[free]/(2*h)}
   H<-(H+t(H))/2;s<-1/sqrt(pmax(abs(diag(H)),1e-8));hs<-H*outer(s,s);ev<-eigen(hs,symmetric=TRUE);dir<-rep(0,length(p));dir[free]<--s*drop(ev$vectors%*%(drop(crossprod(ev$vectors,g[free]*s))/pmax(ev$values,1e-4)))
   dir[(p<=m$lower+1e-7 & dir<0)|(p>=m$upper-1e-7 & dir>0)]<-0
   if(sum(dir*g)>=0)dir<--pg*scales^2
   step<-min(1,if(any(dir>0))min((m$upper[dir>0]-p[dir>0])/dir[dir>0])else 1,if(any(dir<0))min((m$lower[dir<0]-p[dir<0])/dir[dir<0])else 1)
   old<-objective(p);while(step>1e-10){q<-pmin(m$upper,pmax(m$lower,p+step*dir));if(objective(q)<=old+1e-4*step*sum(g*dir))break;step<-step/2}
   if(step<=1e-10)break;p<-q
  };p
 }
 trials<-vector("list",starts)
 for(i in seq_len(starts)){
  start<-if(!is.null(supplied)&&i<=length(supplied))supplied[[i]]else initial
  if(i>1L&&is.null(supplied)){start<-start+.1*scales*sin(seq_along(start)*i);start<-pmin(m$upper,pmax(m$lower,start))}
  result<-stats::optim(start,objective,gradient,method="L-BFGS-B",lower=m$lower,upper=m$upper,control=utils::modifyList(list(maxit=700L,factr=1e3,pgtol=1e-8,lmm=24L,parscale=scales),control))
  if(!is.null(checkpoint))saveRDS(result$par,checkpoint)
  result$par<-polish(result$par);result$value<-objective(result$par);result$projected_gradient<-max(abs(.cran_projected(result$par,gradient(result$par),m$lower,m$upper)));trials[[i]]<-result
  if(!is.null(checkpoint))saveRDS(result$par,checkpoint)
  message(m$law," ",window," start ",i,": ",format(-result$value*n,digits=13),"; projected gradient ",format(result$projected_gradient,digits=3))
 }
 best<-trials[[which.min(vapply(trials,`[[`,0,"value"))]]
 # A Newton step can leave a shallow stationary point and find a better basin.
 # Re-optimise that improved point before treating it as a completed fit.
 for(cycle in seq_len(6L)){
  if(best$projected_gradient<2e-6)break
  raw<-loglik_cran_updates(m,best$par,TRUE)$gradient[take,,drop=FALSE]
  scales<-pmin(100,pmax(1e-5,1/sqrt(colMeans(raw^2)+1e-12)))
  port<-stats::nlminb(best$par,objective,gradient,lower=m$lower,upper=m$upper,scale=1/scales,control=list(eval.max=1500L,iter.max=700L,rel.tol=1e-13,x.tol=1e-9))
  candidate<-polish(port$par);candidate_value<-objective(candidate)
  if(candidate_value<=best$value+1e-12){best$par<-candidate;best$value<-candidate_value}
  best$projected_gradient<-max(abs(.cran_projected(best$par,gradient(best$par),m$lower,m$upper)))
  if(!is.null(checkpoint))saveRDS(best$par,checkpoint)
  message(m$law," ",window," refinement ",cycle,": ",format(-best$value*n,digits=13),"; projected gradient ",format(best$projected_gradient,digits=3))
 }
 if(best$projected_gradient>2e-5)stop("Projected gradient remains above tolerance; checkpoint retained, fit not finalised")
 v<-loglik_cran_updates(m,best$par,TRUE);names(best$par)<-m$names
 rp<-.cran_release_parameters(m,window);r<-if(!is.null(m$release_values[[window]]))m$release_values[[window]] else cran_release_evaluate(m$release_data,rp,m$families)
 bound<-abs(best$par-m$lower)<1e-5 | abs(best$par-m$upper)<1e-5
 G<-v$gradient[take,,drop=FALSE];opg<-crossprod(G);se<-rep(NA_real_,length(best$par));free<-which(!bound)
 if(length(free)){of<-opg[free,free,drop=FALSE];os<-sqrt(pmax(diag(of),1e-300));correlation<-of/outer(os,os);if(rcond(correlation)>1e-10)se[free]<-sqrt(diag(solve(correlation)))/os}
 prediction<-lapply(c(2024L,2025L),function(y){k<-m$year==y;list(year=y,updates=sum(k),edge_loglik=sum(v$loglik[k]),release_loglik=sum(r$loglik[k]),update_loglik=sum(v$loglik[k]+r$loglik[k]))});names(prediction)<-c("2024","2025")
 structure(list(law=m$law,window=window,parameters=best$par,release_parameters=rp,parameter_count=sum(m$lower<m$upper)+15L,edge_loglik=sum(v$loglik[take]),update_loglik=sum(v$loglik[take]+r$loglik[take]),predictive=prediction,standard_error=se,bound_parameters=m$names[bound & m$lower<m$upper],fixed_parameters=m$names[m$lower==m$upper],projected_gradient=best$projected_gradient,multistart_loglik=-n*vapply(trials,`[[`,0,"value"),convergence=vapply(trials,`[[`,0L,"convergence"),contributions=list(year=m$year,edge=v$loglik,release=r$loglik),weights=if(m$law %in% c("TERGM", "TERGM-logdegree"))"logit"else if(m$pa)"exponential"else m$weights,proposal_means=m$proposal_means,term_specification=m$term_specification,migration_degree=isTRUE(m$migration_degree),pair_proposals=m$pair_proposals,pair_max_states=m$pair_max_states,specification=paste("separate addition/deletion coefficients and age decay; no penalty;", m$term_specification, "term specification")),class="cran_update_fit")
}

# Exponential tilt of the migration-candidate subset. The complement remains
# unchanged; derivatives include the complete candidate normalizer.
.cran_degree_add_sum <- function(m,A,beta,age,eta){
  keep <- which(m$degree_history>0)
  group <- m$hist$group[keep];i <- group+1L
  x <- m$XH[keep,,drop=FALSE];b <- m$hist$features[keep,2]
  degree <- m$degree_history[keep]
  w <- exp(drop(x%*%beta)+age*b-A$logsum[i])
  mass <- .cran_sum(group,w,m$ng)
  if(any(mass>1+1e-8))stop("Migration subset exceeds complete normalizer")
  tilted <- w*exp(eta*degree)
  denominator <- pmax(0,1-mass)+.cran_sum(group,tilted,m$ng)
  score <- (A$score+.cran_sum(group,(tilted-w)*x,m$ng))/denominator
  birth <- (A$birth+.cran_sum(group,(tilted-w)*b,m$ng))/denominator
  list(logsum=A$logsum+log(denominator),score=score,birth=birth,
       degree=.cran_sum(group,tilted*degree,m$ng)/denominator)
}
