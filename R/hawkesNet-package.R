#' Likelihood and simulation for self-exciting network updates
#'
#' Models complete simultaneous network updates as marks of a point process.
#' The update law specifies which changes occur, and the ground intensity
#' specifies their timing. A normalized update law makes the compensator an
#' integral of the ground intensity alone.
#'
#' @section Interfaces:
#' \itemize{
#'   \item \code{\link{network_update_model}} and
#'   \code{\link{simulate_network_updates}} provide the small-network
#'   birth/turnover reference model with separate initial- and later-topology
#'   coefficients. \code{\link{fit_network_updates}} fits its likelihood.
#'   \item \code{\link{poisson_proposals}} defines an exact union law for
#'   arbitrary caller-supplied single-edit or grouped Poisson proposals.
#'   \item \code{\link{sim_hawkesNet}} and \code{\link{fit_hawkesNet}} provide
#'   addition-only examples using \code{network} histories. Select a mark law
#'   explicitly; for CS-1 use \code{PMF_mark_CS} with
#'   \code{cs_mode = "independent"}.
#'   \item \code{\link{cran_update_data}} and \code{\link{fit_cran_updates}}
#'   fit complete directed release/dependency updates, including births,
#'   additions, deletions, paired migrations and shared-person covariates.
#'   \item \code{\link{cran_timing_data}} and \code{\link{fit_cran_timing}}
#'   provide the calendar, release-type and topology-dependent timing models.
#'   \item \code{\link{simulate_cran_updates}} and
#'   \code{\link{simulate_cran_joint}} generate forward checks at fixed event
#'   times or with generated event times, respectively.
#' }
#'
#' @section Examples and data:
#' Short runnable scripts are installed in \code{system.file("examples",
#' package = "hawkesNet")}. The CRAN workflow guide is installed in
#' \code{system.file("cran", "README.md", package = "hawkesNet")}; its archive
#' inputs are external to the package. Historical study scripts and the dyadic
#' relational interfaces use different model definitions from the manuscript's
#' primary CRAN laws.
#'
#' @seealso \code{\link{hawkes_timing}}, \code{\link{loglik_hawkesNet}},
#'   \code{\link{loglik_cran_updates}}, \code{\link{loglik_cran_timing}}
#' @keywords package
"_PACKAGE"

## -- Suppress R CMD check notes for non-standard evaluation symbols ------
## ggplot2 aes(), data.table column references, closure variables, etc.
if (getRversion() >= "2.15.1") {
  globalVariables(c(
    # ggplot2 aes() aesthetics (used in gof.R, utils.R)
    "time", "intensity", "density", "value", "x", "type",
    "proportion", "distance", "waiting_time", "statistic",
    "ymin", "ymax", "log_waiting_time", "bin", "count", "level_num",
    # data.table NSE (used in hawkesNet.R)
    ".SD", ".N", ".I",
    # closure variables injected by cond_intensity / cond_intensity_inhom
    "diffs_local", "ldf_local",
    # other closure / loop variables
    "i", "j"
  ))
}

## -- network -------------------------------------------------------------
#' @importFrom network network network.initialize network.size network.edgecount network.copy is.network is.directed add.edge add.edges add.vertices delete.edges delete.vertices set.vertex.attribute get.vertex.attribute delete.vertex.attribute list.vertex.attributes set.edge.attribute get.edge.attribute list.edge.attributes set.network.attribute as.edgelist get.edgeIDs get.dyads.eids as.matrix.network.edgelist "%v%" "%n%" "%e%"
NULL

## -- sna -----------------------------------------------------------------
#' @importFrom sna degree geodist
NULL

## -- ernm ----------------------------------------------------------------
#' @importFrom ernm createCppModel as.BinaryNet calculateStatistics
NULL

## -- data.table ----------------------------------------------------------
#' @importFrom data.table data.table setkey rbindlist
NULL

## -- hash ----------------------------------------------------------------
#' @importFrom hash hash keys
NULL

## -- stats ---------------------------------------------------------------
#' @importFrom stats approx approxfun as.formula ave density dnorm dpois ks.test median optim pexp plogis rpois runif rexp sd setNames
NULL

## -- utils ---------------------------------------------------------------
#' @importFrom utils head tail relist
NULL

## -- parallel ------------------------------------------------------------
#' @importFrom parallel makeCluster stopCluster clusterEvalQ parLapply mclapply detectCores
NULL

## -- graphics ------------------------------------------------------------
#' @importFrom graphics par plot points abline
NULL
