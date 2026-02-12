#' @keywords internal
"_PACKAGE"

## -- Suppress R CMD check notes for non-standard evaluation symbols ------
## ggplot2 aes(), data.table column references, closure variables, etc.
if (getRversion() >= "2.15.1") {
  globalVariables(c(
    # ggplot2 aes() aesthetics (used in gof.R, utils.R)
    "time", "intensity", "density", "value", "x", "type",
    "proportion", "distance", "waiting_time", "statistic",
    "ymin", "ymax",
    # data.table NSE (used in hawkesNet.R)
    ".SD", ".N", ".I",
    # closure variables injected by cond_intensity / cond_intensity_inhom
    "diffs_local", "ldf_local",
    # other closure / loop variables
    "i", "j"
  ))
}

## -- network -------------------------------------------------------------
#' @importFrom network network network.initialize network.size network.copy
#'   is.directed
#'   add.edge add.edges add.vertices delete.edges delete.vertices
#'   set.vertex.attribute get.vertex.attribute delete.vertex.attribute
#'   list.vertex.attributes
#'   set.edge.attribute get.edge.attribute list.edge.attributes
#'   set.network.attribute
#'   as.edgelist get.edgeIDs get.dyads.eids
#'   as.matrix.network.edgelist
#'   "%v%" "%n%" "%e%"
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
#' @importFrom stats approx approxfun as.formula density dnorm dpois ks.test
#'   optim pexp plogis rpois runif rexp sd setNames
NULL

## -- utils ---------------------------------------------------------------
#' @importFrom utils head tail relist
NULL

## -- parallel ------------------------------------------------------------
#' @importFrom parallel makeCluster stopCluster clusterEvalQ parLapply
#'   mclapply detectCores
NULL

## -- graphics ------------------------------------------------------------
#' @importFrom graphics par plot points abline
NULL
