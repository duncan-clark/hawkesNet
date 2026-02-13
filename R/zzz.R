# Called when the package is loaded; runs before the namespace is locked.
# Reserves .ernm_model_cache so mark_PMF.R can use it without assign() failing on locked namespaces (e.g. some clusters).
.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  # Reserve the name before the environment is locked
  cache <- new.env()
  cache[[".cache_pid"]] <- Sys.getpid()
  assign(".ernm_model_cache", cache, envir = ns)
}
