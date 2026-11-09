#' Independent Poisson proposals of unordered edit sets
#'
#' Implements Constructions 1 and 2: each candidate group has an independent
#' Poisson count with the supplied mean, and the update is the union of groups
#' with positive counts. Edit identifiers are arbitrary character strings;
#' callers define the admissible edits against one fixed pre-update graph.
#' For a size-conditional channel use rates
#' `omega * pi[size] * weight / sum(weights_of_that_size)`. Several channels
#' can be combined by concatenating their groups and rates.
#'
#' Probabilities marginalise all proposal combinations, including overlapping
#' groups. Evaluation factors over connected components of overlapping groups.
#' An observed component with more than `max_component` edits raises an error;
#' no support restriction or independent-edit approximation is substituted.
#' Singleton components are evaluated analytically.
#' @param groups List of nonempty character vectors, one per candidate group.
#' @param rates Nonnegative finite mean proposal counts, one per group.
#' @param nonempty Condition on at least one proposed edit.
#' @param law Object returned by `poisson_proposals()`.
#' @param update Character vector of distinct observed edit identifiers.
#' @param max_component Largest connected observed component to enumerate,
#'   an integer from 1 to 20. The cost can be exponential in this value.
#' @return `poisson_proposals()` returns a law; `sample_update()` a character
#'   vector of selected edits; `logprob_update()` a scalar log probability.
#' @examples
#' law <- poisson_proposals(list("delete_old", "add_new",
#'                              c("delete_old", "add_new")),
#'                          c(0.2, 0.05, 0.15))
#' exp(logprob_update(law, c("delete_old", "add_new")))
#' set.seed(1)
#' sample_update(law)
#' @export
poisson_proposals <- function(groups, rates, nonempty = FALSE) {
  if (!is.list(groups) || !is.numeric(rates) || length(groups) != length(rates) ||
      any(!is.finite(rates)) || any(rates < 0) || !is.finite(sum(rates)))
    stop("groups and finite nonnegative rates must have matching lengths")
  if (!is.logical(nonempty) || length(nonempty) != 1L || is.na(nonempty))
    stop("nonempty must be TRUE or FALSE")
  if (any(!vapply(groups, function(g) is.character(g) && length(g) > 0L &&
                 !anyNA(g) && all(nzchar(g)) && !anyDuplicated(g), logical(1))))
    stop("each group must contain distinct, nonempty character edit identifiers")
  if (nonempty && sum(rates) == 0) stop("nonempty conditioning needs a positive rate")
  structure(list(groups = groups, rates = rates, nonempty = nonempty),
            class = "poisson_proposals")
}

# Draw independent Bernoulli(1-exp(-rate)) variables conditional on a success.
# Conditioning on the first success avoids rejection when total rate is tiny.
.proposal_indicators <- function(rates, nonempty = FALSE) {
  if (!nonempty) return(stats::runif(length(rates)) < -expm1(-rates))
  total <- sum(rates)
  first_weights <- exp(-c(0, head(cumsum(rates), -1L))) * (-expm1(-rates))
  first <- sample.int(length(rates), 1L, prob = first_weights / (-expm1(-total)))
  selected <- rep(FALSE, length(rates)); selected[first] <- TRUE
  later <- seq_along(rates) > first
  selected[later] <- stats::runif(sum(later)) < -expm1(-rates[later])
  selected
}

#' @rdname poisson_proposals
#' @export
sample_update <- function(law) {
  if (!inherits(law, "poisson_proposals")) stop("law must come from poisson_proposals()")
  selected <- .proposal_indicators(law$rates, law$nonempty)
  sort(unique(as.character(unlist(law$groups[selected], use.names = FALSE))))
}

#' @rdname poisson_proposals
#' @export
logprob_update <- function(law, update, max_component = 20L) {
  if (!inherits(law, "poisson_proposals")) stop("law must come from poisson_proposals()")
  if (!is.character(update) || anyNA(update) || anyDuplicated(update))
    stop("update must contain distinct character edit identifiers")
  if (length(max_component) != 1L || !is.finite(max_component) ||
      max_component != as.integer(max_component) || max_component < 1L || max_component > 20L)
    stop("max_component must be an integer from 1 to 20")
  total <- sum(law$rates)
  if (!length(update)) return(if (law$nonempty) -Inf else -total)
  eligible <- vapply(law$groups, function(g) all(g %in% update), logical(1)) & law$rates > 0
  groups <- law$groups[eligible]; rates <- law$rates[eligible]
  if (!all(update %in% unlist(groups, use.names = FALSE))) return(-Inf)
  logp <- -sum(law$rates[!eligible])
  remaining <- seq_along(groups)
  while (length(remaining)) {
    component <- remaining[1L]; edits <- groups[[component]]
    repeat {
      joined <- remaining[vapply(groups[remaining], function(g) any(g %in% edits), logical(1))]
      expanded <- unique(unlist(groups[joined], use.names = FALSE))
      if (setequal(expanded, edits)) { component <- joined; break }
      edits <- expanded
    }
    remaining <- setdiff(remaining, component)
    if (length(edits) == 1L) {
      logp <- logp + log(-expm1(-sum(rates[component])))
      next
    }
    if (length(edits) > max_component)
      stop("observed overlap component exceeds max_component; exact union evaluation refused")
    # Log-space dynamic programming retains probabilities of very rare unions.
    states <- 0L; log_probabilities <- 0
    for (j in component) {
      mask <- sum(2L^(match(groups[[j]], edits) - 1L))
      target <- bitwOr(states, as.integer(mask))
      allstates <- c(states, target)
      values <- c(log_probabilities - rates[j], log_probabilities + log(-expm1(-rates[j])))
      combined <- split(values, allstates)
      states <- as.integer(names(combined))
      log_probabilities <- vapply(combined, function(x) {
        m <- max(x); if (is.infinite(m)) m else m + log(sum(exp(x - m)))
      }, numeric(1))
    }
    logp <- logp + log_probabilities[match(2^length(edits) - 1L, states)]
  }
  if (law$nonempty) logp <- logp - log(-expm1(-total))
  as.numeric(logp)
}
