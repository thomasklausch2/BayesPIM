#' Validate posterior probabilities for BayesPIM summaries
#' @noRd
bayespim_validate_probs <- function(probs) {
  if (
    !is.numeric(probs) ||
    length(probs) < 1L ||
    anyNA(probs) ||
    any(!is.finite(probs)) ||
    any(probs < 0 | probs > 1)
  ) {
    stop("`probs` must contain finite probabilities between 0 and 1.", call. = FALSE)
  }
  probs
}

