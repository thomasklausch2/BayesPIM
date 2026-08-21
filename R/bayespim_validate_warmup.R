#' Validate a warmup value for BayesPIM methods
#' @noRd
bayespim_validate_warmup <- function(warmup, n_iter) {
  valid <- is.numeric(warmup) &&
    length(warmup) == 1L &&
    !is.na(warmup) &&
    is.finite(warmup) &&
    warmup >= 0 &&
    abs(warmup - round(warmup)) < sqrt(.Machine$double.eps)

  if (!valid) {
    stop("`warmup` must be a single non-negative whole number.", call. = FALSE)
  }

  warmup <- as.integer(round(warmup))
  if (warmup >= n_iter) {
    stop("`warmup` must leave at least one draw in each chain.", call. = FALSE)
  }

  warmup
}

