#' Log of a difference between two survival probabilities
#' @noRd
logsdiff <- function(logSL, logSR) {
  if (length(logSL) != length(logSR)) {
    stop("`logSL` and `logSR` must have equal lengths.", call. = FALSE)
  }

  out <- rep(NA_real_, length(logSL))

  both_zero <- is.infinite(logSL) & logSL < 0 &
    is.infinite(logSR) & logSR < 0
  out[both_zero] <- -Inf

  valid <- !both_zero &
    is.finite(logSL) &
    !is.na(logSR) &
    logSR <= logSL

  delta <- logSR[valid] - logSL[valid]
  correction <- numeric(length(delta))

  use_log1p <- delta < -log(2)
  correction[use_log1p] <-
    log1p(-exp(delta[use_log1p]))
  correction[!use_log1p] <-
    log(-expm1(delta[!use_log1p]))

  out[valid] <- logSL[valid] + correction
  out
}
