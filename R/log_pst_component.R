#' Evaluate a one-component change in an AFT log-posterior
#' @noRd
log_pst_component <- function(v, j, eta, log_ll, ...) {
  if (j > length(eta)) stop("`j` exceeds the length of `eta`.", call. = FALSE)
  eta[j] <- v
  log_pst(eta = eta, log_ll = log_ll, ...)
}
