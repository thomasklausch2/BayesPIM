#' Update all free AFT parameters with univariate slice sampling
#' @noRd
slice_passthrough <- function(eta, w = 1, log_ll, fix_sigma, fix_q, has_q, ...) {
  sigma_pos <- length(eta) - as.integer(has_q)
  
  for (i in seq_len(sigma_pos - 1L)) {
    eta[i] <- sample_slice_component(
      i, eta = eta, w = w, log_ll = log_ll, ...
    )
  }

  if (!fix_sigma) {
    eta[sigma_pos] <- sample_slice_component(
      sigma_pos, eta = eta, w = w, log_ll = log_ll, ...
    )
  }
  
  if (!fix_q && has_q) {
    eta[sigma_pos+1L] <- sample_slice_component(
      sigma_pos+1L, eta = eta, w = w, log_ll = log_ll, ...
    )
  }

  eta
}
