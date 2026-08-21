#' Map gamma AFT coefficients to shape and rate parameters
#' @noRd
trans_par_gamma <- function(x1, par) {
  p <- length(par)
  eta <- drop(x1 %*% par[seq_len(p - 1L)])
  sigma <- exp(par[p])
  shape <- sigma^(-2)
  rate <- shape * exp(-eta)

  cbind(
    shape = rep(shape, length(eta)),
    rate = rate
  )
}
