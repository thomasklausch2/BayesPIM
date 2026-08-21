#' Draw regression coefficients from their multivariate normal full conditional
#' @noRd
fc_beta <- function(x, y, sig_inv_xt, sig_inv){
  beta_hat = sig_inv_xt %*% y 
  mvrnorm(1, beta_hat, sig_inv)
}
