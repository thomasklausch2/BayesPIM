#' fc_beta
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
fc_beta <- function(X, y, sig_inv_Xt, sig_inv){
  p = ncol(X)
  beta_hat = sig_inv_Xt %*% y 
  mvrnorm(1, beta_hat, sig_inv)
}