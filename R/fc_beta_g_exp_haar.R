#' Draw the parameter-expansion variance for the probit prevalence model
#' @noRd
fc_beta_g_exp_haar <- function(y, x, sig_inv_xt){
  n       = length(y)
  beta_hat = sig_inv_xt %*% y
  rss    = sum( (y - x %*% beta_hat )^2)
  alpha_sq = rss/rchisq(1,df=n)
  alpha_sq
}
