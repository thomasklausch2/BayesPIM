#' Map lognormal AFT coefficients to meanlog and sdlog parameters
#' @noRd
trans_par_ind_norm = function(p, v, x1){
  mu = x1 %*% as.matrix(p)
  sd = exp(v)
  cbind(mu,sd)
}
