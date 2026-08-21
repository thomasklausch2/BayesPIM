#' Draw test sensitivity from its beta full conditional without a prevalence component
#' @noRd
pst_kappa_no_prev = function( v_obs, j_, a = 1, b = 1){
  n_ = sum( sapply( v_obs, function(x) is.finite( x[length(x)] )) )
  m  = sapply(v_obs, length)
  alpha = n_ + a
  beta  = sum(m - j_ - 1) + b
  rbeta(1, shape1= alpha, shape2= beta)
}
