#' Detection probability of each subject's observed screening outcome
#' @noRd
p_v_obs_2 = function(v_obs, kappa, r){
  m = sapply(v_obs, length)-1
  cens =  sapply(v_obs, function(x) is.infinite(x[length(x)]))
  p_vec = geom(m+r, kappa)
  p_vec_inf = geom_inf(m+r, kappa)
  ifelse(cens, p_vec_inf, p_vec)
}
