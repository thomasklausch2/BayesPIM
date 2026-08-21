#' Detection probabilities across each subject's screening intervals
#' @noRd
p_v_obs = function(v_obs, kappa){
  m = sapply(v_obs, length)-1
  cens =  sapply(v_obs, function(x) is.infinite(x[length(x)]))
  m_max = max(m)
  p_vec = geom(1:m_max, kappa)
  p_vec_inf = geom_inf(1:m_max, kappa)
  P = sapply( 1: length(m), function(x) if(!cens[x]) rev(p_vec[1:m[x]])
              else         rev(p_vec_inf[1:m[x]])
  )
  P
}
