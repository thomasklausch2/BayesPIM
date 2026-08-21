#' Draw from a truncated parametric time-to-event distribution by inversion
#' @noRd
r_trdist =  function(par, a = 0, b=Inf, dist , tol = 1e-8){ 
    n = length(a)
    cdf_a = pdist(a, par, dist)
    cdf_b = pdist(b, par, dist)
    dif = cdf_b-cdf_a
    dif = dif < tol
    u   = ifelse(dif, a, qdist( cdf_a + runif(n) * (cdf_b-cdf_a), par, dist ))
    return(u)
  }
