#' Log-likelihood for accelerated failure time models
#' @noRd
ll_aft <- function(eta, y, x, dist) {
  # Guard before `eta` is sliced: an unsupported family would otherwise fail on
  # a length mismatch rather than on its own name. The generalized gamma is
  # deliberately absent; it is fitted through ll_aft_ic() under the collapsed
  # sampler.
  if (!dist %in% c("weibull", "loglog", "lognormal", "gamma")) {
    stop("Unsupported distribution: ", dist)
  }

  if (dist == "gamma") {
    if (any(!is.finite(y)) || any(y <= 0)) {
      return(-Inf)
    }

    gamma_par <- trans_par_gamma(x1 = x, par = eta)
    if (any(!is.finite(gamma_par)) || any(gamma_par <= 0)) {
      return(-Inf)
    }

    LL <- sum(stats::dgamma(
      y,
      shape = gamma_par[, "shape"],
      rate = gamma_par[, "rate"],
      log = TRUE
    ))

    if (!is.finite(LL)) {
      return(-Inf)
    }

    return(LL)
  }

  p = length(eta)
  beta=eta[1:(p-1)]
  sigma=exp(eta[p])
  log_sigma = eta[p]
  n = length(y)
  V = (log(y)- x %*% as.matrix(beta))/sigma

  if (any(!is.finite(V))) {
    return(-Inf)
  }

  if (any(!is.finite(sigma)) || any(sigma <= 0)) {
    return(-Inf)
  }

  if(dist=='weibull') LL = -n*log_sigma + sum( V - exp(V) )

  if(dist=='loglog'){
    logexpV = log(  1 + exp(-V) )
    logexpV[is.infinite(logexpV)]= -V[is.infinite(logexpV)]
    LL = -n*log_sigma - sum( V ) - 2 * sum( logexpV )
  }

  if(dist=='lognormal') LL = -n*log_sigma - 1/2 * sum( V^2 )
  
  LL
}
