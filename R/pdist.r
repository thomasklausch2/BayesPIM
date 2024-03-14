#' pdist
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
pdist <- function(q, par, dist = "exp") {
  if (dist == "exp") {
    return(pexp(q, rate = par[, 1]))
  }
  if (dist == "weibull2") {
    return(pweibull2(q, beta = par[, 2], theta = par[, 1]))
  }
  if (dist == "gamma") {
    return(pgamma(q, shape = par[, 1], rate = par[, 2]))
  }
  if (dist == "weibull") {
    return(pweibull(q, shape = par[, 2], scale = par[, 1]))
  }
  if (dist == "loglog") {
    return(ploglog(q, lambda = par[, 1], gamma = par[, 2]))
  }
  if (dist == "lognormal") {
    return(plnorm(q, meanlog = par[, 1], sdlog = par[, 2]))
  }
  if (dist == "gengamma") {
    return(pggamma(q, a = 1/par[, 1], b = par[, 2], k = par[, 3]))
  }
}