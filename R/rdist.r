#' rdist
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
rdist <- function(n, par, dist = "exp") {
  if (dist == "exp") {
    return(rexp(n, rate = par[, 1]))
  }
  if (dist == "weibull2") {
    return(rweibull2(n, beta = par[, 2], theta = par[, 1]))
  }
  if (dist == "gamma") {
    return(rgamma(n, shape = par[, 1], rate = par[, 2]))
  }
  if (dist == "weibull") {
    return(rweibull(n, shape = par[, 2], scale = par[, 1]))
  }
  if (dist == "loglog") {
    return(rloglog(n, lambda = par[, 1], gamma = par[, 2]))
  }
  if (dist == "lognormal") {
    return(rlnorm(n, meanlog = par[, 1], sdlog = par[, 2]))
  }
  if (dist == "gengamma") {
    return(rggamma(n, a = 1/par[, 1], b = par[, 2], k = par[, 3]))
  }
}