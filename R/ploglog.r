#' ploglog
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
ploglog = function(x,lambda, gamma){ pllogis(x, shape = gamma, scale= lambda ) }
