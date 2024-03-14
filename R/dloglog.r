#' dloglog
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
dloglog = function(x,lambda, gamma){ dllogis(x, shape = gamma, scale= lambda ) }
