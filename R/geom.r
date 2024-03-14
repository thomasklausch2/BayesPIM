#' geom
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
geom = function(j, kappa) kappa * (1 - kappa)^(j - 1)
