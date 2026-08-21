#' Density of the log-logistic distribution
#' @noRd
dloglog = function(x,lambda, gamma){ dllogis(x, shape = gamma, scale= lambda ) }
