#' Geometric probability of first detection at a given screening index
#' @noRd
geom = function(j, kappa) kappa * (1 - kappa)^(j - 1)
