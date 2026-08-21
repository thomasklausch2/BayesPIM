#' Geometric probability of no detection through a given screening index
#' @noRd
geom_inf = function(j, kappa) (1 - kappa)^(j - 1)
