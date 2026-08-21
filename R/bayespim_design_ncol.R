#' Count columns in an optional BayesPIM design matrix
#' @noRd
bayespim_design_ncol <- function(x) {
  if (is.null(x)) {
    return(0L)
  }
  ncol(as.matrix(x))
}

