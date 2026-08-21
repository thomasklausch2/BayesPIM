#' Convert a correlation matrix and standard deviations to a covariance matrix
#' @noRd
cor2cov <- function(R, S) {
  R <- as.matrix(R)
  if (length(S) != ncol(R)) {
    stop(
      sprintf(
        "`S` must supply one standard deviation per variable: %d supplied, %d expected.",
        length(S), ncol(R)
      ),
      call. = FALSE
    )
  }
  # Elementwise scaling by outer(S, S) equals diag(S) %*% R %*% diag(S) but is
  # well defined for a single variable, where diag() would read its argument as
  # a matrix dimension rather than as a diagonal.
  R * tcrossprod(S)
}
