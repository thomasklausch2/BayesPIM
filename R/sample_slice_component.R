#' Update one parameter with univariate slice sampling
#' @noRd
sample_slice_component <- function(j, eta, w, log_ll, ...) {
  v0 <- eta[j]
  logy <- log_pst_component(
    v = v0, j = j, eta = eta, log_ll = log_ll, ...
  ) + log(stats::runif(1))

  l <- v0 - w * stats::runif(1)
  r <- l + w

  logl <- log_pst_component(v = l, j = j, eta = eta, log_ll = log_ll, ...)
  while (logl > logy) {
    l <- l - w
    logl <- log_pst_component(v = l, j = j, eta = eta, log_ll = log_ll, ...)
  }

  logr <- log_pst_component(v = r, j = j, eta = eta, log_ll = log_ll, ...)
  while (logr > logy) {
    r <- r + w
    logr <- log_pst_component(v = r, j = j, eta = eta, log_ll = log_ll, ...)
  }

  v1 <- stats::runif(1, l, r)
  logv1 <- log_pst_component(v = v1, j = j, eta = eta, log_ll = log_ll, ...)
  while (logv1 < logy) {
    if (v1 < v0) l <- v1 else r <- v1
    v1 <- stats::runif(1, l, r)
    logv1 <- log_pst_component(v = v1, j = j, eta = eta, log_ll = log_ll, ...)
  }

  v1
}
