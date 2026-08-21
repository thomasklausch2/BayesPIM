#' Log-likelihood for interval-censored accelerated failure time models
#' @noRd
ll_aft_ic <- function(eta, L, R, x, dist) {
  if (length(L) == 0L) {
    return(0)
  }

  par <- switch(
    dist,
    weibull = trans_par(x, eta),
    loglog = trans_par(x, eta),
    lognormal = trans_par_ind_norm(
      p = eta[-length(eta)],
      v = eta[length(eta)],
      x1 = x
    ),
    gamma = trans_par_gamma(x, eta),
    gengamma = trans_par_gengamma(x, eta),
    stop("Unsupported distribution: ", dist)
  )

  if (any(!is.finite(par))) {
    return(-Inf)
  }

  interval_censored <- L > 0 & is.finite(R)
  right_censored <- L > 0 & is.infinite(R)
  left_censored <- L <= 0

  logsurvL <- pdist(
    L[interval_censored],
    par = par[interval_censored, , drop = FALSE],
    dist = dist,
    lower_tail = FALSE,
    log_p = TRUE
  )
  logsurvR <- pdist(
    R[interval_censored],
    par = par[interval_censored, , drop = FALSE],
    dist = dist,
    lower_tail = FALSE,
    log_p = TRUE
  )
  
  interval_ll <- logsdiff(logsurvL, logsurvR)
  if (anyNA(interval_ll)) return(-Inf)
  ll_interval <- sum(interval_ll)

  ll_left <- sum(pdist(
    R[left_censored],
    par = par[left_censored, , drop = FALSE],
    dist = dist,
    lower_tail = TRUE,
    log_p = TRUE
  ))

  ll_right <- sum(pdist(
    L[right_censored],
    par = par[right_censored, , drop = FALSE],
    dist = dist,
    lower_tail = FALSE,
    log_p = TRUE
  ))

  ll_interval + ll_left + ll_right
}
