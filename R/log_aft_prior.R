#' Log prior for accelerated failure time models
#'
#' Evaluates the log-prior density of the incidence-model parameters in an
#' accelerated failure time (AFT) model. For the Weibull, lognormal,
#' log-logistic, and gamma distributions, the parameter vector is
#'
#' \deqn{\eta = (\beta_0, \beta_1, \ldots, \beta_p, \log(\sigma)).}
#'
#' For the Prentice generalized-gamma distribution, the signed shape parameter
#' \eqn{Q} is appended without transformation:
#'
#' \deqn{\eta =
#'   (\beta_0, \beta_1, \ldots, \beta_p, \log(\sigma), Q).}
#'
#' The same prior-function contract is used for every supported distribution
#' and incidence sampler. A custom function supplied through
#' [bayespim()] must accept the named arguments `eta`, `dist`, `beta_prior`,
#' `tau_t`, `sig_prior`, and `q_prior_sd`, and must return one numeric
#' log-density value. Arguments that are irrelevant to a custom prior may be
#' ignored. When [bayespim()] is called with
#' `standardize_covariates = TRUE`, `eta` contains the internally standardized
#' incidence coefficients; returned posterior draws are transformed back to
#' the original covariate scale only after sampling.
#'
#' The default prior places either independent Student-\eqn{t} priors or
#' independent zero-centered normal priors on the regression coefficients.
#' A zero-centered half-normal prior is placed on the positive family
#' scale/dispersion parameter \eqn{\sigma}. For the gamma model,
#' \eqn{\sigma} is the conditional coefficient of variation; for the other
#' families it is their AFT scale parameter. The log density includes the
#' Jacobian for the transformation from \eqn{\sigma} to
#' \eqn{\log(\sigma)}. For the generalized-gamma model, \eqn{Q} additionally
#' receives a zero-centered normal prior with standard deviation
#' `q_prior_sd`.
#'
#' @param eta Numeric AFT parameter vector. It contains the intercept and slope
#'   coefficients followed by `log(sigma)`. For `dist = "gengamma"`, the
#'   untransformed signed shape parameter `Q` is the final element.
#' @param dist Character string identifying the AFT distribution. Supported
#'   values are `"weibull"`, `"lognormal"`, `"loglog"`, `"gamma"`, and
#'   `"gengamma"`.
#' @param beta_prior Character string specifying the prior family for the
#'   regression coefficients. Supported values are `"norm"` and `"t"`.
#' @param tau_t Numeric prior parameter. For `beta_prior = "t"`, this is the
#'   degrees of freedom of the Student-\eqn{t} prior. For
#'   `beta_prior = "norm"`, this is the standard deviation of the normal prior.
#' @param sig_prior Positive numeric standard deviation of the half-normal prior
#'   on \eqn{\sigma}. For `dist = "gamma"`, \eqn{\sigma} is the conditional
#'   coefficient of variation.
#' @param q_prior_sd Positive numeric standard deviation of the zero-centered
#'   normal prior on \eqn{Q}. It is ignored unless `dist = "gengamma"`.
#'
#' @return A single numeric value giving the log-prior density.
#'
#' @examples
#' # Use the default prior but give the incidence intercept five times its
#' # default prior scale. The function retains the complete prior contract and
#' # therefore works for all supported AFT distributions.
#' log_aft_prior_relaxed_intercept <- function(
#'     eta, dist, beta_prior, tau_t, sig_prior, q_prior_sd
#' ) {
#'   log_prior <- log_aft_prior(
#'     eta = eta,
#'     dist = dist,
#'     beta_prior = beta_prior,
#'     tau_t = tau_t,
#'     sig_prior = sig_prior,
#'     q_prior_sd = q_prior_sd
#'   )
#'
#'   if (beta_prior == "t") {
#'     log_prior -
#'       dt(eta[1], df = tau_t, log = TRUE) +
#'       dt(eta[1] / 5, df = tau_t, log = TRUE) -
#'       log(5)
#'   } else {
#'     log_prior -
#'       dnorm(eta[1], sd = tau_t, log = TRUE) +
#'       dnorm(eta[1], sd = 5 * tau_t, log = TRUE)
#'   }
#' }
#'
#' @export
log_aft_prior <- function(
    eta,
    dist,
    beta_prior,
    tau_t,
    sig_prior,
    q_prior_sd
) {
  has_q <- identical(dist, "gengamma")
  
  sigma_pos <- length(eta) - as.integer(has_q)
  
  if (sigma_pos < 2L) {
    stop("`eta` does not contain an AFT coefficient and log-sigma.")
  }
  
  beta <- eta[seq_len(sigma_pos - 1L)]
  logsigma <- eta[sigma_pos]
  sigma <- exp(logsigma)
  
  log_beta_prior <- switch(
    beta_prior,
    norm = sum(dnorm(beta, mean = 0, sd = tau_t, log = TRUE)),
    t = sum(dt(beta, df = tau_t, log = TRUE)),
    stop("Unsupported `beta_prior`: ", beta_prior)
  )
  
  # Half-normal prior on sigma, transformed to log-sigma.
  log_sigma_prior <-
    log(2) +
    dnorm(sigma, mean = 0, sd = sig_prior, log = TRUE) +
    logsigma
  
  log_q_prior <- 0
  
  if (has_q) {
    Q <- eta[sigma_pos + 1L]
    
    log_q_prior <- dnorm(
      Q,
      mean = 0,
      sd = q_prior_sd,
      log = TRUE
    )
  }
  
  log_beta_prior + log_sigma_prior + log_q_prior
}
