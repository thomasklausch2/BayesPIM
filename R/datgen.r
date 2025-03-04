#' @title gen.dat: Simulate Screening Data for a Prevalence-Incidence Mixture Model
#'
#' @description 
#' Generates synthetic data according to the Bayesian prevalence-incidence mixture (PIM) framework of Klausch et al. (2025) with interval-censored screening outcomes. 
#' The function simulates continuous or discrete baseline covariates, event times from one of several parametric families, and irregular screening schedules, 
#' yielding interval-censored observations suitable for testing or demonstrating PIM-based or other interval-censored survival methods.
#'
#' @details 
#' The data-generating process includes:
#'
#' \enumerate{
#'   \item \strong{Covariates \eqn{Z}:}
#'         Continuous covariates are simulated using a correlation structure specified by \code{r} and a common standard deviation \code{s}. 
#'         If \code{p.discrete = 1}, a single discrete covariate is added, drawn from \eqn{\mathrm{Bernoulli}(0.5)}.
#'
#'   \item \strong{Event Times \eqn{X}:}
#'         An Accelerated Failure Time (AFT) model is used:
#'         \deqn{\log(x_i) = \beta_{x0} + \beta_{x}^\top z_{xi} + \sigma_X \,\epsilon_i,}{
#'               log(x_i) = beta_{x0} + beta_{x}' z_{xi} + sigma_X * epsilon_i,}
#'         where \eqn{\beta_{x0}} is the intercept (set by \code{mu.X}) and \eqn{\beta_{x}} are the other regression coefficients (provided via \code{beta.X}). 
#'         The error term \eqn{\epsilon_i} is drawn from the distribution chosen by \code{dist.X}: 
#'         \code{"weibull"}, \code{"lognormal"}, \code{"loglog"} (log-logistic), or \code{"gengamma"} (generalized gamma). 
#'         For \code{"gengamma"}, the shape parameter \code{k} is additionally used.
#'
#'   \item \strong{Irregular Screening Schedules \eqn{V_i}:}
#'         Each individual has multiple screening times generated randomly between \code{v.min} and \code{v.max}, 
#'         ending in right censoring or the time of detection. 
#'         These screening times (including a 0 for baseline and \code{Inf} for censoring) are returned in \code{Vobs}.
#'
#'   \item \strong{Prevalence Indicator \eqn{g_i}:}
#'         Baseline prevalence is modeled via either a probit or logit link, consistent with:
#'         \deqn{w_i = \beta_{w0} + \beta_{w}^\top z_{wi} + \psi_i,}{
#'               w_i = beta_{w0} + beta_{w}' z_{wi} + noise,}
#'         where \eqn{\beta_{w0}} is determined by \code{theta}, and \eqn{\beta_{w}} by \code{beta.W}. 
#'         Specifically:
#'         \itemize{
#'           \item If \code{sel.mod = "probit"}, then \eqn{\beta_{w0} = \mathrm{qnorm}(\theta)}.
#'           \item If \code{sel.mod = "logit"}, then \eqn{\beta_{w0} = \log(\theta / (1-\theta))}.
#'         }
#'         We set \eqn{g_i = 1} if \eqn{w_i > 0}, and \eqn{g_i = 0} otherwise.
#'
#'   \item \strong{Baseline Test Missingness \eqn{r_i}:}
#'         A baseline test indicator \eqn{r_i \in \{0,1\}} is generated via \eqn{\mathrm{Bernoulli}(\text{prob.r})}, 
#'         so \eqn{r_i = 1} means the baseline test is performed and \eqn{r_i = 0} means it is missing.
#'
#'   \item \strong{Test Sensitivity \eqn{\kappa}:}
#'         A misclassification parameter \eqn{\kappa} (test sensitivity) can be specified via \code{kappa}. 
#'         If \eqn{\kappa < 1}, some truly positive cases are missed.
#' }
#'
#' @param kappa Numeric. Test sensitivity parameter \eqn{\kappa} used when generating misclassification. A value of 1 implies perfect sensitivity.
#' @param n Integer. Sample size.
#' @param p Integer. Number of continuous baseline covariates to simulate.
#' @param p.discrete Integer. If \code{1}, include an additional discrete covariate \eqn{Z_{\mathrm{discrete}}} from \eqn{\mathrm{Bernoulli}(0.5)}; otherwise, none.
#' @param r Numeric. Correlation coefficient(s) used to build the covariance matrix of continuous covariates. If \code{p > 1}, 
#'          off-diagonal entries of the correlation matrix are set to \code{r}.
#' @param s Numeric. Standard deviation(s) of the continuous covariates. If \code{p > 1}, all continuous covariates share the same \code{s}.
#' @param sigma.X Numeric. Scale parameter \eqn{\sigma_X} in the AFT model for \eqn{\log(x_i)}.
#' @param mu.X Numeric. Intercept \eqn{\beta_{x0}} in the AFT model. In the linear predictor, it appears as 
#'             \eqn{\log(x_i) = \beta_{x0} + \beta_{x}^\top Z_i + \sigma_X \epsilon_i}. 
#'             Practically, \code{mu.X} is prepended to \code{beta.X} when forming the full parameter vector.
#' @param beta.X Numeric vector. The coefficients \eqn{\beta_{x}} for the AFT model. 
#'               Combined with \code{mu.X}, the log-scale model is 
#'               \code{cbind(1, Z_i) \%*\% c(mu.X, beta.X)}.
#' @param beta.W Numeric vector. The coefficients \eqn{\beta_{w}} for the prevalence model. 
#'               The intercept \eqn{\beta_{w0}} is derived from \code{theta}.
#' @param theta Numeric. Baseline prevalence parameter on the probability scale. Under:
#' \itemize{
#'   \item \code{sel.mod = "probit"}: \eqn{\beta_{w0} = \mathrm{qnorm}(\theta)}.
#'   \item \code{sel.mod = "logit"}:  \eqn{\beta_{w0} = \log(\theta / (1 - \theta))}.
#' }
#' @param v.min Numeric. Minimum spacing for irregular screening intervals.
#' @param v.max Numeric. Maximum spacing for irregular screening intervals.
#' @param mean.rc Numeric. Mean of the exponential distribution controlling a random right-censoring time \eqn{t_{\mathrm{rc}}} after the first screening.
#' @param dist.X Character. Distribution for survival times \eqn{x_i}: \code{"weibull"}, \code{"lognormal"}, \code{"loglog"} (log-logistic), or \code{"gengamma"} (generalized gamma).
#' @param k Numeric. Shape parameter for \code{"gengamma"} only.
#' @param sel.mod Character. Either \code{"probit"} or \code{"logit"}, specifying the link function for the prevalence model.
#' @param prob.r Numeric. Probability that a baseline test is performed (\eqn{r_i = 1}). If \code{prob.r = 0}, no baseline tests are done.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{Vobs}}{A list of length \code{n}, each entry containing screening times. 
#'                      The first element is 0 (baseline), and \code{Inf} may indicate right censoring.}
#'   \item{\code{X.true}}{Numeric vector of length \code{n} giving the true (latent) event times \eqn{x_i}.}
#'   \item{\code{Z}}{Numeric matrix of dimension \eqn{n \times p} (plus an extra column if \code{p.discrete = 1}) containing the covariates.}
#'   \item{\code{C}}{Binary vector of length \code{n}, indicating whether an individual is truly positive at baseline (\eqn{g_i = 1}).}
#'   \item{\code{r}}{Binary vector of length \code{n}, indicating whether the baseline test was performed (\eqn{r_i = 1}) or missing (\eqn{r_i = 0}).}
#'   \item{\code{p.W}}{Numeric vector of length \code{n} giving the true prevalence probabilities, \eqn{P(g_i = 1)}.}
#' }
#' 
#' @references 
#' T. Klausch, B. I. Lissenberg-Witte, and V. M. Coupé, “A Bayesian prevalence-incidence mixture model for screening outcomes with misclassification,” arXiv:2412.16065.
#' 
#' @examples
#' # Generate a small dataset for testing
#' set.seed(2025)
#' sim_data <- gen.dat(n = 100, p = 1, p.discrete = 1,
#'                     sigma.X = 0.5, mu.X = 2,
#'                     beta.X = c(0.2, 0.2), beta.W = c(0.5, -0.2),
#'                     theta = 0.2,
#'                     dist.X = "weibull", sel.mod = "probit")
#' str(sim_data)
#' 
#' @export
gen.dat <- function(kappa = 0.7,
                    n = 1000,
                    p = 2,
                    p.discrete = 0,
                    r = 0,
                    s = 1,
                    sigma.X = 1/2,
                    mu.X = 4,
                    beta.X = NULL,
                    beta.W = NULL,
                    theta  = 0.15,
                    v.min = 1,
                    v.max = 6,
                    mean.rc = 40,
                    dist.X = 'weibull',
                    k = 1,
                    sel.mod = 'probit',
                    prob.r  = 0) {

  # Sim Z
  R <- matrix(r, p, p)
  diag(R) <- 1
  S <- rep(s, p)
  Sigma <- cor2cov(R, S)
  if (p > 0) {
    Z <- mvrnorm(n, mu = rep(0, p), Sigma)
  } else {
    Z <- NULL
  }

  # Sim discrete Z
  if (p.discrete == 1) {
    Z.discrete <- rbinom(n, 1, 0.5)
    Z <- cbind(Z, Z.discrete)
    colnames(Z) <- paste(1:ncol(Z))
  }

  Z1 <- cbind(as.matrix(rep(1, n)), Z)

  # Sim.X
  if (dist.X == 'weibull' | dist.X == 'loglog' | dist.X == 'lognormal') {
    if (dist.X == 'weibull') {
      e.X <- r.ev(n)
    }
    if (dist.X == 'loglog') {
      e.X <- rlog(n)
    }
    if (dist.X == 'lognormal') {
      e.X <- rnorm(n)
    }
    X <- Z1 %*% c(mu.X, beta.X) + sigma.X * e.X  # log surv times
    X <- exp(X)  # surv times
    X <- as.numeric(X)
  }

  # Simulate X using the generalized gamma distribution
  if (dist.X == 'gengamma') {
    # Parameters for generalized gamma distribution
    lambda <- exp(-Z1 %*% c(mu.X, beta.X))
    gamma <- 1 / sigma.X

    # Generate samples from the generalized gamma distribution
    X <- rggamma(n, a = lambda^(-1), b = gamma, k = k)
  }

  # Screening times
  # Generate screening sequences
  V <- list()
  t.rc <- numeric()
  for (i in 1:n) {
    v_i <- runif(1, v.min, v.max)
    t.rc[i] <- v_i + rexp(1, 1 / mean.rc)
    j <- 1
    while (v_i[j] < t.rc[i]) {
      v_i[j+1] <- runif(1, v_i[j] + v.min, v_i[j] + v.max)
      j <- j + 1
    }
    v_i <- v_i[-length(v_i)]
    v_i <- c(0, v_i, Inf)
    V[[i]] <- v_i
  }

  if(sel.mod == 'probit'){
  mu.W   = Z1 %*% as.matrix(c(qnorm(theta), beta.W)) + rnorm(n)
  g <- as.numeric(mu.W > 0)}
  if(sel.mod == 'logit'){
  g <-rbinom(n,1, 1/(1+exp(-(Z1 %*% as.matrix(c( log(theta/(1-theta)), beta.W))) )))
  }

  r = rbinom(n, 1, prob.r)
  Vobs <- v_to_vobs(V, X, kappa, C=g, baseline.test = r )
  X.true <- X

  ret <- list(Vobs = Vobs, X.true = X.true, Z = Z, C = g, r = r, p.W = pnorm(Z1 %*% as.matrix(c(qnorm(theta), beta.W))))
  ret
}
