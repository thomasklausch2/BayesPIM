#' @title Fitting Bayesian Prevalence-Incidence Mixture Model
#'
#' @description Estimates the prevalence-incidence mixture model of Klausch et
#' al. (2026) using a Bayesian Gibbs sampler. The model is formulated as an
#' interval-censored survival model over successive intervals, with the
#' possibility of missed events due to imperfect test sensitivity. In addition,
#' baseline tests at time zero may fail to detect pre-study events (prevalence).
#'
#' @details 
#' This Bayesian prevalence-incidence mixture model (PIM) characterizes time
#' to incidence through accelerated failure time (AFT) scaling. With
#' \eqn{\eta_i = \bm{x}_{ti}'\bm{\beta}_t}, covariates multiply event times by
#' \eqn{\exp(\eta_i)} relative to the corresponding baseline distribution.
#' For the Weibull, lognormal, and log-logistic families this is the familiar
#' log-location-scale representation
#'
#' \deqn{\log(t_i) = \eta_i + \sigma_t\epsilon_i.}{
#'   log(t_i) = eta_i + sigma_t * epsilon_i.}
#'
#' For \code{dist = "gamma"}, BayesPIM instead uses
#'
#' \deqn{t_i \mid \bm{x}_{ti} \sim
#'   \mathrm{Gamma}\{\sigma_t^{-2},
#'   \sigma_t^{-2}\exp(-\eta_i)\},}{
#'   t_i | x_ti ~ Gamma(shape = sigma_t^-2,
#'   rate = sigma_t^-2 * exp(-eta_i)),}
#'
#' so \eqn{E(t_i \mid \bm{x}_{ti}) = \exp(\eta_i)} and \eqn{\sigma_t} is the
#' conditional coefficient of variation. For \code{dist = "gengamma"},
#' BayesPIM uses the Prentice generalized gamma with location \eqn{\eta_i},
#' scale \eqn{\sigma_t}, and signed shape \eqn{Q}. Generalized-gamma fitting is
#' available only with \code{sampler = "slice_collapsed"}. The covariate matrix
#' for the latent-time model is supplied through \code{x_t}.
#' 
#' Baseline prevalence is modeled using a probit formulation \eqn{Pr(g_i=1 | \bm{x}_{gi}) = Pr(w_i > 0 | \bm{x}_{gi})} with
#' 
#' \deqn{w_i = \bm{x}_{gi}' \bm{\beta}_g + \psi_i}{w_i = x_{gi}' beta_g + psi_i}
#' 
#' where \eqn{\psi_i} follows a standard normal distribution, and the covariate vector \eqn{\bm{x}_{gi}} is given in the \code{x_g} matrix. The latent continuous probit variable \eqn{w_i} determines the modeled binary prevalence variable: \eqn{g_i = 1} if \eqn{w_i > 0} and \eqn{g_i = 0} otherwise.
#' 
#' The argument \code{v_obs} provides the observed testing times for all individuals. It is a list of numeric vectors, where each vector starts with \code{0} (representing the baseline time) and is followed by one or more screening times. The final entry is \code{Inf} in the case of right censoring or indicates the time of a positive test if an event is observed. Specifically:
#' \itemize{
#'   \item If the baseline test is positive, the vector consists solely of \code{c(0)}.
#'   \item If the baseline test is negative and right censoring occurs before the first regular screening, the vector is \code{c(0, Inf)}.
#'   \item Otherwise, the vector ends with \code{Inf} in the case of right censoring (e.g., \code{c(0, 1, 3, 6, Inf)}) or ends at the event time (e.g., \code{c(0, 1, 3, 6)} for an event detected at time \code{6}).
#' }
#' 
#' By convention, every vector in \code{v_obs} starts with \code{0}. However, the binary vector \code{r} of \code{length} \eqn{n} indicates whether the baseline test was conducted (\code{r[i] = 1}) or missing (\code{r[i] = 0}) for each individual \code{i} in \code{v_obs}. For further details on coding, see Section 2 of the main paper.
#' 
#' Test sensitivity can be fixed to a value \code{kappa} by setting
#' \code{update_kappa = FALSE}, or it can be estimated if
#' \code{update_kappa = TRUE}. When estimated with
#' \code{kappa_prior = c(m, s)}, a Beta prior with mean \eqn{m} and standard
#' deviation \eqn{s} is used. Its shape parameters are calculated
#' analytically as
#'
#' \deqn{\nu = \frac{m(1-m)}{s^2} - 1,\qquad
#'   a = m\nu,\qquad b = (1-m)\nu.}{
#'   nu = m * (1 - m) / s^2 - 1, a = m * nu, b = (1 - m) * nu.}
#'
#' The requested mean must lie strictly between zero and one, and the standard
#' deviation must satisfy \eqn{0 < s < \sqrt{m(1-m)}}. Malformed or infeasible
#' specifications stop with an informative error. The calculated shape
#' parameters must be finite and positive. If \code{kappa_prior = NULL}, the
#' function warns and uses an uninformative \eqn{\mathrm{Beta}(1,1)} prior.
#' In general, we advise against using an uninformative prior, but this default
#' avoids favoring any specific informative prior.
#'
#' The Gibbs sampler runs for \code{ndraws} iterations for each of
#' \code{chains} total chains. Incidence parameters can be updated by
#' Metropolis-Hastings, by univariate slice sampling conditional on augmented
#' exact incidence times, or by collapsed univariate slice sampling conditional
#' only on sampled screening intervals. The collapsed slice sampler
#' (\code{sampler = "slice_collapsed"}) is the default because it has shown
#' lower autocorrelation and faster convergence than the two exact-time
#' samplers. This differs from the original implementation described by
#' Klausch et al., which used a Metropolis-Hastings sampler. The collapsed
#' sampler also supports the generalized-gamma model, for which convergence
#' with the other samplers is typically slow; consequently,
#' \code{dist = "gengamma"} is available only with the collapsed sampler.
#'
#' When \code{sampler = "mh"}, the Metropolis step applies a normal proposal
#' distribution with standard deviation \code{prop_sd}, which must be selected
#' by trial and error. An optimal acceptance rate is approximately 23%, which
#' can be computed per MCMC run from the model output. The function
#' \link{search_prop_sd} provides a heuristic for selecting an effective
#' proposal standard deviation.
#' 
#' We recommend running at least two chains to facilitate standard MCMC diagnostics
#' such as the Gelman-Rubin statistic. For larger analyses, users may run more
#' chains, but should ensure that their computing environment permits the requested
#' parallel workers. CRAN examples use small sequential runs. 
#' Additionally, we suggest first running the sampler for a moderate number of iterations to assess its behavior before using the updating functionality in \code{prev_run} to extend sampling (see below).
#' 
#' The option \code{update_till_converge = TRUE} allows \code{bayespim} to run until convergence. Convergence is achieved when the rank-normalized split R-hat calculated by \code{posterior::rhat()} is at most \code{max_rhat} for every sampled parameter and the minimum effective sample size calculated by \code{posterior::ess_mean()} reaches \code{min_effss}. The diagnostics use every stored post-warm-up draw; no additional diagnostic thinning is applied. The current diagnostics are always calculated, stored, and printed. If automated updating is enabled, the sampler continues updating and printing diagnostics until convergence is attained or \code{maxit} is reached.
#'
#' Setting \code{silent = TRUE} suppresses this printing for both forms of updating, the automated loop just described and a manual update requested through \code{prev_run}. Only the progress output is withheld: the diagnostics are still calculated and stored in the \code{convergence} element of the returned object, and warnings and errors are still signaled.
#'
#' The priors for the regression coefficients and distributional parameters
#' can be controlled using \code{beta_prior}, \code{tau_t},
#' \code{sig_prior}, \code{tau_g}, and \code{q_prior_sd}. Specifically:
#' \itemize{
#'   \item \code{beta_prior} determines the prior type for \eqn{\beta_{tj}} (either \code{normal} or Student-\eqn{t} \code{t}).
#'   \item \code{tau_t} specifies either the standard deviation (for normal priors) or degrees of freedom (for Student-\eqn{t} priors). The default is a standard normal prior.
#'   \item A half-normal prior is used for the positive family
#'   scale/dispersion parameter \eqn{\sigma}, with \code{sig_prior}
#'   controlling the standard deviation. In the gamma model, \eqn{\sigma} is
#'   the conditional coefficient of variation.
#'   \item A zero-centered normal prior is assigned to \eqn{\beta_{gj}}, with \code{tau_g} controlling its standard deviation (default: standard normal).
#'   \item In the generalized-gamma model, a zero-centered normal prior is
#'   assigned to the signed shape parameter \eqn{Q}, with standard deviation
#'   \code{q_prior_sd}.
#' }
#' 
#' Sometimes model fitting can be improved by fixing the \eqn{\sigma}
#' parameter to a value, which is achieved through setting
#' \code{fix_sigma = TRUE}. Then, the value specified as \code{sig_prior} is
#' regarded as the fixed value of \eqn{\sigma}. For the gamma model this fixes
#' the conditional coefficient of variation. The functionality can also be
#' used to obtain the exponential distribution, akin to a Markov model. For
#' this choose \code{dist = "weibull"}, \code{sig_prior = 1}, and
#' \code{fix_sigma = TRUE}.
#' 
#' The \code{prev_run} argument allows updating a previous run with additional MCMC draws. The MCMC chain resumes from the last draws, continues, and merges with the original run. If an initial model was fit using \code{mod <- bayespim(...)}, it can be updated using \code{mod_update <- bayespim(prev_run = mod)}. By default, \code{ndraws} additional iterations are added unless otherwise specified via \code{ndraws_update}. If \code{warmup_updated = TRUE}, the warmup increases with each update; otherwise the initial warmup is retained.
#'
#' The Gibbs sampler requires starting values. Incidence-model coefficients are
#' anchored at estimates from an interval-censored survival model, while
#' prevalence-model coefficients are anchored at a penalized probit fit. With
#' multiple chains, the two coefficient vectors are multiplied by equally
#' spaced factors from \code{ini_spread} to \code{1}; a single chain uses the
#' unscaled fitted values. Each complete candidate state, including the latent
#' prevalence indicators and incidence intervals or times, must have a finite
#' posterior. Each chain makes up to ten initialization attempts. The first
#' eight redraw the latent states at that chain's own spread factor; the ninth
#' uses the fitted estimates unscaled; and the tenth uses zero-centered
#' coefficients. If none of the ten yields a finite posterior, initialization
#' stops with an informative error.
#' 
#'    
#' @param v_obs A list of length \eqn{n} of numeric vectors representing screening times. The first element of each vector should always be \code{0} and the last element \code{Inf} in the case of right censoring.
#' @param x_t A numeric matrix of dimension \eqn{n \times p_t} containing
#'   covariates for the AFT latent-time model. Missing values are not allowed.
#'   Categorical variables must be dummy-coded.
#' @param x_g A numeric matrix of dimension \eqn{n \times p_g} containing
#'   covariates for the probit prevalence model. Missing values are not allowed.
#'   Categorical variables must be dummy-coded.
#' @param r A binary vector of length \eqn{n} indicating whether a baseline test was conducted (\code{1} for yes, \code{0} for no / missing baseline test). 
#'
#' @param dist Character. Distribution for the time-to-incidence variable.
#'   Supported choices are \code{"weibull"}, \code{"lognormal"},
#'   \code{"loglog"} (log-logistic), \code{"gamma"}, and \code{"gengamma"}
#'   (Prentice generalized gamma). The generalized-gamma model is supported
#'   only with \code{sampler = "slice_collapsed"}.
#' @param kappa Numeric or \code{NULL}. Fixed test sensitivity used when
#'   \code{update_kappa = FALSE}; it must then be specified and lie in
#'   \eqn{(0,1]}. When \code{update_kappa = TRUE}, this argument is ignored
#'   because each chain receives internally generated starting values for
#'   \eqn{\kappa}.
#' @param update_kappa Logical. If \code{TRUE}, the test sensitivity (\eqn{\kappa}) is updated during the Gibbs sampler.
#' @param kappa_prior A numeric vector \code{c(mean, sd)} specifying the mean
#'   and standard deviation of a Beta prior for \eqn{\kappa}. The mean must lie
#'   strictly between zero and one, and \code{sd} must be positive and smaller
#'   than \eqn{\sqrt{\mathrm{mean}(1-\mathrm{mean})}}. Malformed or infeasible
#'   values produce an error. If \code{NULL}, a warning is issued and the
#'   uniform \eqn{\mathrm{Beta}(1,1)} prior is used.
#'
#' @param ndraws Integer. The total number of MCMC draws for the main Gibbs sampler.
#' @param warmup Integer. The number of initial generated MCMC iterations per
#'   chain omitted before posterior summaries and convergence diagnostics.
#'   This is always interpreted on the generated-iteration scale, independently
#'   of \code{save_every}. Defaults to half of \code{ndraws}, rounded down, and
#'   must be smaller than \code{ndraws} for an initial run.
#' @param warmup_updated Logical. If \code{TRUE}, each update adds the current
#'   warmup increment to the warmup stored in \code{prev_run}. If \code{FALSE},
#'   the stored warmup is retained.
#' @param prop_sd Numeric. The standard deviation for the proposal (jumping)
#'   distribution used when \code{sampler = "mh"}. It can be searched for
#'   heuristically using \link{search_prop_sd} and is not used by either slice
#'   sampler.
#' @param slice_width Numeric. The initial bracket width used by the univariate
#'   slice samplers. It affects computational efficiency but not the stationary
#'   distribution and is not used by \code{sampler = "mh"}.
#' @param chains Integer. The number of MCMC chains to run.
#' @param seed_chains Optional integer vector with one unique seed per chain. For
#'   a new fit, each seed initializes that chain's starting values and Gibbs
#'   sampler. If \code{NULL}, unique seeds are generated randomly. When
#'   \code{prev_run} is supplied, the saved end-of-chain RNG states are restored
#'   and \code{seed_chains} is ignored.
#' @param save_every Positive integer. Store the parameter state from every
#'   \code{save_every}-th generated iteration. The default, \code{1}, stores
#'   every draw. Values greater than one reduce the memory occupied by the
#'   returned chains and by chains accumulated during automatic updating, but
#'   permanently discard intervening parameter draws. Convergence diagnostics
#'   and summaries use every stored post-warm-up draw and never thin again.
#'   This setting is inherited by continued runs.
#' @param update_till_converge Logical. If \code{TRUE}, the model is updated iteratively until R-hat is at most \code{max_rhat} and effective sample size is at least \code{min_effss} for every sampled parameter. Diagnostics are calculated with the \code{posterior} package and printed whether this argument is \code{TRUE} or \code{FALSE}, unless \code{silent = TRUE}.
#' @param maxit A positive whole number or \code{Inf}. The maximum number of
#'   MCMC draws allowed before interrupting automatic convergence updates.
#'   \code{maxit} does not truncate the initially requested \code{ndraws}; a
#'   warning is given if an initial run requests more draws than a finite
#'   \code{maxit}. Default is \code{Inf} (i.e., no automatic interruption).
#' @param max_rhat Numeric. The maximum rank-normalized split R-hat accepted
#'   for every sampled parameter during convergence assessment and automatic
#'   updating. Defaults to \code{1.01}.
#' @param min_effss Integer. The minimum effective sample size required for each parameter before convergence is accepted during iterative updating.
#' @param sampler Character. Incidence-parameter update method. Use \code{"mh"}
#'   for Metropolis-Hastings, \code{"slice"} for univariate slice sampling
#'   conditional on augmented exact incidence times, or \code{"slice_collapsed"}
#'   to sample only the latent screening interval and update the incidence
#'   parameters with the interval-censored likelihood. The default is
#'   \code{"slice_collapsed"}. The generalized-gamma model requires this
#'   sampler.
#'
#' @param beta_prior Character. Specifies the type of prior for the latent-time regression coefficients (\eqn{\beta_{tj}}); options are \code{'norm'} for normal and \code{'t'} for student-t.
#' @param log_prior_fun Function used to evaluate the log-prior for the AFT
#'   incidence parameters. It must accept the named arguments \code{eta},
#'   \code{dist}, \code{beta_prior}, \code{tau_t}, \code{sig_prior}, and
#'   \code{q_prior_sd}, and return one numeric log-density value. For
#'   \code{dist = "gengamma"}, \code{eta} ends in
#'   \code{c(log(sigma), Q)}; otherwise it ends in \code{log(sigma)}.
#' @param tau_t Numeric. The hyperparameter for the prior distribution of the regression coefficients (\eqn{\beta_{tj}}) in the AFT latent-time model. For a normal prior, this is the standard deviation; for a student-t prior, it represents the degrees of freedom. The default produces a standard-normal prior.
#' @param sig_prior Numeric. Positive standard deviation of the half-normal
#'   prior on the family scale/dispersion parameter \eqn{\sigma}. For
#'   \code{dist = "gamma"}, \eqn{\sigma} is the conditional coefficient of
#'   variation. When \code{fix_sigma = TRUE}, this is instead the fixed value
#'   of \eqn{\sigma}.
#' @param tau_g Numeric. The hyperparameter (standard deviation) for the normal prior distribution of the regression coefficients (\eqn{\beta_{gj}}) in the probit prevalence model. The default produces a standard-normal prior.
#' @param fix_sigma Logical. If \code{TRUE}, the family scale/dispersion
#'   parameter \eqn{\sigma} is fixed at the value supplied through
#'   \code{sig_prior}; if \code{FALSE}, it is updated.
#' @param q_prior_sd Positive numeric standard deviation of the zero-centered
#'   normal prior on the signed Prentice generalized-gamma shape parameter
#'   \eqn{Q}. When \code{fix_q = TRUE}, \eqn{Q} is instead fixed at this value.
#'   It is otherwise ignored for distributions without \eqn{Q}.
#' @param fix_q Logical. If \code{TRUE}, fix the generalized-gamma shape
#'   parameter \eqn{Q} at \code{q_prior_sd}; if \code{FALSE}, update \eqn{Q}.
#'   This argument applies only when \code{dist = "gengamma"} and
#'   \code{sampler = "slice_collapsed"}.
#'
#' @param prev_run Optional. An unmodified object of class \code{bayespim}
#'   returned by a previous run. The complete fitted state is required, including
#'   the parameter chains, the terminal parameters, the latent prevalence and
#'   interval states, the covariate scaling, and the saved RNG state; a missing
#'   or malformed component is reported as an error before sampling begins.
#'   Data, sampler, priors, and RNG state are adopted from the previous run, and
#'   each chain continues from its last draw.
#' @param ndraws_update Integer greater than or equal to 2. The number of MCMC
#'   draws for updating a previous run or for convergence updates. If
#'   unspecified, \code{ndraws} is used.
#'
#' @param prev Logical. If \code{TRUE}, prevalence adjustment is applied; if \code{FALSE}, prevalence is assumed to be zero.
#'
#' @param par_exp Logical. If \code{TRUE}, the parameter expansion technique of Liu & Wu (1999) with a Haar prior is employed for updating the regression coefficients (\eqn{\beta_{wj}}) in the prevalence model. Experimental: tests suggest that it does not improve convergence or reduce autocorrelation.
#' @param rescale_times Logical. If \code{TRUE}, screening times are rescaled
#'   internally by the median latest finite observation time to improve numerical
#'   stability; returned times and incidence parameters are restored to the
#'   original scale.
#' @param standardize_covariates Logical. If \code{TRUE}, the default,
#'   non-binary columns of \code{x_t} and \code{x_g} are centered at their
#'   sample means and divided by their sample standard deviations for fitting.
#'   Columns with exactly two observed values are treated as binary/dummy
#'   variables and left unchanged. Returned coefficients are transformed back
#'   to the original covariate scale. The coefficient priors therefore apply
#'   on the standardized scale. This setting and the fitted centers and scales
#'   are inherited by continued runs.
#' @param ini_spread Numeric. Lower endpoint of the deterministic scaling
#'   factors used to disperse fresh-chain starting coefficients. With multiple
#'   chains, factors are equally spaced from \code{ini_spread} to \code{1},
#'   where \code{0} initializes coefficients at zero and \code{1} uses the
#'   fitted initialization estimates. With one chain, the fitted estimates are
#'   used directly. Values must lie between \code{0} and \code{1}.
#' @param silent Logical. If \code{TRUE}, suppress the progress information
#'   printed while the model is fitted and updated, namely the notice that a
#'   previous run is being updated and the convergence diagnostics reported
#'   after each set of draws. Warnings and errors are still signaled, and the
#'   diagnostics remain available in the \code{convergence} element of the
#'   returned object. Unlike most arguments, \code{silent} is never inherited
#'   from \code{prev_run}: it applies only to the current call, so an update of
#'   a silently fitted model prints again unless \code{silent = TRUE} is
#'   supplied afresh.
#'
#' @return A list containing the following elements:
#'
#' \item{ini}{A matrix with one row per chain holding the starting values of the
#' original fit. Coefficients are on the returned covariate and time scale, but
#' the latent-time scale parameter is stored as \eqn{\log(\sigma_t)} rather than
#' \eqn{\sigma_t}. A continued run retains the starting values of the run it
#' continues.}
#' \item{par}{An untrimmed \code{mcmc.list} containing the parameter draws
#' selected by \code{save_every}. With the default \code{save_every = 1}, all
#' generated draws are stored.}
#' \item{terminal_par}{A matrix containing the exact terminal parameter state
#' for each chain. It is used to continue a run when the terminal iteration is
#' not one of the iterations selected by \code{save_every}.}
#' \item{terminal_par_internal}{The corresponding terminal state on the
#' internal standardized-covariate and rescaled-time parameterization. It is
#' retained so updates can continue without reconstructing the sampler state.}
#' \item{times}{The terminal augmented latent times \eqn{t_i} for each chain. This is
#' \code{NULL} for \code{sampler = "slice_collapsed"}, which does not augment exact
#' latent times. For exact-time samplers, entries with terminal \code{g == 1}
#' are inactive stored values and should not be interpreted as current latent
#' incidence-time draws.}
#' \item{k}{The terminal latent interval index for every subject and chain.
#' Entries with terminal \code{g == 0} are the active sampled incidence
#' intervals.}
#' \item{g}{The terminal binary prevalence variables \eqn{g_i} for every subject,
#' concatenated by chain.}
#' \item{ac}{For \code{sampler = "mh"}, a matrix with stored MCMC draws in rows and
#' chains in columns, where each entry indicates acceptance (\code{1}) or
#' rejection (\code{0}). This is \code{NULL} for either slice sampler.}
#' \item{ac_cur}{For a continued run, the acceptance matrix from the most recent
#' update only. This is \code{NULL} for either slice sampler and is absent from
#' a fresh fit.}
#' \item{dat}{A data frame containing the last observed interval.}
#' \item{priors}{A list of prior specifications for the model parameters,
#' including \code{tau_t}, \code{sig_prior}, \code{tau_g}, and
#' \code{q_prior_sd}.}
#' \item{warmup}{The final warmup cutoff used by convergence diagnostics and,
#' by default, the summary and plot methods.}
#' \item{warmup_updated}{Whether warmup is increased during model updates.}
#' \item{seed_chains}{The unique seeds used to initialize the chains in the
#' original fit.}
#' \item{rng_state}{A list containing the RNG state saved at the end of each
#' chain. These states are restored when the fit is updated.}
#' \item{save_every}{The interval at which parameter draws were stored.}
#' \item{total_iterations}{The total number of iterations generated per chain,
#' including iterations not stored when \code{save_every > 1}.}
#' \item{covariate_scaling}{The centers, standard deviations, binary-column
#' indicators, and standardization indicators used for \code{x_t} and
#' \code{x_g}.}
#' \item{runtime}{The total runtime of the MCMC sampler.}
#' \item{max_rhat}{The maximum R-hat threshold used for convergence assessment
#' and automatic updating.}
#' \item{convergence}{A list containing the latest parameter-wise R-hat and ESS values, a printable diagnostics table, the convergence criteria, the number of draws and chains assessed, excluded fixed parameters, and the convergence status.}
#'
#' The returned list has class \code{"bayespim"}. Additionally, most input arguments are returned as part of the output for reference.
#'
#' @references 
#' T. Klausch, B. I. Lissenberg-Witte, and V. M. H. Coupé (2026). "A Bayesian prevalence-incidence mixture model for screening outcomes with misclassification.", Statistics in Medicine, 45(8-9), e70433. doi:10.1002/sim.70433
#' 
#' J. S. Liu and Y. N. Wu, “Parameter Expansion for Data Augmentation,” Journal of the American Statistical Association, vol. 94, no. 448, pp. 1264–1274, 1999, <doi:10.2307/2669940>.
#'
#' @seealso [summary.bayespim()] and [plot.bayespim()]
#' 
#' @examples
#' # A deliberately short fit for illustrating the interface.
#' # Use substantially more draws for scientific inference.
#' set.seed(2025)
#' dat <- gen_data(kappa = 0.7, n = 100, theta = 0.2,
#'                p = 1, p_discrete = 1,
#'                beta_t = c(0.2, 0.2), beta_g = c(0.2, 0.2),
#'                v_min = 20, v_max = 30, mean_rc = 80,
#'                sigma_t = 0.2, mu_t = 5, dist = "weibull",
#'                prob_r = 1)
#' 
#' fit <- bayespim(
#'   v_obs = dat$v_obs,
#'   x_t = dat$x,
#'   x_g = dat$x,
#'   r = dat$r,
#'   kappa = 0.7,
#'   update_kappa = FALSE,
#'   ndraws = 100,
#'   warmup = 10,
#'   chains = 2,
#'   seed_chains = c(202501, 202502),
#'   update_till_converge = FALSE,
#'   sampler = "slice_collapsed",
#'   dist = "weibull",
#'   silent = TRUE
#' )
#' fit$runtime
#' 
#' @name bayespim
#' @export
bayespim <- function(
    # Group 1: Data inputs
  v_obs, 
  x_t = NULL, 
  x_g = NULL, 
  r = NULL,
  
  # Group 2: Basic model settings
  dist = 'weibull', 
  kappa = NULL, 
  update_kappa = FALSE, 
  kappa_prior = NULL, 
  
  # Group 3: Main MCMC sampler settings
  ndraws = 1000, 
  warmup = floor(ndraws / 2),
  warmup_updated = FALSE,
  prop_sd = NULL, 
  slice_width = 1,
  chains , 
  seed_chains = NULL,
  save_every = 1,
  update_till_converge = FALSE, 
  maxit = Inf, 
  max_rhat = 1.01,
  min_effss = chains * 100, 
  sampler = 'slice_collapsed',
  
  # Group 4: Prior specifications
  log_prior_fun = log_aft_prior,
  beta_prior = 'norm', 
  tau_t = 1, 
  sig_prior = 1, 
  tau_g = 1, 
  fix_sigma = FALSE, 
  q_prior_sd = 1,
  fix_q = FALSE,
  
  # Group 5: Updating previous run settings
  prev_run = NULL, 
  ndraws_update = NULL, 
  
  # Group 6: Additional model flags
  prev = TRUE, 
  par_exp = FALSE, 
  rescale_times = TRUE,
  standardize_covariates = TRUE,
  ini_spread = 0.5,

  # Group 7: Console output
  silent = FALSE
) {
  
  t0 <- Sys.time()
  max_rhat_missing <- missing(max_rhat)
  iteration_offset <- 0L
  covariate_scaling <- NULL
  
  # Validate inputs
  validate_bayespim_inputs(
    v_obs = v_obs,
    x_t = x_t,
    x_g = x_g,
    r = r,
    dist = dist,
    kappa = kappa,
    update_kappa = update_kappa,
    kappa_prior = kappa_prior,
    ndraws = ndraws,
    warmup = warmup,
    warmup_updated = warmup_updated,
    prop_sd = prop_sd,
    slice_width = slice_width,
    chains = chains,
    seed_chains = seed_chains,
    save_every = save_every,
    update_till_converge = update_till_converge,
    maxit = maxit,
    max_rhat = max_rhat,
    min_effss = min_effss,
    sampler = sampler,
    log_prior_fun = log_prior_fun,
    beta_prior = beta_prior,
    tau_t = tau_t,
    sig_prior = sig_prior,
    tau_g = tau_g,
    fix_sigma = fix_sigma,
    q_prior_sd = q_prior_sd,
    fix_q = fix_q,
    prev_run = prev_run,
    ndraws_update = ndraws_update,
    prev = prev,
    par_exp = par_exp,
    rescale_times = rescale_times,
    standardize_covariates = standardize_covariates,
    ini_spread = ini_spread,
    silent = silent,
    stage = "initial"
  )
  
  # Import previous run
  if (!is.null(prev_run)) {
    # Stored chain state
    inis_prev = prev_run$ini
    previous_par = prev_run$par
    chains = length(previous_par)
    par_prev <- as.matrix(prev_run$terminal_par_internal)
    iteration_offset <- as.integer(prev_run$total_iterations)
    times_prev = prev_run$times
    k_prev = prev_run$k
    g_prev = prev_run$g
    prev_runtime = prev_run$runtime

    # Group 1: Data inputs
    v_obs = prev_run$v_obs
    x_t = prev_run$x_t
    x_g = prev_run$x_g
    r = prev_run$r

    # Group 2: Basic model settings
    dist = prev_run$dist
    kappa = prev_run$kappa[length(prev_run$kappa)]
    update_kappa = prev_run$update_kappa
    kappa_prior = prev_run$kappa_prior

    # Group 3: Main MCMC sampler settings
    if (is.null(ndraws_update)) {
      ndraws = prev_run$ndraws
    } else {
      ndraws = ndraws_update
    }
    warmup_imported = prev_run$warmup
    warmup_updated = prev_run$warmup_updated
    if (warmup_updated) {
      warmup = warmup + warmup_imported
    } else {
      warmup = warmup_imported
    }
    if (is.null(prop_sd)) {
      prop_sd = prev_run$prop_sd
    }
    slice_width = prev_run$slice_width
    seed_chains = prev_run$seed_chains
    save_every = prev_run$save_every
    if (is.infinite(maxit)) {
      maxit = prev_run$maxit
    }
    if (max_rhat_missing) {
      max_rhat = prev_run$max_rhat
    }
    sampler = prev_run$sampler

    # Group 4: Prior specifications
    log_prior_fun = prev_run$log_prior_fun
    beta_prior = prev_run$beta_prior
    tau_t = prev_run$priors$tau_t
    sig_prior = prev_run$priors$sig_prior
    tau_g = prev_run$priors$tau_g
    fix_sigma = prev_run$fix_sigma
    q_prior_sd = prev_run$priors$q_prior_sd
    fix_q = prev_run$fix_q

    # Group 5 controls the current update and is not inherited.

    # Group 6: Additional model flags
    prev = prev_run$prev
    par_exp = prev_run$par_exp
    rescale_times = prev_run$rescale_times
    standardize_covariates = prev_run$standardize_covariates
    covariate_scaling <- prev_run$covariate_scaling

    # `silent` is deliberately not imported from prev_run: it scopes to this call.
    if (!silent) message('Updating previous MCMC run. \n')
  }
  
  # Validate inputs after import
  validate_bayespim_inputs(
    v_obs = v_obs,
    x_t = x_t,
    x_g = x_g,
    r = r,
    dist = dist,
    kappa = kappa,
    update_kappa = update_kappa,
    kappa_prior = kappa_prior,
    ndraws = ndraws,
    warmup = warmup,
    warmup_updated = warmup_updated,
    prop_sd = prop_sd,
    slice_width = slice_width,
    chains = chains,
    seed_chains = seed_chains,
    save_every = save_every,
    update_till_converge = update_till_converge,
    maxit = maxit,
    max_rhat = max_rhat,
    min_effss = min_effss,
    sampler = sampler,
    log_prior_fun = log_prior_fun,
    beta_prior = beta_prior,
    tau_t = tau_t,
    sig_prior = sig_prior,
    tau_g = tau_g,
    fix_sigma = fix_sigma,
    q_prior_sd = q_prior_sd,
    fix_q = fix_q,
    prev_run = prev_run,
    ndraws_update = ndraws_update,
    prev = prev,
    par_exp = par_exp,
    rescale_times = rescale_times,
    standardize_covariates = standardize_covariates,
    ini_spread = ini_spread,
    silent = silent,
    stage = "effective"
  )
  
  # Determine kappa prior
  if(update_kappa & !is.null(kappa_prior)){
    f = find_ab(m=kappa_prior[1], s=kappa_prior[2])
    kappa_ab = c(f$a,f$b)
  }
  if(update_kappa & is.null(kappa_prior)) kappa_ab = c(1,1)
  
  prep <- prepare_bayespim_data(
    v_obs = v_obs,
    x_t = x_t,
    x_g = x_g,
    kappa = kappa,
    update_kappa = update_kappa,
    rescale_times = rescale_times,
    standardize_covariates = standardize_covariates,
    covariate_scaling = covariate_scaling,
    tau_g = tau_g,
    prev = prev
  )
  
  v_obs_original <- prep$v_obs_original
  v_obs <- prep$v_obs
  g_fixed <- prep$g_fixed
  pobs <- prep$pobs
  L <- prep$L
  R <- prep$R
  times_rescale_factor <- prep$times_rescale_factor
  v_obs_l <- prep$v_obs_l
  v_obs_r <- prep$v_obs_r
  x_t_original <- prep$x_t_original
  x_t <- prep$x_t
  x1_t <- prep$x1_t
  x_g_original <- prep$x_g_original
  x_g <- prep$x_g
  x1_g <- prep$x1_g
  covariate_scaling <- prep$covariate_scaling
  p1_t <- prep$p1_t
  p1_g <- prep$p1_g
  sig_inv <- prep$sig_inv
  sig_inv_xt <- prep$sig_inv_xt
  
  rm(prep)

  # Initialize independent streams for fresh chains.
  if(is.null(seed_chains)){
    seed_chains = sample.int(.Machine$integer.max, chains, replace = FALSE)
  } else seed_chains = as.integer(seed_chains)

  # Search starting values
  inis = NULL
  if(is.null(prev_run)){
    inis = ini_bayespim(v_obs, x_t, x1_t, x1_g, p1_t, r, g_fixed, chains, kappa, update_kappa, 
                        pobs, v_obs_l, v_obs_r, tau_g = tau_g, 
                        kappa_ab = if(update_kappa) kappa_ab else NULL, 
                        fix_sigma = fix_sigma,  seed_chains = seed_chains, 
                        sampler = sampler, prev = prev, log_prior_fun = log_prior_fun,
                        tau_t = tau_t, sig_prior = sig_prior, beta_prior = beta_prior,
                        dist = dist, spread_lower = ini_spread, q_prior_sd = q_prior_sd, fix_q = fix_q)
  }
  
  ## Start-up parallel run
  cl <- makePSOCKcluster(chains)
  cl_open <- TRUE
  
  on.exit({
    if (cl_open) try(stopCluster(cl), silent = TRUE)
    registerDoSEQ()
  }, add = TRUE)
  
  clusterSetRNGStream(cl)
  registerDoParallel(cl)

  batch_iterations <- iteration_offset + seq_len(ndraws)
  saved_iterations <- batch_iterations[batch_iterations %% save_every == 0L]
  save_positions <- match(saved_iterations, batch_iterations)

  run = foreach( j = seq_len(chains), .packages = c('BayesPIM') ) %dopar% {
    if(!is.null(prev_run)){
      assign(".Random.seed", prev_run$rng_state[[j]], envir = .GlobalEnv)
    } else {
      assign(".Random.seed", inis$rng_state[[j]], envir = .GlobalEnv)
    }
    
    init <- ini_gibbs(
      update_run = !is.null(prev_run),
      j = j,
      inis = inis,
      p1_t = p1_t,
      p1_g = p1_g,
      v_obs = v_obs,
      times_prev = if (!is.null(prev_run)) times_prev else NULL,
      g_prev = if (!is.null(prev_run)) g_prev else NULL,
      k_prev = if (!is.null(prev_run)) k_prev else NULL,
      par_prev = if (!is.null(prev_run)) par_prev else NULL,
      update_kappa = update_kappa,
      kappa = kappa,
      dist = dist,
      fix_sigma = fix_sigma,
      sig_prior = sig_prior,
      fix_q = fix_q,
      q_prior_sd = q_prior_sd,
      prev = prev,
      sampler = sampler,
      prop_sd = prop_sd,
      slice_width = slice_width,
      rescale_factor = times_rescale_factor
    )
    
    cur_par_t <- init$cur_par_t
    beta_g <- init$beta_g
    kappa_current <- init$kappa
    times <- init$times
    g_aug <- init$g_aug
    k <- init$k
    incidence_data <- init$incidence_data
    incidence_step <- init$incidence_step
    incidence_control <- init$incidence_control
    incidence_log_ll <- init$incidence_log_ll
    initial_par_t <- cur_par_t
    initial_beta_g <- beta_g
    n_saved <- length(saved_iterations)
    saved_par_t <- matrix(NA_real_, nrow = n_saved, ncol = ncol(cur_par_t))
    saved_beta_g <- matrix(NA_real_, nrow = n_saved, ncol = ncol(beta_g))
    saved_kappa <- numeric(n_saved)
    saved_ac <- if (sampler == "mh") numeric(n_saved) else NULL
    save_row <- 0L
    next_save_position <- if (n_saved > 0L) save_positions[1L] else Inf
    
    # Beginning of Gibbs
    for(i in seq_len(ndraws)){
      
      # Parameter update sweep
      incidence_update <- incidence_step(
        eta = cur_par_t,
        data = incidence_data,
        x = x1_t,
        dist = dist,
        log_ll = incidence_log_ll,
        log_prior_fun = log_prior_fun,
        tau_t = tau_t,
        sig_prior = sig_prior,
        beta_prior = beta_prior,
        q_prior_sd = q_prior_sd,
        fix_sigma = fix_sigma,
        control = incidence_control
      )

      cur_par_t[1L, ] <- incidence_update$eta
      accepted_current <- if (sampler == "mh") incidence_update$accepted else NULL
      
      # Transform latent-time parameters.
      if(dist == 'weibull' | dist == 'loglog') cur_dist_par_t = trans_par(x1_t, par = cur_par_t[1L,])
      if(dist == 'lognormal') cur_dist_par_t = trans_par_ind_norm(x1 = x1_t, p = cur_par_t[1L,1:p1_t], v= cur_par_t[1L,(p1_t+1)])
      if(dist == 'gamma') cur_dist_par_t = trans_par_gamma(x1 = x1_t, par = cur_par_t[1L,])
      if(dist == 'gengamma')  cur_dist_par_t = trans_par_gengamma(x1 = x1_t, par = cur_par_t[1L,])
      
      # Update pobs
      if(update_kappa) pobs = p_v_obs_rcpp(v_obs, kappa = kappa_current)

      # Calculate interval weights
      omega = interval_probs_rcpp(pobs_vec = unlist(pobs), v_obs, v_obs_l, v_obs_r,
                                  cur_dist_par_t, dist)
      interval_sums          = omega$interval_sums
      interval_probabilities = omega$pobs_norm
      
      # Calculate prob_g
      if(prev){
        mu_w    = as.numeric(x1_g %*% as.matrix(as.numeric(beta_g[1L,])))
        prob_g  = pnorm( mu_w )
        
        # Update g | theta_x, theta_g, D, collapsing over k; g_aug = 0 for all i if prev == FALSE
        g_aug <- augment_g_collapsed_rcpp(interval_sums = interval_sums, v_obs = v_obs, kappa = kappa_current,
                                          prob_g = prob_g, r = r, g_fixed = g_fixed )
        
        # Update w | g, x_g' beta_g.
        w_aug   = augment_w(g = g_aug, mu_w )
        
        # Update beta_g | w.
        if(par_exp){
          alpha_sq     = fc_beta_g_exp_haar(y = w_aug, x = x1_g, sig_inv_xt = sig_inv_xt) 
          beta_g[1L,] = fc_beta(x = x1_g, w_aug/sqrt(alpha_sq), sig_inv_xt = sig_inv_xt, sig_inv=sig_inv)
        } else{
          beta_g[1L,] = fc_beta(x = x1_g, w_aug, sig_inv_xt = sig_inv_xt, sig_inv=sig_inv)
        }
      }
      
      # Update k | theta_x, g, D
      incidence_rows <- which(g_aug == 0L)
      k_active <- sample_k_rcpp( interval_probabilities[incidence_rows] )
      phi_active <- look_up_mat_rcpp( v_obs[incidence_rows], k_active )
      k[incidence_rows] <- k_active
      
      # Update exact latent times for the conditional samplers.
      if (sampler == "slice_collapsed") {
        incidence_data <- list( L = phi_active[, 1L], R = phi_active[, 2L],rows = incidence_rows )
      } else if (sampler %in% c("mh", "slice")) {
        if (length(incidence_rows) > 0L) {
          times_active <- r_trdist( par = cur_dist_par_t[incidence_rows, , drop = FALSE],
                                         a = phi_active[, 1L], b = phi_active[, 2L], dist = dist ) 
          times_active[times_active == 0] <- 1e-300
          times[incidence_rows] <- times_active
        }
        
        incidence_data <- list(times = times, rows = incidence_rows)
      }
      
      # Update kappa
      if(update_kappa) {
        if(prev) kappa_current = fc_kappa_rcpp( v_obs, j_ = k, a=kappa_ab[1], b=kappa_ab[2], g=g_aug, r = r, g_fixed= g_fixed)
        if(!prev) kappa_current = pst_kappa_no_prev( v_obs, j_ = k, a=kappa_ab[1], b=kappa_ab[2])
      }

      if (i == next_save_position) {
        save_row <- save_row + 1L
        saved_par_t[save_row, ] <- cur_par_t[1L, ]
        if (prev) saved_beta_g[save_row, ] <- beta_g[1L, ]
        saved_kappa[save_row] <- kappa_current
        if (sampler == "mh") saved_ac[save_row] <- accepted_current
        next_save_position <- if (save_row < n_saved) {
          save_positions[save_row + 1L]
        } else {
          Inf
        }
      }
    }
    
    out=list()
    out$k = k
    if(sampler == "slice_collapsed") out["times"] = list(NULL) else out$times = times
    if(prev) {
      out$ini = cbind(initial_par_t, initial_beta_g)
      out$par = cbind(saved_par_t, saved_beta_g)
      out$terminal_par = cbind(cur_par_t, beta_g)
      out$g_aug = g_aug
    }
    if(!prev) {
      out$ini = initial_par_t
      out$par = saved_par_t
      out$terminal_par = cur_par_t
      out$g_aug = g_aug
    }
    out$kappa = saved_kappa
    if (update_kappa) {
      out$terminal_par <- cbind(out$terminal_par, kappa = kappa_current)
    }
    out$saved_iterations <- saved_iterations
    out["ac"] <- list( if (sampler == "mh") matrix(saved_ac, ncol = 1L) else NULL
    )
    out$rng_state = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    out
  }
  stopCluster(cl)
  cl_open <- FALSE
  registerDoSEQ()
  # end of gibbs
  
  # unwrap mcmc chains
  rng_state = lapply(run, `[[`, "rng_state")
  unwrapped <- unwrap_bayespim_chains(
    run = run,
    times_rescale_factor = times_rescale_factor,
    p1_t = p1_t,
    x_t = x_t,
    x_g = x_g,
    prev = prev,
    has_q = dist == 'gengamma',
    update_kappa = update_kappa,
    sampler = sampler,
    saved_iterations = saved_iterations,
    save_every = save_every,
    covariate_scaling = covariate_scaling,
    prev_run = prev_run,
    previous_par = if (!is.null(prev_run)) previous_par else NULL
  )
  
  inis <- unwrapped$inis
  par <- unwrapped$par
  terminal_par <- unwrapped$terminal_par
  terminal_par_internal <- unwrapped$terminal_par_internal
  ac <- unwrapped$ac
  times_draw <- unwrapped$times_draw
  k_draw <- unwrapped$k_draw
  g_draw <- unwrapped$g_draw
  
  if (!is.null(prev_run)) {
    ac_cur <- unwrapped$ac_cur
  }
  
  rm(unwrapped)

  # Undo v_obs recoding
  v_obs <- v_obs_original
  L <- L * times_rescale_factor
  R <- R * times_rescale_factor
  if(sampler != "slice_collapsed") times_draw <- times_draw * times_rescale_factor
  
  # Runtime stamp
  t1 = Sys.time()
  runtime = t1-t0
  if(!is.null(prev_run)) runtime = runtime + prev_runtime
  
  # Export fitted state
  ret = list()
  if (is.null(prev_run)) {
    ret$ini = inis
  } else {
    ret$ini = inis_prev
  }
  ret$par = par
  ret$terminal_par = terminal_par
  ret$terminal_par_internal = terminal_par_internal
  if (sampler == "slice_collapsed") {
    ret["times"] = list(NULL)
  } else {
    ret$times = times_draw
  }
  ret$k = k_draw
  ret$g = g_draw
  ret["ac"] <- list(ac)
  if (!is.null(prev_run)) ret["ac_cur"] <- list(ac_cur)

  # Group 1: Data inputs
  ret$v_obs = v_obs
  ret["x_t"] = list(x_t_original)
  ret["x_g"] = list(x_g_original)
  ret["r"] = list(r)

  # Derived data summary
  ret$dat = data.frame(L=L, R=R)

  # Group 2: Basic model settings
  ret$dist = dist
  ret["kappa"] = list(kappa)
  ret$update_kappa = update_kappa
  ret["kappa_prior"] = list(kappa_prior)

  # Group 3: Main MCMC sampler settings
  ret$ndraws = ndraws
  ret$warmup = warmup
  ret$warmup_updated = warmup_updated
  ret["prop_sd"] = list(prop_sd)
  ret$slice_width = slice_width
  ret$seed_chains = seed_chains
  ret$rng_state = rng_state
  ret$save_every = save_every
  ret$total_iterations = iteration_offset + ndraws
  ret$maxit = maxit
  ret$max_rhat = max_rhat
  ret$sampler = sampler

  # Group 4: Prior specifications
  ret$log_prior_fun = log_prior_fun
  ret$beta_prior = beta_prior
  ret$priors = list(
    tau_t = tau_t,
    sig_prior = sig_prior,
    tau_g = tau_g,
    q_prior_sd = q_prior_sd
  )
  ret$fix_sigma = fix_sigma
  ret$fix_q = fix_q

  # Group 5 arguments are call-only and are not stored in the fit.

  # Group 6: Additional model flags
  ret$prev = prev
  ret$par_exp = par_exp
  ret$rescale_times = rescale_times
  ret$standardize_covariates = standardize_covariates
  ret$covariate_scaling = covariate_scaling

  # Other fit metadata
  ret$runtime = runtime
  class(ret) <- c("bayespim", class(ret))

  internal_sampler_call <- !is.null(prev_run) && isTRUE(attr(prev_run, "bayespim_internal_convergence_update"))

  if (!internal_sampler_call) {
    ndraws_update_convergence <- if (is.null(ndraws_update)) {
      ndraws
    } else {
      ndraws_update
    }

    convergence_state <- handle_bayespim_convergence(
      ret = ret,
      min_effss = min_effss,
      maxit = maxit,
      ndraws_update = ndraws_update_convergence,
      max_rhat = max_rhat,
      update_till_converge = update_till_converge,
      silent = silent
    )
    ret <- convergence_state$ret

    while (isTRUE(convergence_state$needs_update)) {
      internal_prev_run <- ret
      attr(
        internal_prev_run,
        "bayespim_internal_convergence_update"
      ) <- TRUE

      ret <- bayespim(
        prev_run = internal_prev_run,
        ndraws_update = convergence_state$ndraws_next,
        update_till_converge = FALSE,
        maxit = maxit,
        max_rhat = max_rhat,
        silent = silent
      )

      convergence_state <- handle_bayespim_convergence(
        ret = ret,
        min_effss = min_effss,
        maxit = maxit,
        ndraws_update = ndraws_update_convergence,
        max_rhat = max_rhat,
        update_till_converge = TRUE,
        silent = silent
      )
      ret <- convergence_state$ret
    }
  }
  return(ret)
}
