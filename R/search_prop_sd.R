#' Automated Heuristic Search of a Proposal Standard Deviation for \code{bayespim}
#'
#' When \link{bayespim} is fitted with \code{sampler = "mh"}, it uses a
#' Metropolis step for sampling the incidence-model parameters and requires a
#' standard deviation for the normal proposal distribution. This function uses
#' a heuristic algorithm to find a proposal standard deviation such that the
#' Metropolis sampler accepts proposed draws at a rate within the user-defined
#' interval (by default around 20--25%). The default sampler in
#' \code{bayespim()} is \code{"slice_collapsed"}, so the Metropolis sampler
#' must be requested explicitly before using this function.
#'
#' @details
#' Starting from an initial \code{bayespim} model object \code{m}, the function
#' attempts to calibrate the standard deviation of the proposal distribution.
#' Specifically, it:
#' \enumerate{
#'   \item Runs an initial update of \code{ndraws} iterations and computes an
#'     acceptance rate.
#'   \item If the acceptance rate lies within \code{acc_bounds}, the number
#'     of MCMC draws \code{ndraws} is doubled, and the process repeats.
#'   \item Otherwise, the proposal standard deviation \eqn{\sigma} is adjusted
#'     based on whether the acceptance rate \eqn{p} is below the lower bound \eqn{a}
#'     or above the upper bound \eqn{b} of \code{acc_bounds}.
#'   \item The formula for adjustment is:
#'     \deqn{\sigma \leftarrow \sigma \times (1 - \frac{ (a-p)}{a}) \quad\text{if } p < a, \quad
#'           \sigma \leftarrow \sigma \times (1 + \frac{ (p-b)}{b}) \quad\text{if } p > b.}
#' }
#' By default, if the acceptance rate falls within \eqn{[0.2, 0.25]}, that \eqn{\sigma}
#' is considered acceptable, and the process continues until \code{succ_min} consecutive
#' successes (doubles) are achieved. 
#'
#' @param m A model object of class \code{bayespim}.
#' @param ndraws Starting number of MCMC iterations after which the acceptance rate
#'   is first evaluated. Defaults to 1000.
#' @param succ_min The algorithm doubles the number of MCMC draws \code{succ_min} times
#'   (each time the acceptance rate is within \code{acc_bounds}), ensuring
#'   stability. Defaults to 3.
#' @param acc_bounds A numeric vector of length two specifying the lower and upper
#'   bounds for the acceptable acceptance rate. Defaults to \code{c(0.2, 0.25)}.
#'
#' @return A \code{list} with the following elements:
#' \describe{
#'   \item{\code{prop_sd}}{The final (adjusted) proposal standard deviation.}
#'   \item{\code{ac}}{The acceptance rate in the last iteration.}
#' }
#'
#' @examples
#' \dontrun{
#' # search_prop_sd() requires an MH fit. This deliberately short initial fit
#' # illustrates the interface; use more draws for an actual calibration.
#' set.seed(2025)
#' dat <- gen_data(
#'   kappa = 0.7, n = 100, theta = 0.2,
#'   p = 1, p_discrete = 1,
#'   beta_t = c(0.2, 0.2), beta_g = c(0.2, 0.2),
#'   v_min = 20, v_max = 30, mean_rc = 80,
#'   sigma_t = 0.2, mu_t = 5, dist = "weibull",
#'   prob_r = 1
#' )
#'
#' fit_mh <- bayespim(
#'   v_obs = dat$v_obs, x_t = dat$x, x_g = dat$x, r = dat$r,
#'   kappa = 0.7, update_kappa = FALSE,
#'   ndraws = 100, warmup = 10, chains = 2,
#'   seed_chains = c(202501, 202502),
#'   update_till_converge = FALSE,
#'   sampler = "mh", prop_sd = 0.005, dist = "weibull",
#'   silent = TRUE
#' )
#'
#' search_sd <- search_prop_sd(
#'   m = fit_mh,
#'   ndraws = 100,
#'   succ_min = 1
#' )
#' print(search_sd)
#' }
#'
#' @export
search_prop_sd = function(m, ndraws = 1000, succ_min = 3, acc_bounds =c(0.2,0.25)){
  if (!identical(m$sampler, "mh")) {
    stop("`search_prop_sd()` is only available for models fitted with `sampler = \"mh\"`.", call. = FALSE)
  }
  found = FALSE
  it = 1; succ = 0;
  while(succ!=succ_min){
    message(paste('Iteration',it,'\n'))
    if(it == 1) { ac_cur = mean(m$ac)
    prop_sd = m$prop_sd}
    if(it >1) {
      m = bayespim( prev_run = m, ndraws_update = ndraws, prop_sd = prop_sd,
                    silent = TRUE)
      ac_cur = mean(m$ac_cur)
    }
    message( paste('Acceptance rate was:', round(ac_cur, 3),'\n' ))
    ac = ac_cur
    if((ac > acc_bounds[1] & ac < acc_bounds[2]) ){
      found = TRUE} else{
      if(ac < acc_bounds[1]) {
        dif =  1-(acc_bounds[1] - ac)/acc_bounds[1]
        prop_sd=prop_sd * dif
      }
      if(ac > acc_bounds[2]) {
        dif =  1+(ac - acc_bounds[2])/acc_bounds[2]
        prop_sd=prop_sd * dif
      }
      found = FALSE
      message(paste('prop_sd is set to', round(prop_sd,3),'\n'))
    }
    it=it+1
    if(found ){
      ndraws = ndraws*2
      succ = succ +1
      found = FALSE
      if(succ!=succ_min) message(paste('Success. Doubling number of MCMC draws:',ndraws,'\n'))
      }
  }
  message('Finished calibrating proposal variance. \n')
  ret= list()
  ret$prop_sd = prop_sd
  ret$ac      = ac
  ret
}
