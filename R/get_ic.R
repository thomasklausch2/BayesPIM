#' Compute Information Criteria for a Bayesian Prevalence-Incidence Mixture Model
#'
#' Computes and returns information criteria for a fitted Bayesian prevalence-incidence mixture model, including the Widely Applicable Information Criterion 1 (WAIC-1), WAIC-2, and the Deviance Information Criterion (DIC). These criteria are commonly used for model comparison and evaluation in Bayesian analysis. See Gelman et al. (2014) for further details on these criteria.
#'
#' @details
#' This function calculates information criteria for a fitted Bayesian prevalence-incidence mixture model (`bayespim`). The information criteria include:
#'
#' - **WAIC-1**: Based on the sum of posterior variances of log-likelihood contributions.
#' - **WAIC-2**: Similar to WAIC-1 but incorporates an alternative variance estimate.
#' - **DIC**: Measures model fit by penalizing complexity via the effective number of parameters.
#'
#' The computation is performed by evaluating log-likelihood values for MCMC samples. By default, all MCMC samples after burn-in are used, though a subset can be specified via the `samples` argument. All incidence-time distributions supported by [bayespim()] (Weibull, log-logistic, log-normal, gamma, and generalized gamma) are handled.
#' For a model fitted with covariate standardization, likelihood calculations
#' reconstruct the standardized parameterization and operate on standardized
#' covariates to reduce cancellation on large original scales; this
#' is algebraically equivalent to using the returned original-scale
#' coefficients.
#'
#' Parallelization is available via the `foreach` package, utilizing multiple cores if `cores` is set accordingly. If `cores = NULL`, all available cores will be used.
#'
#' @param mod A fitted prevalence-incidence mixture model of class \code{bayespim}.
#' @param samples The number of MCMC samples to use. By default, all draws
#' available after the model's warm-up are used. If the model was fitted with
#' \code{save_every > 1}, this means all stored post-warm-up draws.
#' @param cores The number of cores for parallel processing using \code{foreach}. If \code{NULL} (default), all available cores will be used.
#'
#' @return A \code{matrix} containing WAIC-1, WAIC-2, and DIC values for the model.
#'
#' @references
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. Stat Comput, 24(6), 997–1016.
#'
#' @examples
#' data(mod)
#' set.seed(2025)
#' get_ic(mod, samples = 20, cores = 1)
#'
#' @export
get_ic = function(mod, samples = NULL, cores = NULL){
  update_kappa = mod$update_kappa
  p1_t = bayespim_design_ncol(mod$x_t) + 1L
  p1_g = bayespim_design_ncol(mod$x_g) + 1L
  covariate_scaling <- mod$covariate_scaling

  if(is.null(cores)) cores = detectCores()
  par_matrix = bayespim_analysis_chains(mod)
  par_matrix = as.matrix(par_matrix)
  if (is.null(samples)) samples = nrow(par_matrix)
  if(!update_kappa) par_matrix = cbind(par_matrix, mod$kappa)
  if (samples > nrow(par_matrix)) {
    stop(
      "`samples` cannot exceed the number of stored post-warm-up posterior draws (",
      nrow(par_matrix),
      ").",
      call. = FALSE
    )
  }

  # `l_obs_2s()` expects log(sigma); the chains store sigma on its natural
  # scale. The generalized-gamma shape Q (column p1_t + 2) is left untouched.
  # Work with standardized covariates and their matching coefficients. This
  # avoids evaluating large original-scale covariate/intercept combinations.
  par_matrix <- bayespim_transform_coefficients(
    par = par_matrix,
    p1_t = p1_t,
    p1_g = p1_g,
    prev = mod$prev,
    has_q = mod$dist == "gengamma",
    covariate_scaling = covariate_scaling,
    direction = "public_to_internal"
  )
  mod_computational <- mod
  mod_computational$x_t <- bayespim_standardize_matrix(
    mod$x_t,
    covariate_scaling$x_t,
    "x_t"
  )
  mod_computational$x_g <- bayespim_standardize_matrix(
    mod$x_g,
    covariate_scaling$x_g,
    "x_g"
  )
  par_matrix[,(p1_t + 1)] = log(par_matrix[,(p1_t + 1)])
  m_s = par_matrix[sample(nrow(par_matrix), samples, replace = FALSE), , drop = FALSE]

  # Spread the selected draws across the workers, covering every draw exactly
  # once (the previous split silently dropped the final draw).
  n_workers = min(cores, samples)
  chunks = split(seq_len(samples), (seq_len(samples) - 1L) %% n_workers + 1L)

  cl    = makePSOCKcluster(n_workers)
  cl_open <- TRUE

  # Release the workers and the parallel backend even if the loop below fails.
  on.exit({
    if (cl_open) try(stopCluster(cl), silent = TRUE)
    registerDoSEQ()
  }, add = TRUE)

  clusterSetRNGStream(cl)
  registerDoParallel(cl)
  posterior_mean = apply(par_matrix, 2, mean)

  # Only the BayesPIM namespace needs to be loaded on each worker; it exposes
  # `l_obs_2s()` and every distribution helper it dispatches to.
  run = foreach(idx = chunks, .packages = c('BayesPIM')) %dopar% {
    m_s_cores = m_s[idx, , drop = FALSE]
    apply(m_s_cores, 1, function(x) l_obs_2s(est = x, mod = mod_computational, log_scale = FALSE, sumup = FALSE))
  }
  stopCluster(cl)
  cl_open <- FALSE
  registerDoSEQ()

  run = do.call(cbind, run)

  lppd = sum (log(apply(run,1, mean)))
  q1   = sum (apply( log(run), 1, mean))
  q2   = sum( apply( log(run), 1, var) )
  q3   = mean( apply( log(run),2, sum) )
  q4   = sum(l_obs_2s(posterior_mean, mod=mod_computational, log_scale = TRUE, sumup=FALSE))
  DIC  = -2* ( q4 - 2*(q4-q3) )
  WAIC1    = -2*(-lppd + 2*q1)
  WAIC2    = -2*(lppd - q2)
  mat=matrix(nrow=1, ncol=3)
  mat[1,]=c(WAIC1, WAIC2, DIC)
  colnames(mat) = c('WAIC1', 'WAIC2', 'DIC')
  mat
}
