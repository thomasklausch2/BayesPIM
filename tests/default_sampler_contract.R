library(BayesPIM)

stopifnot(
  identical(formals(bayespim)$sampler, "slice_collapsed"),
  identical(
    formals(BayesPIM:::validate_bayespim_inputs)$sampler,
    "slice_collapsed"
  )
)

set.seed(44)
sampler_data <- gen_data(
  n = 20,
  p = 1,
  beta_t = 0.1,
  beta_g = 0.1,
  theta = 0.1,
  mu_t = 1,
  sigma_t = 0.5,
  v_min = 0.5,
  v_max = 1,
  mean_rc = 4,
  prob_r = 1
)

default_fit <- bayespim(
  v_obs = sampler_data$v_obs,
  x_t = sampler_data$x,
  x_g = sampler_data$x,
  r = sampler_data$r,
  kappa = 0.8,
  ndraws = 4,
  warmup = 1,
  chains = 1,
  q_prior_sd = 0.75,
  min_effss = 1,
  silent = TRUE
)
stopifnot(
  identical(default_fit$sampler, "slice_collapsed"),
  identical(
    names(default_fit$priors),
    c("tau_t", "sig_prior", "tau_g", "q_prior_sd")
  ),
  identical(default_fit$priors$q_prior_sd, 0.75),
  !("q_prior_sd" %in% names(default_fit)),
  default_fit$save_every == 1L,
  isTRUE(default_fit$standardize_covariates),
  is.list(default_fit$covariate_scaling),
  default_fit$total_iterations == 4L,
  nrow(as.matrix(default_fit$par[[1L]])) == 4L,
  nrow(default_fit$terminal_par) == 1L,
  nrow(default_fit$terminal_par_internal) == 1L,
  is.null(default_fit$times),
  is.null(default_fit$ac)
)

default_summary <- NULL
summary_output <- capture.output(
  default_summary <- summary(default_fit)
)
stopifnot(
  identical(default_summary$sampler, "slice_collapsed"),
  any(summary_output == "Incidence sampler: slice_collapsed")
)

continued_fit <- bayespim(
  prev_run = default_fit,
  ndraws_update = 2,
  min_effss = 1,
  silent = TRUE
)
stopifnot(
  identical(continued_fit$sampler, "slice_collapsed"),
  identical(continued_fit$priors$q_prior_sd, 0.75),
  isTRUE(continued_fit$standardize_covariates),
  identical(continued_fit$covariate_scaling, default_fit$covariate_scaling),
  !("q_prior_sd" %in% names(continued_fit)),
  is.null(continued_fit$times),
  is.null(continued_fit$ac)
)
