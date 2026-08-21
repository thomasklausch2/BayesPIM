library(BayesPIM)

expected_exports <- c(
  "bayespim",
  "gen_data",
  "get_ic",
  "ppCIF",
  "log_aft_prior",
  "search_prop_sd",
  "trim_mcmc"
)
stopifnot(setequal(getNamespaceExports("BayesPIM"), expected_exports))

bayespim_args <- names(formals(bayespim))
validator_args <- names(formals(BayesPIM:::validate_bayespim_inputs))
stopifnot(
  identical(formals(bayespim)$sampler, "slice_collapsed"),
  identical(
    formals(BayesPIM:::validate_bayespim_inputs)$sampler,
    "slice_collapsed"
  ),
  identical(
    validator_args[seq_along(bayespim_args)],
    bayespim_args
  ),
  identical(
    validator_args[-seq_along(bayespim_args)],
    "stage"
  ),
  match("fix_sigma", bayespim_args) < match("q_prior_sd", bayespim_args),
  match("q_prior_sd", bayespim_args) < match("fix_q", bayespim_args),
  match("fix_q", bayespim_args) < match("prev_run", bayespim_args)
)
# `validate_bayespim_inputs()` mirrors the `bayespim()` signature, so shared
# defaults must not drift apart. The four arguments below are the deliberate
# exception: the validator defaults them to NULL and uses NULL as a "not
# supplied" sentinel, which `bayespim()` cannot express for an argument that is
# either required or whose default depends on another argument.
sentinel_defaults <- c("v_obs", "chains", "warmup", "min_effss")
bayespim_formals <- formals(bayespim)
validator_formals <- formals(BayesPIM:::validate_bayespim_inputs)
deparse_default <- function(x) paste(deparse(x), collapse = "")
shared_defaults <- setdiff(
  intersect(names(bayespim_formals), names(validator_formals)),
  sentinel_defaults
)
drifted_defaults <- shared_defaults[vapply(
  shared_defaults,
  function(nm) {
    !identical(
      deparse_default(bayespim_formals[[nm]]),
      deparse_default(validator_formals[[nm]])
    )
  },
  logical(1)
)]
stopifnot(
  length(shared_defaults) > 0L,
  length(drifted_defaults) == 0L,
  identical(bayespim_formals$ini_spread, validator_formals$ini_spread),
  all(vapply(
    sentinel_defaults,
    function(nm) is.null(validator_formals[[nm]]),
    logical(1)
  ))
)

required_args <- c(
  "v_obs", "x_t", "x_g", "update_kappa", "kappa_prior",
  "prop_sd", "slice_width", "save_every", "update_till_converge", "max_rhat",
  "beta_prior", "tau_t", "tau_g", "fix_sigma", "prev_run",
  "ndraws_update", "par_exp", "rescale_times", "standardize_covariates",
  "q_prior_sd", "fix_q"
)
stopifnot(all(required_args %in% bayespim_args))
stopifnot(!any(grepl("[.]", setdiff(bayespim_args, "..."))))

invalid_sampler_fit <- try(
  bayespim(
    v_obs = replicate(6, c(0, 1, Inf), simplify = FALSE),
    x_t = matrix(0, nrow = 6, ncol = 1),
    r = rep(1, 6),
    kappa = 1,
    ndraws = 2,
    warmup = 1,
    chains = 1,
    sampler = c("mh", "slice"),
    prev = FALSE
  ),
  silent = TRUE
)
stopifnot(
  inherits(invalid_sampler_fit, "try-error"),
  grepl(
    "`sampler` must be one of 'mh', 'slice', or 'slice_collapsed'.",
    as.character(invalid_sampler_fit),
    fixed = TRUE
  )
)

expected_prior_args <- c(
  "eta", "dist", "beta_prior", "tau_t", "sig_prior", "q_prior_sd"
)
stopifnot(identical(names(formals(log_aft_prior)), expected_prior_args))

expected_step_prior_args <- c(
  "log_prior_fun", "tau_t", "sig_prior", "beta_prior", "q_prior_sd"
)
for (step_fun in list(
  BayesPIM:::step_mh,
  BayesPIM:::step_slice,
  BayesPIM:::step_slice_collapsed
)) {
  stopifnot(all(expected_step_prior_args %in% names(formals(step_fun))))
}

stopifnot("fix_q" %in% names(formals(BayesPIM:::ini_bayespim)))

set.seed(42)
generated <- gen_data(
  n = 30,
  p = 1,
  beta_t = 0.1,
  beta_g = 0.1,
  theta = 0.1,
  mu_t = 2,
  sigma_t = 0.3,
  v_min = 1,
  v_max = 2,
  mean_rc = 8,
  prob_r = 1
)
expected_generated_names <- c("v_obs", "times_true", "x", "g", "r", "prob_g")
stopifnot(setequal(names(generated), expected_generated_names))

stopifnot(identical(
  names(formals(BayesPIM:::augment_g_collapsed_rcpp)),
  c("interval_sums", "v_obs", "kappa", "prob_g", "r", "g_fixed")
))
stopifnot(identical(
  names(formals(BayesPIM:::look_up_mat_rcpp)),
  c("v_obs", "interval_indices")
))
