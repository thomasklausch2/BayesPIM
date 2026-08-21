library(BayesPIM)

beta_shapes <- BayesPIM:::find_ab(m = 0.7, s = 0.1)
stopifnot(
  isTRUE(all.equal(beta_shapes$a, 14)),
  isTRUE(all.equal(beta_shapes$b, 6)),
  isTRUE(all.equal(beta_shapes$check_mean, 0.7)),
  isTRUE(all.equal(beta_shapes$check_sd, 0.1)),
  all(is.finite(c(beta_shapes$a, beta_shapes$b))),
  all(c(beta_shapes$a, beta_shapes$b) > 0)
)

validation_args <- list(
  v_obs = replicate(6, c(0, 1, Inf), simplify = FALSE),
  x_t = matrix(0, nrow = 6, ncol = 1),
  r = rep(1, 6),
  update_kappa = TRUE,
  ndraws = 4,
  warmup = 1,
  chains = 1,
  min_effss = 1,
  sampler = "slice_collapsed",
  prev = FALSE
)

warnings <- character()
validated <- withCallingHandlers(
  do.call(BayesPIM:::validate_bayespim_inputs, validation_args),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  isTRUE(validated),
  length(warnings) == 1L,
  grepl("Uniform(0, 1) prior", warnings, fixed = TRUE)
)

stopifnot(isTRUE(do.call(
  BayesPIM:::validate_bayespim_inputs,
  c(validation_args, list(kappa_prior = c(0.7, 0.1)))
)))

expect_error_message <- function(expr, pattern) {
  error <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = identity
  )
  stopifnot(
    inherits(error, "error"),
    grepl(pattern, conditionMessage(error), fixed = TRUE)
  )
}

expect_error_message(
  do.call(
    BayesPIM:::validate_bayespim_inputs,
    c(validation_args, list(kappa_prior = "invalid"))
  ),
  "`kappa_prior` must be NULL or a numeric vector of length 2"
)
expect_error_message(
  do.call(
    BayesPIM:::validate_bayespim_inputs,
    c(validation_args, list(kappa_prior = c(1, 0.1)))
  ),
  "prior mean for kappa, must lie strictly between 0 and 1"
)
expect_error_message(
  do.call(
    BayesPIM:::validate_bayespim_inputs,
    c(validation_args, list(kappa_prior = c(0.7, -0.1)))
  ),
  "prior standard deviation for kappa, must be positive"
)
expect_error_message(
  do.call(
    BayesPIM:::validate_bayespim_inputs,
    c(validation_args, list(kappa_prior = c(0.7, sqrt(0.7 * (1 - 0.7)))))
  ),
  "`kappa_prior` is not feasible for a Beta prior"
)
expect_error_message(
  BayesPIM:::find_ab(m = 0.7, s = .Machine$double.xmin),
  "not finite and positive"
)
