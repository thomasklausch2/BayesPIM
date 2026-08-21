library(BayesPIM)

set.seed(904)
scaling_data <- gen_data(
  n = 30,
  p = 1,
  p_discrete = 1,
  beta_t = c(0.02, 0.2),
  beta_g = c(0.03, 0.2),
  theta = 0.15,
  mu_t = 2,
  sigma_t = 0.4,
  v_min = 0.5,
  v_max = 1,
  mean_rc = 5,
  prob_r = 1
)
x_original <- scaling_data$x
x_original[, 1L] <- 50 + 20 * x_original[, 1L]
colnames(x_original) <- c("age", "exposed")

common_args <- list(
  v_obs = scaling_data$v_obs,
  x_t = x_original,
  x_g = x_original,
  r = scaling_data$r,
  kappa = 0.8,
  ndraws = 12,
  warmup = 4,
  chains = 1,
  seed_chains = 905,
  min_effss = 1,
  silent = TRUE
)

automatic_fit <- do.call(bayespim, common_args)
scaling <- automatic_fit$covariate_scaling
stopifnot(
  isTRUE(automatic_fit$standardize_covariates),
  identical(automatic_fit$x_t, x_original),
  identical(automatic_fit$x_g, x_original),
  identical(unname(scaling$x_t$binary), c(FALSE, TRUE)),
  identical(unname(scaling$x_t$standardized), c(TRUE, FALSE)),
  isTRUE(all.equal(scaling$x_t$center, c(age = mean(x_original[, 1L]), exposed = 0))),
  isTRUE(all.equal(scaling$x_t$scale, c(age = stats::sd(x_original[, 1L]), exposed = 1))),
  identical(scaling$x_t, scaling$x_g)
)

x_standardized <- BayesPIM:::bayespim_standardize_matrix(
  x_original,
  scaling$x_t,
  "x_t"
)
stopifnot(
  isTRUE(all.equal(x_standardized[, 1L], as.numeric(scale(x_original[, 1L])))),
  identical(x_standardized[, 2L], x_original[, 2L])
)

# Fitting a manually standardized design with standardization disabled must
# produce the same internal trajectory as automatic standardization.
manual_args <- common_args
manual_args$x_t <- x_standardized
manual_args$x_g <- x_standardized
manual_args$standardize_covariates <- FALSE
manual_fit <- do.call(bayespim, manual_args)
stopifnot(
  !manual_fit$standardize_covariates,
  !any(manual_fit$covariate_scaling$x_t$standardized),
  identical(automatic_fit$terminal_par_internal, manual_fit$terminal_par_internal),
  identical(automatic_fit$rng_state, manual_fit$rng_state)
)

p1_t <- ncol(x_original) + 1L
p1_g <- ncol(x_original) + 1L
manual_public <- BayesPIM:::bayespim_transform_coefficients(
  par = as.matrix(manual_fit$par[[1L]]),
  p1_t = p1_t,
  p1_g = p1_g,
  prev = TRUE,
  has_q = FALSE,
  covariate_scaling = scaling,
  direction = "internal_to_public"
)
stopifnot(isTRUE(all.equal(
  unname(as.matrix(automatic_fit$par[[1L]])),
  unname(manual_public),
  tolerance = 1e-12
)))

# Original- and standardized-scale linear predictors must agree. The stable
# helper deliberately performs the calculation on the standardized scale.
automatic_matrix <- as.matrix(automatic_fit$par[[1L]])
manual_matrix <- as.matrix(manual_fit$par[[1L]])
lp_t_automatic <- BayesPIM:::bayespim_scaled_linear_predictor(
  x_original,
  automatic_matrix[, seq_len(p1_t), drop = FALSE],
  scaling$x_t
)
lp_t_manual <- cbind(1, x_standardized) %*%
  t(manual_matrix[, seq_len(p1_t), drop = FALSE])
first_g <- p1_t + 2L
lp_g_automatic <- BayesPIM:::bayespim_scaled_linear_predictor(
  x_original,
  automatic_matrix[, seq.int(first_g, length.out = p1_g), drop = FALSE],
  scaling$x_g
)
lp_g_manual <- cbind(1, x_standardized) %*%
  t(manual_matrix[, seq.int(first_g, length.out = p1_g), drop = FALSE])
stopifnot(
  isTRUE(all.equal(lp_t_automatic, lp_t_manual, tolerance = 1e-11)),
  isTRUE(all.equal(lp_g_automatic, lp_g_manual, tolerance = 1e-11))
)

# IC and posterior prediction use the standardized computational scale even
# though the automatic fit exposes original-scale covariates and coefficients.
set.seed(906)
ic_automatic <- get_ic(automatic_fit, samples = 6, cores = 1)
set.seed(906)
ic_manual <- get_ic(manual_fit, samples = 6, cores = 1)
stopifnot(isTRUE(all.equal(ic_automatic, ic_manual, tolerance = 1e-10)))

set.seed(907)
cif_automatic <- ppCIF(
  automatic_fit,
  pst_samples = 6,
  ppd_type = "percentiles",
  quant = c(0, 1, 2)
)
set.seed(907)
cif_manual <- ppCIF(
  manual_fit,
  pst_samples = 6,
  ppd_type = "percentiles",
  quant = c(0, 1, 2)
)
stopifnot(
  isTRUE(all.equal(cif_automatic$mixture, cif_manual$mixture, tolerance = 1e-10)),
  isTRUE(all.equal(
    cif_automatic$nonprevalent,
    cif_manual$nonprevalent,
    tolerance = 1e-10
  ))
)

# Both divisible and non-divisible storage/update paths continue from the exact
# internal terminal state; scaling constants are inherited, not recomputed.
continued_fit <- bayespim(
  prev_run = automatic_fit,
  ndraws_update = 4,
  min_effss = 1,
  silent = TRUE
)
uninterrupted_args <- common_args
uninterrupted_args$ndraws <- 16
uninterrupted_fit <- do.call(bayespim, uninterrupted_args)
stopifnot(
  identical(continued_fit$covariate_scaling, automatic_fit$covariate_scaling),
  identical(continued_fit$terminal_par_internal, uninterrupted_fit$terminal_par_internal),
  isTRUE(all.equal(
    as.matrix(continued_fit$par[[1L]]),
    as.matrix(uninterrupted_fit$par[[1L]]),
    tolerance = 1e-12
  ))
)

zero_variance_error <- try(
  do.call(
    bayespim,
    within(common_args, {
      x_t <- cbind(constant = rep(1, nrow(x_original)))
      x_g <- x_original
    })
  ),
  silent = TRUE
)
stopifnot(
  inherits(zero_variance_error, "try-error"),
  grepl("zero-variance column(s): constant", as.character(zero_variance_error), fixed = TRUE)
)

invalid_switch <- try(
  do.call(bayespim, c(common_args, list(standardize_covariates = 1))),
  silent = TRUE
)
stopifnot(
  inherits(invalid_switch, "try-error"),
  grepl(
    "`standardize_covariates` must be a single TRUE or FALSE value.",
    as.character(invalid_switch),
    fixed = TRUE
  )
)
