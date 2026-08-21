library(BayesPIM)

x <- cbind(1, c(-0.5, 0, 0.5))
beta <- c(0.2, -0.1)
logsigma <- log(0.8)

for (Q in c(-0.7, 0, 0.7)) {
  eta <- c(beta, logsigma, Q)
  par <- BayesPIM:::trans_par_gengamma(x, eta)

  stopifnot(
    isTRUE(all.equal(par[, "mu"], drop(x %*% beta))),
    isTRUE(all.equal(par[, "sigma"], rep(exp(logsigma), nrow(x)))),
    isTRUE(all.equal(par[, "Q"], rep(Q, nrow(x))))
  )

  evaluation_times <- c(0.5, 1, 2)
  probabilities <- BayesPIM:::pdist(
    evaluation_times,
    par = par,
    dist = "gengamma"
  )
  stopifnot(isTRUE(all.equal(
    probabilities,
    flexsurv::pgengamma(
      evaluation_times,
      mu = par[, "mu"],
      sigma = par[, "sigma"],
      Q = par[, "Q"]
    )
  )))

  stopifnot(isTRUE(all.equal(
    BayesPIM:::qdist(probabilities, par = par, dist = "gengamma"),
    evaluation_times
  )))
  stopifnot(isTRUE(all.equal(
    BayesPIM:::ddist(evaluation_times, par = par, dist = "gengamma"),
    flexsurv::dgengamma(
      evaluation_times,
      mu = par[, "mu"],
      sigma = par[, "sigma"],
      Q = par[, "Q"]
    )
  )))

  L <- c(0, 1, 1.5)
  R <- c(0.7, 2, Inf)
  log_likelihood <- BayesPIM:::ll_aft_ic(
    eta = eta,
    L = L,
    R = R,
    x = x,
    dist = "gengamma"
  )

  expected_log_likelihood <-
    flexsurv::pgengamma(
      R[1],
      mu = par[1, "mu"],
      sigma = par[1, "sigma"],
      Q = par[1, "Q"],
      log.p = TRUE
    ) +
    log(
      flexsurv::pgengamma(
        L[2],
        mu = par[2, "mu"],
        sigma = par[2, "sigma"],
        Q = par[2, "Q"],
        lower.tail = FALSE
      ) -
        flexsurv::pgengamma(
          R[2],
          mu = par[2, "mu"],
          sigma = par[2, "sigma"],
          Q = par[2, "Q"],
          lower.tail = FALSE
        )
    ) +
    flexsurv::pgengamma(
      L[3],
      mu = par[3, "mu"],
      sigma = par[3, "sigma"],
      Q = par[3, "Q"],
      lower.tail = FALSE,
      log.p = TRUE
    )

  stopifnot(isTRUE(all.equal(
    log_likelihood,
    expected_log_likelihood,
    tolerance = 1e-12
  )))

  prior <- log_aft_prior(
    eta = eta,
    dist = "gengamma",
    beta_prior = "norm",
    tau_t = 1.2,
    sig_prior = 0.9,
    q_prior_sd = 0.75
  )
  expected_prior <-
    sum(dnorm(beta, sd = 1.2, log = TRUE)) +
    log(2) +
    dnorm(exp(logsigma), sd = 0.9, log = TRUE) +
    logsigma +
    dnorm(Q, sd = 0.75, log = TRUE)
  stopifnot(isTRUE(all.equal(prior, expected_prior)))
}

for (Q in c(-0.5, 0, 0.5)) {
  set.seed(16)
  random_par <- cbind(
    mu = c(-0.2, 0, 0.2),
    sigma = rep(0.8, 3),
    Q = rep(Q, 3)
  )
  random_draws <- BayesPIM:::rdist(3, random_par, "gengamma")
  stopifnot(
    length(random_draws) == 3L,
    all(is.finite(random_draws)),
    all(random_draws > 0)
  )
}

# Posterior prediction must read the public chain columns as beta, positive
# sigma, and signed Q, and pass them to the Prentice generalized gamma.
x_ppd <- matrix(c(-0.5, 0, 0.5), ncol = 1)
chain_ppd <- coda::mcmc(rbind(
  c(0.2, -0.1, 0.8, -0.5),
  c(-0.1, 0.3, 0.6, 0.4)
))

set.seed(1601)
gengamma_ppd <- BayesPIM:::sample_ppd_nonprevalent(
  par_list = coda::mcmc.list(chain_ppd),
  x_t = x_ppd,
  dist = "gengamma",
  sampled_draws = c(1, 2)
)

linear_predictor_ppd <- cbind(1, x_ppd) %*%
  t(as.matrix(chain_ppd)[, 1:2, drop = FALSE])
sigma_ppd <- as.matrix(chain_ppd)[, 3]
q_ppd <- as.matrix(chain_ppd)[, 4]

set.seed(1601)
expected_gengamma_ppd <- vapply(
  seq_len(2),
  function(j) {
    flexsurv::rgengamma(
      nrow(x_ppd),
      mu = linear_predictor_ppd[, j],
      sigma = rep(sigma_ppd[j], nrow(x_ppd)),
      Q = rep(q_ppd[j], nrow(x_ppd))
    )
  },
  numeric(nrow(x_ppd))
)

stopifnot(
  identical(dim(gengamma_ppd), dim(expected_gengamma_ppd)),
  isTRUE(all.equal(gengamma_ppd, expected_gengamma_ppd)),
  all(is.finite(gengamma_ppd)),
  all(gengamma_ppd > 0)
)

for (Q in c(-0.5, 0, 0.5)) {
  set.seed(17)
  generated <- gen_data(
    n = 20,
    p = 1,
    beta_t = 0.1,
    beta_g = 0.1,
    theta = 0.1,
    mu_t = 1,
    sigma_t = 0.8,
    v_min = 0.5,
    v_max = 1,
    mean_rc = 4,
    dist = "gengamma",
    q = Q,
    prob_r = 1
  )
  stopifnot(
    length(generated$times_true) == 20L,
    all(is.finite(generated$times_true)),
    all(generated$times_true > 0)
  )
}

fake_run <- list(list(
  ini = matrix(c(0, log(0.8), -0.5), nrow = 1L),
  par = matrix(
    c(
      0.1, log(0.8), -0.5,
      0.2, log(0.9), 0.4
    ),
    nrow = 2L,
    byrow = TRUE
  ),
  terminal_par = matrix(c(0.2, log(0.9), 0.4), nrow = 1L),
  ac = NULL,
  times = NULL,
  k = 1L,
  g_aug = 0L,
  kappa = c(1, 1)
))

unwrapped <- BayesPIM:::unwrap_bayespim_chains(
  run = fake_run,
  times_rescale_factor = 1,
  p1_t = 1,
  x_t = NULL,
  x_g = NULL,
  prev = FALSE,
  has_q = TRUE,
  update_kappa = FALSE,
  sampler = "slice_collapsed"
)
unwrapped_matrix <- as.matrix(unwrapped$par[[1L]])
stopifnot(
  identical(colnames(unwrapped_matrix), c(
    "beta_t_intercept", "sigma_t", "q_t"
  )),
  isTRUE(all.equal(unwrapped_matrix[, "sigma_t"], c(0.8, 0.9))),
  isTRUE(all.equal(unwrapped_matrix[, "q_t"], c(-0.5, 0.4)))
)

fixed_q <- 0.75
initial_state <- list(
  ini_par_t = matrix(c(0.1, log(0.8), fixed_q), nrow = 1L),
  ini_par_g = matrix(0, nrow = 1L, ncol = 1L),
  ini_kappa = 0.8,
  ini_incidence_data = list(list(L = 1, R = 2, rows = 1L)),
  ini_g = matrix(0L, nrow = 1L, ncol = 1L),
  ini_k = matrix(1L, nrow = 1L, ncol = 1L)
)
gibbs_state <- BayesPIM:::ini_gibbs(
  update_run = FALSE,
  j = 1,
  inis = initial_state,
  p1_t = 1,
  p1_g = 1,
  v_obs = list(c(0, 1, 2)),
  times_prev = NULL,
  g_prev = NULL,
  k_prev = NULL,
  par_prev = NULL,
  update_kappa = FALSE,
  kappa = 0.8,
  dist = "gengamma",
  fix_sigma = FALSE,
  sig_prior = 1,
  fix_q = TRUE,
  q_prior_sd = fixed_q,
  prev = FALSE,
  sampler = "slice_collapsed",
  prop_sd = NULL,
  slice_width = 0.5,
  rescale_factor = 1
)
stopifnot(
  identical(gibbs_state$cur_par_t[1, 3], fixed_q),
  identical(gibbs_state$incidence_control$fix_q, TRUE)
)

validation_args <- list(
  v_obs = list(c(0, 1, Inf)),
  dist = "gengamma",
  kappa = 0.8,
  ndraws = 4,
  warmup = 1,
  chains = 1,
  sampler = "slice_collapsed",
  q_prior_sd = 0.75,
  fix_q = FALSE,
  prev = FALSE
)

expect_validation_error <- function(args, message) {
  error <- tryCatch(
    do.call(BayesPIM:::validate_bayespim_inputs, args),
    error = identity
  )
  stopifnot(
    inherits(error, "error"),
    grepl(message, conditionMessage(error), fixed = TRUE)
  )
}

# Both sampled and fixed Q are valid for the collapsed slice sampler.
stopifnot(isTRUE(invisible(do.call(
  BayesPIM:::validate_bayespim_inputs,
  validation_args
))))
stopifnot(isTRUE(invisible(do.call(
  BayesPIM:::validate_bayespim_inputs,
  modifyList(validation_args, list(fix_q = TRUE))
))))

# Generalized gamma is deliberately unavailable for the exact-time samplers.
for (unsupported_sampler in c("mh", "slice")) {
  args <- modifyList(
    validation_args,
    list(
      sampler = unsupported_sampler,
      prop_sd = if (unsupported_sampler == "mh") 0.1 else NULL
    )
  )
  expect_validation_error(
    args,
    '`dist = "gengamma"` requires `sampler = "slice_collapsed"`.'
  )
}

expect_validation_error(
  modifyList(validation_args, list(dist = "weibull", fix_q = TRUE)),
  '`fix_q = TRUE` requires `dist = "gengamma"` and `sampler = "slice_collapsed"`.'
)
expect_validation_error(
  modifyList(validation_args, list(fix_q = NA)),
  "`fix_q` must be a single TRUE or FALSE value."
)

for (invalid_q_prior_sd in list(0, -1, Inf, NA_real_, c(1, 2), "1")) {
  expect_validation_error(
    modifyList(
      validation_args,
      list(q_prior_sd = invalid_q_prior_sd)
    ),
    "`q_prior_sd` must be a single positive finite number."
  )
}
