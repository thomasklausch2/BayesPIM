library(BayesPIM)

x <- cbind(1, c(-0.5, 0, 0.5))
beta <- c(0.2, -0.1)
sigma <- 0.8
eta <- drop(x %*% beta)

par <- BayesPIM:::trans_par_gamma(
  x1 = x,
  par = c(beta, log(sigma))
)

stopifnot(
  identical(colnames(par), c("shape", "rate")),
  isTRUE(all.equal(
    par[, "shape"],
    rep(sigma^(-2), nrow(x))
  )),
  isTRUE(all.equal(
    par[, "rate"],
    sigma^(-2) * exp(-eta)
  )),
  isTRUE(all.equal(
    par[, "shape"] / par[, "rate"],
    exp(eta)
  )),
  isTRUE(all.equal(
    1 / sqrt(par[, "shape"]),
    rep(sigma, nrow(x))
  ))
)

y <- c(0.5, 1, 2)
gamma_log_likelihood <- BayesPIM:::ll_aft(
  eta = c(beta, log(sigma)),
  y = y,
  x = x,
  dist = "gamma"
)
expected_gamma_log_likelihood <- sum(stats::dgamma(
  y,
  shape = par[, "shape"],
  rate = par[, "rate"],
  log = TRUE
))

stopifnot(
  isTRUE(all.equal(
    gamma_log_likelihood,
    expected_gamma_log_likelihood
  )),
  identical(
    BayesPIM:::ll_aft(
      eta = c(beta, log(sigma)),
      y = c(0, 1, 2),
      x = x,
      dist = "gamma"
    ),
    -Inf
  )
)

# Posterior prediction must use the gamma mean/CV parameterization rather than
# the Weibull-style scale/shape transformation used by the other positive-time
# families.
x_ppd <- matrix(c(-0.5, 0, 0.5), ncol = 1)
chain_ppd <- coda::mcmc(rbind(
  c(0.2, -0.1, 0.8),
  c(-0.1, 0.3, 0.6)
))

set.seed(3301)
gamma_ppd <- BayesPIM:::sample_ppd_nonprevalent(
  par_list = coda::mcmc.list(chain_ppd),
  x_t = x_ppd,
  dist = "gamma",
  sampled_draws = c(1, 2)
)

linear_predictor_ppd <- cbind(1, x_ppd) %*%
  t(as.matrix(chain_ppd)[, 1:2, drop = FALSE])
sigma_ppd <- as.matrix(chain_ppd)[, 3]

set.seed(3301)
expected_gamma_ppd <- vapply(
  seq_len(2),
  function(j) {
    shape_j <- sigma_ppd[j]^(-2)
    stats::rgamma(
      nrow(x_ppd),
      shape = shape_j,
      rate = shape_j * exp(-linear_predictor_ppd[, j])
    )
  },
  numeric(nrow(x_ppd))
)

stopifnot(
  identical(dim(gamma_ppd), dim(expected_gamma_ppd)),
  isTRUE(all.equal(gamma_ppd, expected_gamma_ppd))
)

# Information-criterion likelihood contributions must use the same gamma
# transformation. With perfect sensitivity, the two examples reduce to an
# interval probability and a right-tail probability, respectively.
ic_mod <- list(
  prev = FALSE,
  v_obs = list(c(0, 1), c(0, 2, Inf)),
  dist = "gamma",
  x_t = matrix(c(-0.5, 0.5), ncol = 1)
)
ic_est <- c(beta, log(sigma), 1)
ic_likelihood <- BayesPIM:::l_obs_2s(
  est = ic_est,
  mod = ic_mod,
  log_scale = FALSE,
  sumup = FALSE
)
ic_par <- BayesPIM:::trans_par_gamma(
  cbind(1, ic_mod$x_t),
  c(beta, log(sigma))
)
expected_ic_likelihood <- c(
  stats::pgamma(1, shape = ic_par[1, "shape"], rate = ic_par[1, "rate"]),
  stats::pgamma(
    2,
    shape = ic_par[2, "shape"],
    rate = ic_par[2, "rate"],
    lower.tail = FALSE
  )
)

stopifnot(
  isTRUE(all.equal(unname(ic_likelihood), expected_ic_likelihood))
)

# Exercise the public information-criterion and summary paths using a minimal
# fitted-object contract. This also verifies that the PSOCK worker receives the
# gamma transformation helper.
make_gamma_chain <- function(offset) {
  draws <- cbind(
    beta_t_intercept = seq(0.15, 0.25, length.out = 20) + offset,
    beta_t_x_1 = seq(-0.14, -0.06, length.out = 20) - offset,
    sigma_t = seq(0.72, 0.88, length.out = 20)
  )
  coda::mcmc(draws)
}

gamma_fit <- structure(
  list(
    par = coda::mcmc.list(
      make_gamma_chain(-0.01),
      make_gamma_chain(0.01)
    ),
    warmup = 0,
    save_every = 1,
    total_iterations = 20,
    update_kappa = FALSE,
    prev = FALSE,
    v_obs = ic_mod$v_obs,
    x_t = ic_mod$x_t,
    x_g = NULL,
    covariate_scaling = list(
      x_t = BayesPIM:::bayespim_covariate_scaling(
        ic_mod$x_t, standardize = FALSE, name = "x_t"
      ),
      x_g = BayesPIM:::bayespim_covariate_scaling(
        NULL, standardize = FALSE, name = "x_g"
      )
    ),
    dist = "gamma",
    sampler = "slice_collapsed",
    kappa = 1,
    fix_sigma = FALSE,
    fix_q = FALSE,
    max_rhat = 1.1,
    convergence = list(max_rhat = 1.1, min_ess = 1)
  ),
  class = "bayespim"
)

set.seed(3302)
gamma_ic <- get_ic(gamma_fit, samples = 10, cores = 1)
stopifnot(
  identical(colnames(gamma_ic), c("WAIC1", "WAIC2", "DIC")),
  all(is.finite(gamma_ic))
)

gamma_summary <- NULL
gamma_summary_output <- capture.output(
  gamma_summary <- summary(gamma_fit)
)
stopifnot(
  identical(gamma_summary$distribution, "gamma"),
  identical(gamma_summary$sampler, "slice_collapsed"),
  any(gamma_summary_output == "Latent-time distribution: gamma"),
  any(gamma_summary_output == "Incidence sampler: slice_collapsed")
)

set.seed(3303)
gamma_public_ppd <- ppCIF(
  gamma_fit,
  pst_samples = 4,
  quant = c(0, 1, 2)
)
stopifnot(
  inherits(gamma_public_ppd, "ppCIF"),
  identical(gamma_public_ppd$quant, c(0, 1, 2)),
  length(gamma_public_ppd$mixture$med_cdf) == 3L,
  all(is.finite(gamma_public_ppd$mixture$med_cdf)),
  all(diff(gamma_public_ppd$mixture$med_cdf) >= 0),
  identical(
    gamma_public_ppd$mixture$med_cdf,
    gamma_public_ppd$nonprevalent$med_cdf
  )
)
