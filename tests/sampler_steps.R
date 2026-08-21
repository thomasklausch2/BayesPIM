library(BayesPIM)

x <- cbind(1, c(-0.5, 0, 0.5, 1))
times <- c(1.5, 2.0, 2.5, 3.0)
eta <- c(log(2), 0.1, log(0.8))
data <- list(times = times, rows = seq_along(times))
q_prior_sd <- 0.75

strict_prior <- function(
    eta, dist, beta_prior, tau_t, sig_prior, q_prior_sd
) {
  stopifnot(
    dist %in% c("weibull", "gengamma"),
    identical(beta_prior, "norm"),
    isTRUE(all.equal(tau_t, 1)),
    isTRUE(all.equal(sig_prior, 1)),
    isTRUE(all.equal(q_prior_sd, 0.75))
  )
  log_aft_prior(
    eta = eta,
    dist = dist,
    beta_prior = beta_prior,
    tau_t = tau_t,
    sig_prior = sig_prior,
    q_prior_sd = q_prior_sd
  )
}

# Both log-likelihoods must reject an unsupported family by name, rather than
# failing downstream on a length mismatch or an unassigned result. The exact-time
# likelihood does not cover the generalized gamma; the interval-censored one does.
unsupported_message <- function(expr) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"))
  conditionMessage(error)
}
stopifnot(
  identical(
    unsupported_message(
      BayesPIM:::ll_aft(eta = eta, y = times, x = x, dist = "nonesuch")
    ),
    "Unsupported distribution: nonesuch"
  ),
  identical(
    unsupported_message(
      BayesPIM:::ll_aft_ic(
        eta = eta, L = c(1, 1.5, 2, 2.5), R = c(2, 2.5, 3, Inf),
        x = x, dist = "nonesuch"
      )
    ),
    "Unsupported distribution: nonesuch"
  ),
  identical(
    unsupported_message(
      BayesPIM:::ll_aft(eta = c(eta, 0.4), y = times, x = x, dist = "gengamma")
    ),
    "Unsupported distribution: gengamma"
  )
)
# The interval-censored likelihood does support the generalized gamma.
stopifnot(is.finite(BayesPIM:::ll_aft_ic(
  eta = c(eta, 0.4), L = c(1, 1.5, 2, 2.5), R = c(2, 2.5, 3, Inf),
  x = x, dist = "gengamma"
)))

set.seed(11)
mh <- BayesPIM:::step_mh(
  eta = eta,
  data = data,
  x = x,
  dist = "weibull",
  log_ll = BayesPIM:::ll_aft,
  log_prior_fun = strict_prior,
  beta_prior = "norm",
  tau_t = 1,
  sig_prior = 1,
  q_prior_sd = q_prior_sd,
  control = list(prop_var = diag(0.01, length(eta)))
)
stopifnot(identical(names(mh), c("eta", "log_post", "accepted")))
stopifnot(length(mh$eta) == length(eta), is.finite(mh$log_post))

set.seed(12)
slice <- BayesPIM:::step_slice(
  eta = eta,
  data = data,
  x = x,
  dist = "weibull",
  log_ll = BayesPIM:::ll_aft,
  log_prior_fun = strict_prior,
  beta_prior = "norm",
  tau_t = 1,
  sig_prior = 1,
  q_prior_sd = q_prior_sd,
  control = list(width = 0.5)
)
stopifnot(identical(names(slice), c("eta", "log_post", "accepted")))
stopifnot(length(slice$eta) == length(eta), is.finite(slice$log_post))

interval_data <- list(
  L = c(1, 1.5, 2, 2.5),
  R = c(2, 2.5, 3, Inf),
  rows = seq_len(nrow(x))
)
set.seed(13)
collapsed <- BayesPIM:::step_slice_collapsed(
  eta = eta,
  data = interval_data,
  x = x,
  dist = "weibull",
  log_ll = BayesPIM:::ll_aft_ic,
  log_prior_fun = strict_prior,
  beta_prior = "norm",
  tau_t = 1,
  sig_prior = 1,
  q_prior_sd = q_prior_sd,
  control = list(width = 0.5, fix_q = FALSE)
)
stopifnot(identical(names(collapsed), c("eta", "log_post", "accepted")))
stopifnot(length(collapsed$eta) == length(eta), is.finite(collapsed$log_post))

eta_gengamma <- c(eta, 0.4)

set.seed(14)
collapsed_gengamma <- BayesPIM:::step_slice_collapsed(
  eta = eta_gengamma,
  data = interval_data,
  x = x,
  dist = "gengamma",
  log_ll = BayesPIM:::ll_aft_ic,
  log_prior_fun = strict_prior,
  beta_prior = "norm",
  tau_t = 1,
  sig_prior = 1,
  q_prior_sd = q_prior_sd,
  control = list(width = 0.5, fix_q = FALSE)
)
stopifnot(
  identical(names(collapsed_gengamma), c("eta", "log_post", "accepted")),
  length(collapsed_gengamma$eta) == length(eta_gengamma),
  is.finite(collapsed_gengamma$log_post)
)

eta_gengamma_fixed <- eta_gengamma
eta_gengamma_fixed[length(eta_gengamma_fixed)] <- q_prior_sd

set.seed(15)
collapsed_gengamma_fixed <- BayesPIM:::step_slice_collapsed(
  eta = eta_gengamma_fixed,
  data = interval_data,
  x = x,
  dist = "gengamma",
  log_ll = BayesPIM:::ll_aft_ic,
  log_prior_fun = strict_prior,
  beta_prior = "norm",
  tau_t = 1,
  sig_prior = 1,
  q_prior_sd = q_prior_sd,
  control = list(width = 0.5, fix_q = TRUE)
)
stopifnot(
  identical(
    collapsed_gengamma_fixed$eta[length(eta_gengamma_fixed)],
    q_prior_sd
  ),
  is.finite(collapsed_gengamma_fixed$log_post)
)

# The distribution wrappers must reject an unrecognised family by name too,
# rather than falling through and returning NULL, which would propagate as
# numeric(0) and silently sum to zero in a likelihood.
for (wrapper in c("pdist", "qdist", "rdist", "ddist")) {
  stopifnot(identical(
    unsupported_message(
      get(wrapper, envir = asNamespace("BayesPIM"))(0.5, cbind(1, 1, 0.5), "nonesuch")
    ),
    "Unsupported distribution: nonesuch"
  ))
}
# Every supported family must still round-trip and match its reference density.
wrapper_par <- list(
  exp = cbind(0.7), gamma = cbind(2, 1.5), weibull = cbind(1.3, 2),
  loglog = cbind(1.1, 2.2), lognormal = cbind(0.2, 0.6),
  gengamma = cbind(0.2, 0.6, 0.4)
)
wrapper_x <- c(0.3, 1, 2.5)
for (family in names(wrapper_par)) {
  pars <- wrapper_par[[family]][rep(1L, length(wrapper_x)), , drop = FALSE]
  probs <- BayesPIM:::pdist(wrapper_x, pars, family)
  stopifnot(
    all(is.finite(probs)), all(probs >= 0 & probs <= 1),
    isTRUE(all.equal(BayesPIM:::qdist(probs, pars, family), wrapper_x)),
    all(is.finite(BayesPIM:::ddist(wrapper_x, pars, family))),
    length(BayesPIM:::rdist(length(wrapper_x), pars, family)) == length(wrapper_x)
  )
}

# The complete prior contract must also be honored during fresh-run
# initialization, before the sampler-specific update steps begin.
set.seed(18)
initialization_data <- gen_data(
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
initialized_fit <- bayespim(
  v_obs = initialization_data$v_obs,
  x_t = initialization_data$x,
  r = initialization_data$r,
  kappa = 0.8,
  ndraws = 4,
  warmup = 1,
  chains = 1,
  sampler = "mh",
  prop_sd = 0.01,
  log_prior_fun = strict_prior,
  q_prior_sd = q_prior_sd,
  prev = FALSE,
  silent = TRUE
)
stopifnot(inherits(initialized_fit, "bayespim"))

# `search_prop_sd()` must report its own calibration trace but not the progress
# output of the `bayespim()` updates it drives. `capture.output()` cannot do this
# because that progress is emitted with message(), on stderr; `silent = TRUE` is
# the supported switch and leaves warnings and errors intact.
#
# The acceptance check is strict (`ac > lower & ac < upper`), so bounds must lie
# strictly outside [0, 1] for every possible rate to count as a success. With
# acc_bounds = c(-1, 2) each pass succeeds, so succ_min = 2 gives exactly two
# iterations: the first reads the supplied fit, the second performs one update.
# That keeps the otherwise unbounded search loop deterministic and short.
search_messages <- character()
search_result <- withCallingHandlers(
  search_prop_sd(
    m = initialized_fit,
    ndraws = 4,
    succ_min = 2,
    acc_bounds = c(-1, 2)
  ),
  message = function(m) {
    search_messages <<- c(search_messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)
search_messages <- paste(search_messages, collapse = "")
stopifnot(
  identical(names(search_result), c("prop_sd", "ac")),
  is.finite(search_result$prop_sd),
  search_result$prop_sd > 0,
  is.finite(search_result$ac),
  # its own trace is kept
  grepl("Iteration 1", search_messages, fixed = TRUE),
  grepl("Iteration 2", search_messages, fixed = TRUE),
  !grepl("Iteration 3", search_messages, fixed = TRUE),
  grepl("Acceptance rate was:", search_messages, fixed = TRUE),
  grepl("Finished calibrating proposal variance.", search_messages, fixed = TRUE),
  # the driven bayespim() updates stay quiet
  !grepl("Updating previous MCMC run", search_messages, fixed = TRUE),
  !grepl("Convergence diagnostics", search_messages, fixed = TRUE),
  !grepl("Convergence criteria", search_messages, fixed = TRUE)
)
