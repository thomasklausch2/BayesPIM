library(BayesPIM)

# `cor2cov()` must scale a correlation matrix by its standard deviations for any
# number of variables. `diag()` reads a length-one argument as a matrix
# dimension, so the single-covariate case is checked explicitly.
set.seed(5150)
for (p in 2:4) {
  correlation <- matrix(0.3, p, p)
  diag(correlation) <- 1
  standard_deviations <- runif(p, 0.2, 3)
  stopifnot(isTRUE(all.equal(
    BayesPIM:::cor2cov(correlation, standard_deviations),
    diag(standard_deviations) %*% correlation %*% diag(standard_deviations)
  )))
}
stopifnot(
  isTRUE(all.equal(
    BayesPIM:::cor2cov(matrix(1, 1, 1), 2),
    matrix(4, 1, 1)
  )),
  isTRUE(all.equal(
    BayesPIM:::cor2cov(matrix(1, 1, 1), 0.5),
    matrix(0.25, 1, 1)
  )),
  identical(dim(BayesPIM:::cor2cov(matrix(0, 0, 0), numeric(0))), c(0L, 0L))
)

length_mismatch <- try(
  BayesPIM:::cor2cov(matrix(1, 2, 2), c(1, 2, 3)),
  silent = TRUE
)
stopifnot(
  inherits(length_mismatch, "try-error"),
  grepl(
    "must supply one standard deviation per variable: 3 supplied, 2 expected",
    as.character(length_mismatch),
    fixed = TRUE
  )
)

simulate <- function(p, p_discrete = 0, s = 1, rho = 0, n = 4000) {
  k <- p + p_discrete
  gen_data(
    n = n, p = p, p_discrete = p_discrete, s = s, rho = rho,
    beta_t = rep(0.1, k), beta_g = rep(0.1, k),
    theta = 0.15, mu_t = 3, sigma_t = 0.3,
    v_min = 20, v_max = 30, mean_rc = 80, prob_r = 1
  )
}

# A single continuous covariate with a non-unit standard deviation.
set.seed(5151)
one_covariate <- simulate(p = 1, s = 2)
stopifnot(
  identical(dim(one_covariate$x), c(4000L, 1L)),
  identical(colnames(one_covariate$x), "x_1"),
  abs(stats::sd(one_covariate$x[, 1]) - 2) < 0.15
)

set.seed(5152)
small_sd <- simulate(p = 1, p_discrete = 1, s = 0.5)
stopifnot(
  identical(dim(small_sd$x), c(4000L, 2L)),
  abs(stats::sd(small_sd$x[, 1]) - 0.5) < 0.05,
  all(small_sd$x[, 2] %in% c(0, 1))
)

# One standard deviation per continuous covariate, plus the requested
# correlation between them.
set.seed(5153)
vector_sd <- simulate(p = 3, s = c(0.5, 1, 2), rho = 0.4)
stopifnot(
  identical(dim(vector_sd$x), c(4000L, 3L)),
  max(abs(apply(vector_sd$x, 2, stats::sd) - c(0.5, 1, 2))) < 0.1,
  max(abs(stats::cor(vector_sd$x)[lower.tri(diag(3))] - 0.4)) < 0.05
)

invalid_sd <- list(0, -1, c(1, 2), c(1, NA), Inf)
for (bad in invalid_sd) {
  bad_fit <- try(simulate(p = 3, s = bad, n = 20), silent = TRUE)
  stopifnot(
    inherits(bad_fit, "try-error"),
    grepl(
      "`s` must be a single positive number or one positive number per continuous covariate.",
      as.character(bad_fit),
      fixed = TRUE
    )
  )
}

# `rho` is the correlation between the continuous covariates; the baseline-test
# indicator is a separate quantity returned as `$r`. Equicorrelation is positive definite
# only strictly inside (-1/(p - 1), 1).
set.seed(5156)
negative_correlation <- simulate(p = 3, rho = -0.45)
stopifnot(
  max(abs(stats::cor(negative_correlation$x)[lower.tri(diag(3))] + 0.45)) < 0.05,
  all(negative_correlation$r %in% c(0, 1)),
  length(negative_correlation$r) == 4000L
)

for (p_covariates in c(2L, 3L, 5L)) {
  lower_bound <- -1 / (p_covariates - 1)
  for (bad_rho in c(lower_bound, lower_bound - 0.05, 1, 1.2)) {
    bad_fit <- try(
      simulate(p = p_covariates, rho = bad_rho, n = 20),
      silent = TRUE
    )
    stopifnot(
      inherits(bad_fit, "try-error"),
      grepl(
        sprintf(
          "`rho` must be a single correlation in (%s, 1) for p = %d continuous covariates",
          format(lower_bound, digits = 4), p_covariates
        ),
        as.character(bad_fit),
        fixed = TRUE
      )
    )
  }
  # Just inside the bound remains admissible.
  stopifnot(!inherits(
    try(simulate(p = p_covariates, rho = lower_bound + 0.02, n = 20), silent = TRUE),
    "try-error"
  ))
}

non_scalar_rho <- try(simulate(p = 3, rho = c(0.2, 0.5), n = 20), silent = TRUE)
stopifnot(
  inherits(non_scalar_rho, "try-error"),
  grepl(
    "`rho` must be a single finite correlation between the continuous covariates.",
    as.character(non_scalar_rho),
    fixed = TRUE
  )
)

# With at most one continuous covariate there is nothing to correlate, so `rho`
# has no effect and is not constrained.
set.seed(5157)
single_covariate <- simulate(p = 1, rho = -5, n = 20)
no_continuous_covariate <- simulate(p = 0, p_discrete = 1, rho = 99, n = 20)
stopifnot(
  identical(dim(single_covariate$x), c(20L, 1L)),
  identical(dim(no_continuous_covariate$x), c(20L, 1L)),
  all(no_continuous_covariate$x[, 1] %in% c(0, 1))
)

# `prob_g` must be the prevalence probability actually used to draw `g`, i.e. the
# link of that branch's own linear predictor. It was previously reconstructed on
# return with probit constants regardless of `sel_mod`, which was wrong for
# `sel_mod = "logit"` in both the link and the intercept.
for (link in c("probit", "logit")) {
  for (theta_value in c(0.1, 0.25)) {
    set.seed(6201)
    linked <- gen_data(
      n = 500, p = 1, beta_t = 0.2, beta_g = 0.8, theta = theta_value,
      mu_t = 2, sigma_t = 0.5, v_min = 1, v_max = 3, mean_rc = 8,
      prob_r = 1, sel_mod = link
    )
    linear_predictor <- drop(cbind(1, linked$x) %*% c(
      if (link == "probit") stats::qnorm(theta_value) else log(theta_value / (1 - theta_value)),
      0.8
    ))
    expected <- if (link == "probit") {
      stats::pnorm(linear_predictor)
    } else {
      1 / (1 + exp(-linear_predictor))
    }
    stopifnot(
      isTRUE(all.equal(drop(linked$prob_g), expected, tolerance = 0)),
      all(linked$g %in% c(0, 1)),
      length(linked$prob_g) == 500L
    )
  }
}

# With no prevalence slope, `theta` is the prevalence probability under both links.
set.seed(6202)
for (link in c("probit", "logit")) {
  flat <- gen_data(
    n = 20000, p = 1, beta_t = 0.2, beta_g = 0, theta = 0.25,
    mu_t = 2, sigma_t = 0.5, v_min = 1, v_max = 3, mean_rc = 8,
    prob_r = 1, sel_mod = link
  )
  stopifnot(
    isTRUE(all.equal(unique(round(drop(flat$prob_g), 12)), 0.25)),
    abs(mean(flat$g) - 0.25) < 0.02
  )
}

# The "every remaining test missed" probability is the geometric survival term
# (1 - kappa)^N. Computing it as 1 - sum(p) loses all precision for long
# screening series and eventually goes negative, which rmultinom() rejects. The
# closed form must be exact and non-negative across the whole range.
for (kappa_value in c(0.6, 0.8, 0.95, 0.99)) {
  for (n_tests in c(1L, 5L, 20L, 30L, 60L, 200L)) {
    closed_form <- (1 - kappa_value)^n_tests
    stopifnot(
      closed_form >= 0,
      isTRUE(all.equal(
        sum(c(kappa_value * (1 - kappa_value)^(seq_len(n_tests) - 1L), closed_form)),
        1
      ))
    )
  }
}

# gen_data must therefore survive high sensitivity combined with long screening
# series, which previously failed with "negative probability".
for (settings in list(
  list(kappa = 0.80, mean_rc = 8,   prob_r = 1),
  list(kappa = 0.90, mean_rc = 60,  prob_r = 1),
  list(kappa = 0.95, mean_rc = 100, prob_r = 0.5),
  list(kappa = 0.99, mean_rc = 200, prob_r = 0.5)
)) {
  set.seed(6301)
  long_series <- gen_data(
    n = 300, p = 1, beta_t = 0.2, beta_g = 0.2, theta = 0.2,
    mu_t = 2, sigma_t = 0.5, v_min = 1, v_max = 2,
    kappa = settings$kappa, mean_rc = settings$mean_rc, prob_r = settings$prob_r
  )
  stopifnot(
    length(long_series$v_obs) == 300L,
    !any(vapply(long_series$v_obs, anyNA, logical(1))),
    all(vapply(long_series$v_obs, function(z) z[1] == 0, logical(1))),
    all(vapply(long_series$v_obs, function(z) !is.unsorted(z), logical(1)))
  )
}

# An intercept-only design must still produce one row per individual, without
# relying on array recycling.
set.seed(5154)
warnings_seen <- character()
intercept_only <- withCallingHandlers(
  gen_data(
    n = 50, p = 0, p_discrete = 0,
    theta = 0.15, mu_t = 3, sigma_t = 0.3,
    v_min = 20, v_max = 30, mean_rc = 80, prob_r = 1
  ),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  length(warnings_seen) == 0L,
  is.null(intercept_only$x),
  length(intercept_only$times_true) == 50L,
  length(intercept_only$g) == 50L,
  length(intercept_only$r) == 50L,
  length(intercept_only$prob_g) == 50L,
  length(intercept_only$v_obs) == 50L
)

# A single drawn row must keep its matrix shape.
set.seed(5155)
single_row <- simulate(p = 2, s = c(1, 3), n = 1)
stopifnot(identical(dim(single_row$x), c(1L, 2L)))
