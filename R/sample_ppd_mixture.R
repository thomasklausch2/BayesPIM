#' Draw posterior predictive event times from the prevalence-incidence mixture
#' @noRd
sample_ppd_mixture <- function(
    par_list,
    x_t,
    x_g,
    dist,
    sampled_draws,
    covariate_scaling = NULL
) {
  n <- nrow(x_t)
  times <- sample_ppd_nonprevalent(
    par_list = par_list,
    x_t = x_t,
    dist = dist,
    sampled_draws = sampled_draws,
    covariate_scaling = covariate_scaling$x_t
  )

  par <- as.matrix(par_list[[1L]])
  p1_t <- ncol(x_t) + 1L
  p1_g <- ncol(x_g) + 1L
  incidence_parameter_count <- p1_t + if (dist == "gengamma") 2L else 1L
  beta_g_columns <- seq.int(
    incidence_parameter_count + 1L,
    incidence_parameter_count + p1_g
  )
  beta_g <- par[, beta_g_columns, drop = FALSE]
  if (length(par_list) > 1L) {
    for (i in seq.int(2L, length(par_list))) {
      par <- as.matrix(par_list[[i]])
      beta_g <- rbind(beta_g, par[, beta_g_columns, drop = FALSE])
    }
  }

  linear_predictor <- bayespim_scaled_linear_predictor(
    x = x_g,
    beta = beta_g[sampled_draws, , drop = FALSE],
    scaling = covariate_scaling$x_g
  )
  prob_g <- stats::pnorm(linear_predictor)
  g <- matrix(
    stats::rbinom(length(prob_g), size = 1L, prob = as.vector(prob_g)),
    nrow = n,
    ncol = length(sampled_draws)
  )

  mixture <- times
  mixture[g == 1L] <- 0
  nonprevalent <- times
  nonprevalent[g == 1L] <- NA_real_

  list(
    mixture = as.matrix(mixture),
    nonprevalent = as.matrix(nonprevalent)
  )
}
