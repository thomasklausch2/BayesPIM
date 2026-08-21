#' Draw posterior predictive event times from the incidence model
#' @noRd
sample_ppd_nonprevalent <- function(
    par_list,
    x_t,
    dist,
    sampled_draws,
    covariate_scaling = NULL
) {
  par <- as.matrix(par_list[[1L]])
  p1_t <- ncol(x_t) + 1L
  beta_t <- par[, seq_len(p1_t), drop = FALSE]
  sigma <- par[, p1_t + 1L]
  q <- if (dist == "gengamma") par[, p1_t + 2L] else NULL
  if (length(par_list) > 1L) {
    for (i in seq.int(2L, length(par_list))) {
      par <- as.matrix(par_list[[i]])
      beta_t <- rbind(
        beta_t,
        par[, seq_len(p1_t), drop = FALSE]
      )
      sigma <- c(sigma, par[, p1_t + 1L])
      if (dist == "gengamma") {
        q <- c(q, par[, p1_t + 2L])
      }
    }
  }

  linear_predictor <- bayespim_scaled_linear_predictor(
    x = x_t,
    beta = beta_t[sampled_draws, , drop = FALSE],
    scaling = covariate_scaling
  )
  sigma_draws <- sigma[sampled_draws]

  if (dist == "gengamma") {
    q_draws <- q[sampled_draws]
    times <- vapply(
      seq_along(sampled_draws),
      function(j) {
        rdist(
          n = nrow(linear_predictor),
          par = cbind(
            mu = linear_predictor[, j],
            sigma = rep(sigma_draws[j], nrow(linear_predictor)),
            Q = rep(q_draws[j], nrow(linear_predictor))
          ),
          dist = dist
        )
      },
      numeric(nrow(linear_predictor))
    )
    return(as.matrix(times))
  }

  if (dist == "gamma") {
    times <- vapply(
      seq_along(sampled_draws),
      function(j) {
        shape <- sigma_draws[j]^-2
        rate <- shape * exp(-linear_predictor[, j])
        rdist(
          n = nrow(linear_predictor),
          par = cbind(shape = rep(shape, nrow(linear_predictor)), rate = rate),
          dist = dist
        )
      },
      numeric(nrow(linear_predictor))
    )
    return(as.matrix(times))
  }

  if (dist != "lognormal") {
    a <- sigma_draws^-1
    b <- exp(linear_predictor)
  }
  if (dist == "lognormal") {
    a <- sigma_draws
    b <- linear_predictor
  }
  times <- apply(
    rbind(a, b),
    2,
    function(x) {
      rdist(
        n = length(x) - 1L,
        par = cbind(x[2:length(x)], x[1]),
        dist = dist
      )
    }
  )
  as.matrix(times)
}
