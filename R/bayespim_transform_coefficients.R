#' Transform BayesPIM regression coefficients between covariate scales
#' @noRd
bayespim_transform_coefficients <- function(
    par,
    p1_t,
    p1_g,
    prev,
    has_q,
    covariate_scaling,
    direction = c("internal_to_public", "public_to_internal")
) {
  direction <- match.arg(direction)
  par <- as.matrix(par)
  if (is.null(covariate_scaling)) return(par)

  transform_block <- function(par, indices, scaling) {
    if (length(indices) <= 1L) return(par)
    slope_indices <- indices[-1L]
    slopes <- par[, slope_indices, drop = FALSE]
    if (direction == "internal_to_public") {
      slopes <- sweep(slopes, 2L, scaling$scale, "/")
      par[, indices[1L]] <- par[, indices[1L]] -
        drop(slopes %*% scaling$center)
      par[, slope_indices] <- slopes
    } else {
      par[, indices[1L]] <- par[, indices[1L]] +
        drop(slopes %*% scaling$center)
      par[, slope_indices] <- sweep(slopes, 2L, scaling$scale, "*")
    }
    par
  }

  par <- transform_block(
    par,
    indices = seq_len(p1_t),
    scaling = covariate_scaling$x_t
  )
  if (isTRUE(prev)) {
    first_g <- p1_t + 2L + as.integer(has_q)
    par <- transform_block(
      par,
      indices = seq.int(first_g, length.out = p1_g),
      scaling = covariate_scaling$x_g
    )
  }
  par
}

#' Form a stable linear predictor using public BayesPIM coefficients
#' @noRd
bayespim_scaled_linear_predictor <- function(x, beta, scaling) {
  x <- as.matrix(x)
  beta <- as.matrix(beta)
  if (is.null(scaling) || ncol(x) == 0L) {
    return(cbind(1, x) %*% t(beta))
  }
  x_scaled <- bayespim_standardize_matrix(x, scaling, "prediction covariates")
  slopes <- beta[, -1L, drop = FALSE]
  centered_intercept <- beta[, 1L] + drop(slopes %*% scaling$center)
  scaled_slopes <- sweep(slopes, 2L, scaling$scale, "*")
  cbind(1, x_scaled) %*% t(cbind(centered_intercept, scaled_slopes))
}
