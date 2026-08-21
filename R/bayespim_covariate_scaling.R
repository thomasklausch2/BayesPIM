#' Create the scaling record for one BayesPIM design matrix
#' @noRd
bayespim_covariate_scaling <- function(x, standardize, name) {
  if (is.null(x)) {
    return(list(
      center = numeric(0L),
      scale = numeric(0L),
      binary = logical(0L),
      standardized = logical(0L)
    ))
  }

  x <- as.matrix(x)
  binary <- apply(x, 2L, function(column) {
    length(unique(column)) == 2L
  })
  zero_variance <- apply(x, 2L, function(column) {
    length(unique(column)) < 2L
  })
  if (any(zero_variance)) {
    stop(
      sprintf(
        "`%s` contains zero-variance column(s): %s.",
        name,
        paste(colnames(x)[zero_variance], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  standardized <- rep(isTRUE(standardize), ncol(x)) & !binary
  center <- rep(0, ncol(x))
  scale <- rep(1, ncol(x))
  center[standardized] <- colMeans(x[, standardized, drop = FALSE])
  scale[standardized] <- apply(
    x[, standardized, drop = FALSE],
    2L,
    stats::sd
  )
  if (any(!is.finite(scale) | scale <= 0)) {
    stop(
      sprintf("`%s` contains a column that cannot be standardized.", name),
      call. = FALSE
    )
  }

  names(center) <- colnames(x)
  names(scale) <- colnames(x)
  names(binary) <- colnames(x)
  names(standardized) <- colnames(x)
  list(
    center = center,
    scale = scale,
    binary = binary,
    standardized = standardized
  )
}

#' Apply a stored BayesPIM covariate transformation
#' @noRd
bayespim_standardize_matrix <- function(x, scaling, name) {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  if (
    is.null(scaling) ||
      length(scaling$center) != ncol(x) ||
      length(scaling$scale) != ncol(x)
  ) {
    stop(
      sprintf("Stored covariate-scaling metadata do not match `%s`.", name),
      call. = FALSE
    )
  }
  sweep(sweep(x, 2L, scaling$center, "-"), 2L, scaling$scale, "/")
}
