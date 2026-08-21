#' Posterior predictive cumulative incidence functions
#'
#' Compute posterior predictive cumulative incidence functions (CIFs) from a
#' fitted \code{bayespim} model.
#'
#' @param mod A fitted model returned by \code{bayespim()}; objects without
#'   class \code{"bayespim"} are rejected.
#' @param fix_x_t Either \code{NULL} for a marginal CIF or a numeric vector with
#'   one entry per incidence-model covariate. Numeric entries fix covariates at
#'   the supplied values; \code{NA} entries are integrated over the empirical
#'   covariate distribution.
#' @param fix_x_g The corresponding vector for prevalence-model covariates.
#'   It cannot be supplied for a model fitted with \code{prev = FALSE}. When it
#'   is \code{NULL}, \code{fix_x_t} is also applied to the prevalence model if
#'   the two design matrices have the same number of columns.
#' @param pst_samples Positive integer giving the number of posterior draws used
#'   for prediction. It cannot exceed the number of stored post-warm-up draws.
#' @param perc Numeric vector of cumulative probabilities in \eqn{[0, 1]} for
#'   which event times are returned when \code{ppd_type = "quantiles"}.
#' @param ppd_type Character string selecting the returned representation.
#'   The default, \code{"percentiles"}, returns cumulative probabilities at the
#'   times in \code{quant}. \code{"quantiles"} returns event times at the
#'   cumulative probabilities in \code{perc}.
#' @param quant Numeric vector of non-negative time points at which cumulative
#'   probabilities are evaluated when \code{ppd_type = "percentiles"}. If
#'   \code{NULL}, a grid from zero to the maximum observed follow-up is used.
#'
#' @details
#' A prevalence-incidence mixture model defines two related CIFs. One
#' \code{ppCIF()} call always computes and stores both; the \code{type}
#' argument of \code{plot.ppCIF()} controls which one is displayed. The
#' \code{mixture} CIF includes prevalent cases as a point mass at time zero.
#' The \code{nonprevalent} CIF describes event times among individuals who are
#' non-prevalent at baseline. Both are computed from the same posterior draws
#' and posterior predictive replicates and retained in the returned object.
#' The candidate draws comprise every stored post-warm-up draw; no additional
#' thinning is applied before the explicit random subsample of
#' \code{pst_samples} draws.
#' Fixed covariate values are supplied on their original scale. For a model
#' fitted with covariate standardization, linear predictors are evaluated on
#' the standardized computational scale to reduce numerical cancellation.
#'
#' For a model fitted with \code{prev = FALSE}, the two CIFs are identical.
#' Posterior prediction supports every incidence distribution available in
#' \code{bayespim()}, including the Prentice generalized gamma. For
#' \code{dist = "gengamma"}, event times are drawn using the fitted location
#' \eqn{\mu = \mathbf{x}_t^\top\boldsymbol{\beta}_t}, positive scale
#' \eqn{\sigma_t}, and signed shape \eqn{Q}.
#'
#' @return An object of class \code{"ppCIF"} with:
#' \describe{
#'   \item{\code{mixture}, \code{nonprevalent}}{Lists containing
#'     \code{med_cdf}, the pointwise posterior predictive median, and
#'     \code{med_cdf_ci}, a two-row matrix with the pointwise 2.5 and 97.5
#'     percent posterior predictive quantiles.}
#'   \item{\code{ppd_type}}{The selected representation.}
#'   \item{\code{quant}, \code{perc}}{The applicable time or probability grid;
#'     the unused grid is \code{NULL}.}
#'   \item{\code{pst_samples}}{The number of posterior draws used.}
#'   \item{\code{distribution}}{The fitted incidence distribution.}
#'   \item{\code{prevalence_model}}{Whether prevalence was fitted.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @seealso \code{\link{plot.ppCIF}}, \code{\link{bayespim}}
#'
#' @examples
#' data(mod)
#' set.seed(2025)
#' cif <- ppCIF(
#'   mod,
#'   pst_samples = 50,
#'   ppd_type = "percentiles",
#'   quant = seq(0, 300, length.out = 51)
#' )
#' plot(cif)
#' plot(cif, type = "nonprevalent")
#' plot(cif, type = "both")
#'
#' @export
ppCIF <- function(
    mod,
    fix_x_t = NULL,
    fix_x_g = NULL,
    pst_samples = 1e3,
    perc = seq(0, 1, 0.01),
    ppd_type = c("percentiles", "quantiles"),
    quant = NULL
) {
  if (!inherits(mod, "bayespim")) {
    stop("`mod` must be an object of class 'bayespim'.", call. = FALSE)
  }

  ppd_type <- match.arg(ppd_type)
  covariate_scaling <- mod$covariate_scaling

  par_list <- bayespim_analysis_chains(mod)
  retained_draws <- sum(vapply(par_list, nrow, integer(1L)))
  if (
    length(pst_samples) != 1L ||
      !is.numeric(pst_samples) ||
      is.na(pst_samples) ||
      !is.finite(pst_samples) ||
      pst_samples <= 0 ||
      pst_samples != as.integer(pst_samples)
  ) {
    stop("`pst_samples` must be a positive integer.", call. = FALSE)
  }
  pst_samples <- as.integer(pst_samples)
  if (pst_samples > retained_draws) {
    stop(
      "`pst_samples` cannot exceed the number of stored post-warm-up posterior draws (",
      retained_draws,
      ").",
      call. = FALSE
    )
  }

  if (
    !is.numeric(perc) ||
      length(perc) == 0L ||
      anyNA(perc) ||
      any(!is.finite(perc)) ||
      any(perc < 0 | perc > 1) ||
      is.unsorted(perc)
  ) {
    stop(
      "`perc` must be a non-decreasing numeric vector with values in [0, 1].",
      call. = FALSE
    )
  }

  v_obs <- mod$v_obs
  n <- length(v_obs)
  x_t <- if (is.null(mod$x_t)) {
    matrix(numeric(0), nrow = n, ncol = 0L)
  } else {
    as.matrix(mod$x_t)
  }
  x_g <- if (is.null(mod$x_g)) {
    matrix(numeric(0), nrow = n, ncol = 0L)
  } else {
    as.matrix(mod$x_g)
  }

  validate_fixed_covariates <- function(x, expected_length, name) {
    if (is.null(x)) {
      return(invisible(NULL))
    }
    if (
      !is.numeric(x) ||
        length(x) != expected_length ||
        any(!is.finite(x[!is.na(x)]))
    ) {
      stop(
        "`",
        name,
        "` must be NULL or a numeric vector of length ",
        expected_length,
        " containing finite values or NA.",
        call. = FALSE
      )
    }
    invisible(NULL)
  }

  validate_fixed_covariates(fix_x_t, ncol(x_t), "fix_x_t")
  validate_fixed_covariates(fix_x_g, ncol(x_g), "fix_x_g")
  if (!isTRUE(mod$prev) && !is.null(fix_x_g)) {
    stop(
      "`fix_x_g` cannot be used for a model fitted with `prev = FALSE`.",
      call. = FALSE
    )
  }

  if (!is.null(fix_x_t)) {
    fixed_t <- !is.na(fix_x_t)
    x_t[, fixed_t] <- fix_x_t[fixed_t]
  }
  if (
    isTRUE(mod$prev) &&
      is.null(fix_x_g) &&
      !is.null(fix_x_t) &&
      ncol(x_g) == ncol(x_t)
  ) {
    fix_x_g <- fix_x_t
  }
  if (!is.null(fix_x_g)) {
    fixed_g <- !is.na(fix_x_g)
    x_g[, fixed_g] <- fix_x_g[fixed_g]
  }

  if (ppd_type == "percentiles") {
    if (is.null(quant)) {
      finite_followup <- vapply(
        v_obs,
        function(visits) {
          visits <- visits[is.finite(visits)]
          if (length(visits) == 0L) NA_real_ else visits[length(visits)]
        },
        numeric(1L)
      )
      if (all(is.na(finite_followup))) {
        stop(
          "A default prediction grid cannot be formed because follow-up times are not finite.",
          call. = FALSE
        )
      }
      max_followup <- max(finite_followup, na.rm = TRUE)
      quant <- if (max_followup == 0) {
        0
      } else {
        seq(0, max_followup, length.out = 1001L)
      }
    }
    if (
      !is.numeric(quant) ||
        length(quant) == 0L ||
        anyNA(quant) ||
        any(!is.finite(quant)) ||
        any(quant < 0) ||
        is.unsorted(quant)
    ) {
      stop(
        "`quant` must be a non-decreasing vector of finite, non-negative times.",
        call. = FALSE
      )
    }
  }

  sampled_draws <- sample.int(
    n = retained_draws,
    size = pst_samples,
    replace = FALSE
  )
  if (isTRUE(mod$prev)) {
    ppd <- sample_ppd_mixture(
      par_list = par_list,
      x_t = x_t,
      x_g = x_g,
      dist = mod$dist,
      sampled_draws = sampled_draws,
      covariate_scaling = covariate_scaling
    )
  } else {
    nonprevalent_ppd <- sample_ppd_nonprevalent(
      par_list = par_list,
      x_t = x_t,
      dist = mod$dist,
      sampled_draws = sampled_draws,
      covariate_scaling = covariate_scaling$x_t
    )
    ppd <- list(
      mixture = nonprevalent_ppd,
      nonprevalent = nonprevalent_ppd
    )
  }

  summarize_ppd <- function(draws, estimand) {
    all_missing <- vapply(
      seq_len(ncol(draws)),
      function(j) all(is.na(draws[, j])),
      logical(1L)
    )
    if (any(all_missing)) {
      draws <- draws[, !all_missing, drop = FALSE]
      warning(
        "Dropped ",
        sum(all_missing),
        " posterior predictive ",
        estimand,
        " draw(s) containing no non-prevalent individuals.",
        call. = FALSE
      )
    }
    if (ncol(draws) == 0L) {
      stop(
        "No posterior predictive ",
        estimand,
        " draws contain a non-prevalent individual.",
        call. = FALSE
      )
    }

    if (ppd_type == "percentiles") {
      draw_values <- vapply(
        seq_len(ncol(draws)),
        function(j) {
          values <- draws[, j]
          values <- values[is.finite(values)]
          stats::ecdf(values)(quant)
        },
        numeric(length(quant))
      )
      draw_values <- matrix(
        draw_values,
        nrow = length(quant),
        ncol = ncol(draws)
      )
    } else {
      draw_values <- vapply(
        seq_len(ncol(draws)),
        function(j) {
          values <- draws[, j]
          values <- values[is.finite(values)]
          stats::quantile(values, probs = perc, names = FALSE)
        },
        numeric(length(perc))
      )
      draw_values <- matrix(
        draw_values,
        nrow = length(perc),
        ncol = ncol(draws)
      )
    }

    list(
      med_cdf = vapply(
        seq_len(nrow(draw_values)),
        function(i) stats::median(draw_values[i, ]),
        numeric(1L)
      ),
      med_cdf_ci = vapply(
        seq_len(nrow(draw_values)),
        function(i) {
          stats::quantile(
            draw_values[i, ],
            probs = c(0.025, 0.975),
            names = FALSE
          )
        },
        numeric(2L)
      )
    )
  }

  ret <- list(
    mixture = summarize_ppd(ppd$mixture, "mixture"),
    nonprevalent = summarize_ppd(ppd$nonprevalent, "non-prevalent"),
    ppd_type = ppd_type,
    quant = if (ppd_type == "percentiles") quant else NULL,
    perc = if (ppd_type == "quantiles") perc else NULL,
    pst_samples = pst_samples,
    distribution = mod$dist,
    prevalence_model = isTRUE(mod$prev),
    call = match.call()
  )
  class(ret) <- c("ppCIF", "list")
  ret
}
