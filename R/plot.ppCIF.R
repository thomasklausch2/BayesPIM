#' Plot posterior predictive cumulative incidence functions
#'
#' Plot the mixture CIF, the non-prevalent CIF, or both CIFs from an object
#' returned by \code{\link{ppCIF}}.
#'
#' The plotting contract is independent of the fitted incidence distribution
#' and therefore also applies to generalized-gamma posterior predictions.
#'
#' @param x An object of class \code{"ppCIF"}.
#' @param y Ignored.
#' @param type Character string selecting the CIF to draw:
#'   \code{"mixture"} (the default), \code{"nonprevalent"}, or \code{"both"}.
#'   The last choice produces horizontally adjacent mixture and non-prevalent
#'   panels on a common time scale.
#' @param ci Logical. If \code{TRUE}, show pointwise 95 percent credible
#'   regions.
#' @param main Optional plot title. For \code{type = "both"}, this may be a
#'   character vector of length two.
#' @param xlab,ylab Axis labels.
#' @param xlim,ylim Optional increasing numeric vectors of length two giving
#'   axis limits. A common \code{xlim} is used for both panels when
#'   \code{type = "both"}. The default probability limits are \code{c(0, 1)}.
#' @param col Line color.
#' @param ci_col Credible-region fill color. If \code{NULL}, a transparent
#'   version of \code{col} is used.
#' @param lwd Line width.
#' @param ... Further graphical arguments passed to \code{\link[graphics]{plot.default}}.
#'
#' @details
#' For \code{ppd_type = "percentiles"}, time is shown on the x-axis and
#' uncertainty is represented vertically. For \code{ppd_type = "quantiles"},
#' the stored inverse-CDF representation is drawn with uncertainty in time.
#' Pointwise regions can be omitted with \code{ci = FALSE}.
#'
#' @return Invisibly returns the supplied \code{"ppCIF"} object. The method is
#'   called for its base-graphics plotting side effect.
#'
#' @examples
#' data(mod)
#' set.seed(2025)
#' cif <- ppCIF(
#'   mod,
#'   pst_samples = 50,
#'   quant = seq(0, 300, length.out = 51)
#' )
#' plot(cif, type = "both")
#'
#' @method plot ppCIF
#' @export
plot.ppCIF <- function(
    x,
    y = NULL,
    type = c("mixture", "nonprevalent", "both"),
    ci = TRUE,
    main = NULL,
    xlab = "Time",
    ylab = "Cumulative incidence",
    xlim = NULL,
    ylim = c(0, 1),
    col = "#0072B2",
    ci_col = NULL,
    lwd = 2,
    ...
) {
  if (!inherits(x, "ppCIF")) {
    stop("`x` must be an object of class 'ppCIF'.", call. = FALSE)
  }
  type <- match.arg(type)
  if (
    is.null(x$ppd_type) ||
      !x$ppd_type %in% c("percentiles", "quantiles") ||
      is.null(x$mixture) ||
      is.null(x$nonprevalent)
  ) {
    stop("`x` is not a valid object returned by `ppCIF()`.", call. = FALSE)
  }
  if (
    !is.logical(ci) ||
      length(ci) != 1L ||
      is.na(ci)
  ) {
    stop("`ci` must be TRUE or FALSE.", call. = FALSE)
  }
  validate_limits <- function(limits, name) {
    if (
      !is.null(limits) &&
        (
          !is.numeric(limits) ||
            length(limits) != 2L ||
            anyNA(limits) ||
            any(!is.finite(limits)) ||
            limits[1L] >= limits[2L]
        )
    ) {
      stop(
        "`",
        name,
        "` must be NULL or an increasing finite numeric vector of length two.",
        call. = FALSE
      )
    }
  }
  validate_limits(xlim, "xlim")
  validate_limits(ylim, "ylim")
  if (!is.null(main) && !length(main) %in% c(1L, 2L)) {
    stop("`main` must be NULL or a character vector of length one or two.", call. = FALSE)
  }
  if (is.null(ci_col)) {
    ci_col <- grDevices::adjustcolor(col, alpha.f = 0.2)
  }

  curve_coordinates <- function(curve) {
    if (
      is.null(curve$med_cdf) ||
        is.null(curve$med_cdf_ci) ||
        nrow(curve$med_cdf_ci) != 2L
    ) {
      stop("`x` contains an invalid CIF component.", call. = FALSE)
    }
    if (x$ppd_type == "percentiles") {
      list(
        time = x$quant,
        probability = curve$med_cdf,
        lower_time = NULL,
        upper_time = NULL,
        lower_probability = curve$med_cdf_ci[1L, ],
        upper_probability = curve$med_cdf_ci[2L, ]
      )
    } else {
      list(
        time = curve$med_cdf,
        probability = x$perc,
        lower_time = curve$med_cdf_ci[1L, ],
        upper_time = curve$med_cdf_ci[2L, ],
        lower_probability = NULL,
        upper_probability = NULL
      )
    }
  }

  selected <- if (type == "both") {
    c("mixture", "nonprevalent")
  } else {
    type
  }
  coordinates <- lapply(selected, function(name) curve_coordinates(x[[name]]))
  names(coordinates) <- selected

  if (is.null(xlim)) {
    time_values <- unlist(
      lapply(
        coordinates,
        function(coords) {
          c(coords$time, coords$lower_time, coords$upper_time)
        }
      ),
      use.names = FALSE
    )
    time_values <- time_values[is.finite(time_values)]
    if (length(time_values) == 0L) {
      stop("The selected CIF contains no finite time values.", call. = FALSE)
    }
    xlim <- range(time_values)
    if (xlim[1L] == xlim[2L]) {
      padding <- if (xlim[1L] == 0) 0.5 else abs(xlim[1L]) * 0.04
      xlim <- xlim + c(-padding, padding)
    }
  }

  default_titles <- c(
    mixture = "Posterior predictive mixture CIF",
    nonprevalent = "Posterior predictive non-prevalent CIF"
  )
  if (is.null(main)) {
    panel_titles <- unname(default_titles[selected])
  } else if (length(main) == 1L) {
    panel_titles <- rep(main, length(selected))
  } else {
    panel_titles <- main
  }

  draw_curve <- function(coords, title) {
    graphics::plot(
      coords$time,
      coords$probability,
      type = "n",
      xlim = xlim,
      ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      main = title,
      ...
    )
    if (ci) {
      if (x$ppd_type == "percentiles") {
        graphics::polygon(
          x = c(coords$time, rev(coords$time)),
          y = c(
            coords$lower_probability,
            rev(coords$upper_probability)
          ),
          border = NA,
          col = ci_col
        )
      } else {
        graphics::polygon(
          x = c(coords$lower_time, rev(coords$upper_time)),
          y = c(coords$probability, rev(coords$probability)),
          border = NA,
          col = ci_col
        )
      }
    }
    graphics::lines(
      coords$time,
      coords$probability,
      col = col,
      lwd = lwd
    )
  }

  if (type == "both") {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(1L, 2L))
  }
  for (i in seq_along(coordinates)) {
    draw_curve(coordinates[[i]], panel_titles[i])
  }

  invisible(x)
}
