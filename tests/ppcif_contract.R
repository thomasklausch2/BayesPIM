library(BayesPIM)

make_ppcif_chain <- function(offset = 0) {
  coda::mcmc(cbind(
    beta_t_intercept = seq(0.15, 0.25, length.out = 20) + offset,
    beta_t_x_1 = seq(-0.14, -0.06, length.out = 20),
    sigma_t = seq(0.72, 0.88, length.out = 20),
    beta_g_intercept = seq(-0.1, 0.1, length.out = 20),
    beta_g_x_1 = seq(0.1, 0.2, length.out = 20)
  ))
}

ppcif_fit <- structure(
  list(
    par = coda::mcmc.list(
      make_ppcif_chain(-0.01),
      make_ppcif_chain(0.01)
    ),
    warmup = 0,
    save_every = 1,
    total_iterations = 20,
    prev = TRUE,
    v_obs = list(c(0, 1), c(0, 2, Inf), c(0, 3), c(0, 4, Inf)),
    x_t = matrix(c(-1, -0.5, 0.5, 1), ncol = 1),
    x_g = matrix(c(-1, -0.5, 0.5, 1), ncol = 1),
    covariate_scaling = list(
      x_t = BayesPIM:::bayespim_covariate_scaling(
        matrix(c(-1, -0.5, 0.5, 1), ncol = 1), standardize = FALSE, name = "x_t"
      ),
      x_g = BayesPIM:::bayespim_covariate_scaling(
        matrix(c(-1, -0.5, 0.5, 1), ncol = 1), standardize = FALSE, name = "x_g"
      )
    ),
    dist = "gamma"
  ),
  class = "bayespim"
)

not_a_fit <- try(ppCIF(list()), silent = TRUE)
stopifnot(
  inherits(not_a_fit, "try-error"),
  grepl(
    "`mod` must be an object of class 'bayespim'.",
    as.character(not_a_fit),
    fixed = TRUE
  )
)

set.seed(4401)
ppcif_percentiles <- ppCIF(
  ppcif_fit,
  pst_samples = 20,
  quant = c(0, 0.5, 1, 2)
)
stopifnot(
  inherits(ppcif_percentiles, "ppCIF"),
  identical(
    names(ppcif_percentiles)[1:2],
    c("mixture", "nonprevalent")
  ),
  identical(ppcif_percentiles$ppd_type, "percentiles"),
  identical(ppcif_percentiles$quant, c(0, 0.5, 1, 2)),
  is.null(ppcif_percentiles$perc),
  identical(dim(ppcif_percentiles$mixture$med_cdf_ci), c(2L, 4L)),
  ppcif_percentiles$mixture$med_cdf[1] > 0,
  identical(ppcif_percentiles$nonprevalent$med_cdf[1], 0)
)

set.seed(4402)
ppcif_quantiles <- ppCIF(
  ppcif_fit,
  fix_x_t = 0,
  pst_samples = 20,
  ppd_type = "quantiles",
  perc = c(0, 0.5, 1)
)
stopifnot(
  identical(ppcif_quantiles$ppd_type, "quantiles"),
  identical(ppcif_quantiles$perc, c(0, 0.5, 1)),
  is.null(ppcif_quantiles$quant),
  identical(dim(ppcif_quantiles$nonprevalent$med_cdf_ci), c(2L, 3L))
)

too_many_draws <- try(
  ppCIF(ppcif_fit, pst_samples = 41, quant = c(0, 1)),
  silent = TRUE
)
stopifnot(
  inherits(too_many_draws, "try-error"),
  grepl(
    "`pst_samples` cannot exceed the number of stored post-warm-up posterior draws (40).",
    as.character(too_many_draws),
    fixed = TRUE
  )
)

make_gengamma_ppcif_chain <- function(offset = 0) {
  coda::mcmc(cbind(
    beta_t_intercept = seq(0.15, 0.25, length.out = 20) + offset,
    beta_t_x_1 = seq(-0.14, -0.06, length.out = 20),
    sigma_t = seq(0.72, 0.88, length.out = 20),
    q_t = seq(-0.6, 0.6, length.out = 20),
    beta_g_intercept = seq(-1.6, -1.4, length.out = 20),
    beta_g_x_1 = seq(0.1, 0.2, length.out = 20)
  ))
}

gengamma_fit <- ppcif_fit
gengamma_fit$par <- coda::mcmc.list(
  make_gengamma_ppcif_chain(-0.01),
  make_gengamma_ppcif_chain(0.01)
)
gengamma_fit$dist <- "gengamma"

set.seed(4403)
gengamma_ppcif <- ppCIF(
  gengamma_fit,
  pst_samples = 20,
  quant = c(0, 0.5, 1, 2)
)
stopifnot(
  inherits(gengamma_ppcif, "ppCIF"),
  identical(gengamma_ppcif$distribution, "gengamma"),
  identical(gengamma_ppcif$ppd_type, "percentiles"),
  all(is.finite(gengamma_ppcif$mixture$med_cdf)),
  all(is.finite(gengamma_ppcif$nonprevalent$med_cdf)),
  all(diff(gengamma_ppcif$mixture$med_cdf) >= 0),
  all(diff(gengamma_ppcif$nonprevalent$med_cdf) >= 0)
)

plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
plot_return_mixture <- plot(ppcif_percentiles)
plot_return_nonprevalent <- plot(
  ppcif_percentiles,
  type = "nonprevalent",
  ci = FALSE
)
plot_return_both <- plot(
  ppcif_quantiles,
  type = "both",
  main = c("Mixture", "Non-prevalent")
)
plot_return_gengamma <- plot(
  gengamma_ppcif,
  type = "both"
)
grDevices::dev.off()
unlink(plot_file)

stopifnot(
  identical(plot_return_mixture, ppcif_percentiles),
  identical(plot_return_nonprevalent, ppcif_percentiles),
  identical(plot_return_both, ppcif_quantiles),
  identical(plot_return_gengamma, gengamma_ppcif)
)
