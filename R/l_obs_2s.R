#' Observed-data likelihood contributions for a single posterior draw
#' @noRd
l_obs_2s = function (est, mod, log_scale = T, sumup = T) 
{
  prev = mod$prev
  v_obs = mod$v_obs
  g_fixed = numeric(length(v_obs))
  for (i in 1:length(g_fixed)) g_fixed[i] = as.numeric(length(v_obs[[i]]) == 
                                                         1)
  for (i in 1:length(v_obs)) if (as.logical(g_fixed[i])) 
    v_obs[[i]] = c(0, Inf)
  dist = mod$dist
  x1_t = if (is.null(mod$x_t)) {
    matrix(1, nrow = length(v_obs), ncol = 1L)
  } else {
    cbind(1, mod$x_t)
  }
  p1_t = ncol(x1_t)
  m = sapply(v_obs, length) - 1
  # Number of incidence shape parameters beyond the AFT coefficients and
  # log-sigma. The generalized gamma adds the signed shape parameter Q, which
  # shifts the position of the prevalence coefficients within `est`.
  n_extra_t = as.integer(dist == "gengamma")
  if (prev) {
    x1_g = if (is.null(mod$x_g)) {
      matrix(1, nrow = length(v_obs), ncol = 1L)
    } else {
      cbind(1, mod$x_g)
    }
    p1_g = ncol(x1_g)
    kappa = est[(length(est))]
    est_t = est[1:(p1_t + 1 + n_extra_t)]
    est_g = est[(p1_t + 2 + n_extra_t):(p1_t + p1_g + 1 + n_extra_t)]
  }
  if (!prev) {
    kappa = est[(length(est))]
    est_t = est[1:(length(est) - 1)]
  }
  pobs_vec = unlist(p_v_obs(v_obs, kappa = kappa))
  v_obs_l = unlist(lapply(v_obs, function(x) x[1:(length(x) - 
                                                  1)]))
  v_obs_r = unlist(lapply(v_obs, function(x) x[2:(length(x))]))
  cur_dist_par_t = switch(
    dist,
    weibull   = trans_par(x1_t, par = est_t),
    loglog    = trans_par(x1_t, par = est_t),
    gamma     = trans_par_gamma(x1_t, par = est_t),
    lognormal = trans_par_ind_norm(
      x1 = x1_t,
      p = est_t[seq_len(length(est_t) - 1L)],
      v = est_t[length(est_t)]
    ),
    gengamma  = trans_par_gengamma(x1_t, par = est_t),
    stop("Unsupported distribution: ", dist)
  )
  # Expand each subject's distribution parameters to one row per screening
  # interval, so `par` aligns with the flattened interval endpoints. This is
  # agnostic to the number of parameter columns (2 for most families, 3 for
  # the generalized gamma).
  par = cur_dist_par_t[rep(seq_len(nrow(cur_dist_par_t)), times = m), , drop = FALSE]
  Fxl = pdist(v_obs_l, par = par, dist = dist)
  Fxr = pdist(v_obs_r, par = par, dist = dist)
  pobs_ = (Fxr - Fxl) * pobs_vec
  pobs_ = split(x = pobs_, f = rep(1:length(m), m))
  if (prev) {
    mu_w = as.numeric(x1_g %*% as.matrix(as.numeric(est_g)))
    theta = pnorm(mu_w)
    q0 = sapply(pobs_, sum) * (1 - theta) * (1 - g_fixed)
    q1 = p_v_obs_2(v_obs, kappa, r = mod$r) * (theta) * (1 - g_fixed)
    S = q0 + q1
    S[g_fixed == 1] = (kappa * theta)[g_fixed == 1]
  }
  if (!prev) {
    S = sapply(pobs_, sum)
  }
  if (log_scale & sumup) 
    r = sum(log(S))
  if (!log_scale & !sumup) 
    r = S
  if (log_scale & !sumup) 
    r = log(S)
  r
}
