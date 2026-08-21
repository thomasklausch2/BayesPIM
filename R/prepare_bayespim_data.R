#' Recode and rescale the screening data into internal sampler inputs
#' @noRd
prepare_bayespim_data <- function(
    v_obs,
    x_t,
    x_g,
    kappa,
    update_kappa,
    rescale_times,
    standardize_covariates,
    covariate_scaling = NULL,
    tau_g,
    prev
) {
  v_obs_original <- v_obs
  
  # Recode baseline-positive observations for internal calculations.
  g_fixed <- lengths(v_obs) == 1L
  if (any(g_fixed)) {
    v_obs[g_fixed] <- rep(list(c(0, Inf)), sum(g_fixed))
  }
  
  pobs <- if (update_kappa) {
    NULL
  } else {
    p_v_obs(v_obs, kappa = kappa)
  }
  
  d <- vapply(
    v_obs,
    function(x) as.integer(is.finite(x[length(x)])) + 1L,
    integer(1)
  )
  L <- vapply(v_obs, function(x) x[length(x) - 1L], numeric(1))
  R <- vapply(v_obs, function(x) x[length(x)], numeric(1))
  
  if (rescale_times) {
    max_time <- R
    max_time[is.infinite(max_time)] <- L[is.infinite(max_time)]
    max_time <- max_time[max_time > 0]
    
    if (!length(max_time)) {
      stop(
        "Cannot rescale screening times: no positive finite time is available.",
        call. = FALSE
      )
    }
    
    times_rescale_factor <- median(max_time)
  } else {
    times_rescale_factor <- 1
  }
  
  v_obs <- lapply(v_obs, function(x) x / times_rescale_factor)
  L <- vapply(v_obs, function(x) x[length(x) - 1L], numeric(1))
  R <- vapply(v_obs, function(x) x[length(x)], numeric(1))
  
  v_obs_l <- unlist(lapply(
    v_obs,
    function(x) x[seq_len(length(x) - 1L)]
  ))
  v_obs_r <- unlist(lapply(
    v_obs,
    function(x) x[seq.int(2L, length(x))]
  ))
  
  if (is.null(x_t)) {
    x_t_original <- NULL
    x1_t <- matrix(1, nrow = length(d), ncol = 1L)
  } else {
    x_t_original <- as.matrix(x_t)
    
    if (is.null(colnames(x_t_original))) {
      colnames(x_t_original) <- paste0("x_t_", seq_len(ncol(x_t_original)))
    }
  }
  
  if (is.null(x_g)) {
    x_g_original <- NULL
    x1_g <- matrix(1, nrow = length(d), ncol = 1L)
  } else {
    x_g_original <- as.matrix(x_g)
    
    if (is.null(colnames(x_g_original))) {
      colnames(x_g_original) <- paste0("x_g_", seq_len(ncol(x_g_original)))
    }
  }

  if (is.null(covariate_scaling)) {
    covariate_scaling <- list(
      x_t = bayespim_covariate_scaling(
        x_t_original,
        standardize = standardize_covariates,
        name = "x_t"
      ),
      x_g = bayespim_covariate_scaling(
        x_g_original,
        standardize = standardize_covariates,
        name = "x_g"
      )
    )
  }

  x_t <- bayespim_standardize_matrix(
    x_t_original,
    covariate_scaling$x_t,
    "x_t"
  )
  x_g <- bayespim_standardize_matrix(
    x_g_original,
    covariate_scaling$x_g,
    "x_g"
  )
  if (!is.null(x_t)) x1_t <- cbind(1, x_t)
  if (!is.null(x_g)) x1_g <- cbind(1, x_g)
  
  p1_t <- ncol(x1_t)
  p1_g <- ncol(x1_g)
  
  # Needed for beta_g updates whenever the prevalence model is active.
  sig_inv <- NULL
  sig_inv_xt <- NULL
  
  if (prev) {
    sig_inv <- solve(
      crossprod(x1_g) +
        diag(rep(tau_g^-2, ncol(x1_g)))
    )
    sig_inv_xt <- sig_inv %*% t(x1_g)
  }
  
  list(
    v_obs_original = v_obs_original,
    v_obs = v_obs,
    g_fixed = g_fixed,
    pobs = pobs,
    d = d,
    L = L,
    R = R,
    times_rescale_factor = times_rescale_factor,
    v_obs_l = v_obs_l,
    v_obs_r = v_obs_r,
    x_t_original = x_t_original,
    x_t = x_t,
    x1_t = x1_t,
    x_g_original = x_g_original,
    x_g = x_g,
    x1_g = x1_g,
    covariate_scaling = covariate_scaling,
    p1_t = p1_t,
    p1_g = p1_g,
    sig_inv = sig_inv,
    sig_inv_xt = sig_inv_xt
  )
}
