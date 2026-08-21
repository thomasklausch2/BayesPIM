#' Assemble the per-chain state passed to one Gibbs sampler run
#' @noRd
ini_gibbs = function(
    update_run,
    j,
    inis,
    p1_t,
    p1_g,
    v_obs,
    times_prev,
    g_prev,
    k_prev,
    par_prev,
    update_kappa,
    kappa,
    dist,
    fix_sigma,
    sig_prior,
    fix_q,
    q_prior_sd,
    prev,
    sampler,
    prop_sd,
    slice_width,
    rescale_factor
){
  n_par_t = p1_t + 1L + as.integer(dist == "gengamma")
  n = length(v_obs)
  chain_rows = ((n * (j - 1L)) + 1L):(n * j)
  # Only the current state is needed to advance the chain. Iterations selected
  # by `save_every` are accumulated separately in `bayespim()`.
  cur_par_t = matrix(NA_real_, nrow = 1L, ncol = n_par_t)
  beta_g = matrix(NA_real_, nrow = 1L, ncol = p1_g)
  
  if(update_run){
    # Import the exact terminal state on the sampler's internal scale.
    cur_par_t[1,1:p1_t] = par_prev[j,1:p1_t]
    cur_par_t[1,p1_t + 1L] = par_prev[j,p1_t + 1L]
    if(dist == "gengamma") cur_par_t[1,p1_t + 2L] = par_prev[j,p1_t + 2L]
    # Import prevalence model parameters.
    if(prev) {
      first_g = n_par_t + 1L
      beta_g[1,] = par_prev[j,first_g:(first_g + p1_g - 1L)]
    }
    # Import kappa parameter
    if(update_kappa) kappa = par_prev[j,ncol(par_prev)]
    # Import augmented prevalence states.
    g_aug = g_prev[chain_rows]
    # Import augmented latent times.
    if(sampler == "slice_collapsed") times = NULL else{
      times = times_prev[chain_rows] / rescale_factor
      incidence_data = list(times = times, rows = which(g_aug == 0))
    }
    # Import augmented k
    k = k_prev[chain_rows]
    # Create last selected intervals data.frame
    if (sampler == "slice_collapsed") {
      incidence_rows <- which(g_aug == 0L)
      k_active <- k[incidence_rows]
      phi_active <- look_up_mat_rcpp( v_obs[incidence_rows], k_active )
      incidence_data <- list( L = phi_active[, 1L], R = phi_active[, 2L],rows = incidence_rows )
    }
  } 
  if(!update_run){
    # Import initialized latent-time model parameters.
    cur_par_t[1,] = inis$ini_par_t[j,]
    # Import initialized prevalence model parameters.
    if(prev) beta_g[1,] = inis$ini_par_g[j,]
    # Import initialized kappa
    kappa = inis$ini_kappa[j]
    # Import initialized latent times.
    if(sampler == "slice_collapsed") times = NULL else times = inis$ini_incidence_data[[j]]$times
    # Import initialized prevalence states.
    g_aug = inis$ini_g[j,]
    # Import initialized k
    k = inis$ini_k[j,]
    # Import incidence_data
    incidence_data = inis$ini_incidence_data[[j]]
    
  } 

  if(fix_sigma) cur_par_t[1,p1_t + 1L] = log(sig_prior)

  prop_sd_mat = NULL
  if(sampler == "mh"){
    prop_sd_mat = diag(prop_sd^2, n_par_t)
  }
  
  incidence_step <- switch(
    sampler,
    mh = step_mh,
    slice = step_slice,
    slice_collapsed = step_slice_collapsed,
    stop("Unsupported incidence sampler: ", sampler)
  )
  
  incidence_control <- switch(
    sampler,
    mh = list(prop_var = prop_sd_mat),
    slice = list(width = slice_width),
    slice_collapsed = list(width = slice_width, fix_q = fix_q)
  )
  
  incidence_log_ll <- if(sampler == "slice_collapsed") ll_aft_ic else ll_aft
  
  list(
    cur_par_t = cur_par_t, 
    beta_g = beta_g,
    kappa = kappa,
    times = times,
    g_aug = g_aug,
    k = k,
    incidence_data = incidence_data,
    incidence_step = incidence_step,
    incidence_control = incidence_control,
    incidence_log_ll  = incidence_log_ll
    )
}
