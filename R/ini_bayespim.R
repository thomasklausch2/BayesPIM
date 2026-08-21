#' Construct per-chain starting values for the Gibbs sampler
#' @noRd
ini_bayespim = function(v_obs, x_t, x1_t, x1_g, p1_t, r, g_fixed, chains, kappa, update_kappa, 
                        pobs, v_obs_l, v_obs_r, tau_g = 1, kappa_ab = NULL, fix_sigma = FALSE, 
                        seed_chains, sampler, prev, log_prior_fun, tau_t, sig_prior, beta_prior,
                        dist, spread_lower = 0, q_prior_sd, fix_q){
  
 fit_ini_par_g = function(x1_g, g_fixed, r, kappa, tau_g){
    tested = r == 1L
    x = x1_g[tested, , drop = FALSE]
    y = as.numeric(g_fixed[tested])
    p_g = ncol(x1_g)
    if(nrow(x) == 0L) return(list(coef = rep(0, p_g), vcov = diag(tau_g^2, p_g), convergence = 0L))
    eps = 1e-10
    kappa_fit = pmin(pmax(kappa, eps), 1 - eps)
    objective = function(beta){
      theta = pnorm(drop(x %*% beta))
      p_detected = pmin(pmax(kappa_fit * theta, eps), 1 - eps)
      -sum(stats::dbinom(y, size = 1, prob = p_detected, log = TRUE)) + sum(beta^2) / (2 * tau_g^2)
    }
    p_detected_start = (sum(y) + 0.5) / (length(y) + 1)
    theta_start = pmin(pmax(p_detected_start / kappa_fit, 1e-4), 1 - 1e-4)
    beta_start = c(qnorm(theta_start), rep(0, p_g - 1L))
    fit = optim(par = beta_start, fn = objective, method = "BFGS", hessian = TRUE, control = list(maxit = 1000, reltol = 1e-10))
    H = (fit$hessian + t(fit$hessian)) / 2
    eig = eigen(H, symmetric = TRUE)
    vcov_fit = eig$vectors %*% diag(1 / pmax(eig$values, 1e-8), nrow = length(eig$values)) %*% t(eig$vectors)
    list(coef = fit$par, vcov = vcov_fit, convergence = fit$convergence)
  }

  #prep
  v = v_obs[!g_fixed]
  v = t(vapply(v, function(x) x[c(length(x) - 1L, length(x))], numeric(2)))
  v[v[,1] == 0,1] = -Inf
  dist_survreg = switch(dist, weibull = "weibull", 
                        loglog = "loglogistic", 
                        lognormal = "lognormal", 
                        gamma = 'lognormal',
                        gengamma = 'lognormal')

  # Deterministic quantities shared by the chain-specific initializations
  x_t_nonfixed = if(is.null(x_t)) NULL else x_t[!g_fixed, , drop = FALSE]
  survival_formula = if(is.null(x_t_nonfixed)) survival::Surv(time = v[,1], time2 = v[,2], type = "interval2") ~ 1 else survival::Surv(time = v[,1], time2 = v[,2], type = "interval2") ~ x_t_nonfixed
  mod = survival::survreg(survival_formula, dist = dist_survreg)
  betas = stats::coef(mod)
  logsigma = log(mod$scale)
  
  # Rescaling in case of gamma
  if(dist == 'gamma'){
    surrogate_variance <- mod$scale^2
    betas[1L] <- betas[1L] + surrogate_variance / 2
    logsigma <- 0.5 * log(expm1(surrogate_variance))
  }
  
  if(update_kappa && is.null(kappa_ab)) stop("kappa_ab must be supplied when update_kappa = TRUE.")

  # Each chain's seed governs all of its random starting values. The state after
  # initialization is returned so that the worker continues the same RNG stream.
  ini_par_t = matrix(NA_real_, nrow = chains, ncol = length(c(betas, logsigma)) + 
                     as.integer(dist == 'gengamma'))
  sigma_pos = length(c(betas, logsigma))
  ini_kappa = numeric(chains)
  ini_par_g = matrix(NA_real_, nrow = chains, ncol = ncol(x1_g))
  ini_incidence_data = list()
  ini_g     = matrix(nrow = chains, ncol = nrow(x1_t))
  ini_k     = matrix(nrow = chains, ncol = nrow(x1_t))
  rng_state = vector("list", chains)
  times = numeric(nrow(x1_t))

  spread <- if (chains == 1L) 1 else seq(spread_lower, 1, length.out = chains)
  par_ml_t <- c(betas, logsigma)
 
  for(j in seq_len(chains)){
    set.seed(seed_chains[j])
    tries_left = 10
    try_ini = TRUE
    while(try_ini & tries_left > 0){
      tries_left = tries_left - 1
      # Latent-time model parameters.
      ini_par_t[j,1:sigma_pos] = spread[j] * (par_ml_t )
      if(fix_sigma) ini_par_t[j,sigma_pos] = log(sig_prior)
      
      # Generalized gamma Q
      if(dist == 'gengamma'){
        ini_par_t[j, ncol(ini_par_t)] <- 1 # add jitter later
        if(fix_q) ini_par_t[j, ncol(ini_par_t)] <- q_prior_sd
      }
  
      # kappa
      if(update_kappa){
        if(j == chains) ini_kappa[j] = 1 else ini_kappa[j] = runif(1, 0.5, 1)
      } else ini_kappa[j] = kappa
      
      # Prevalence model parameters.
      fit_w = fit_ini_par_g(x1_g = x1_g, g_fixed = g_fixed, r = r, kappa = ini_kappa[j], tau_g = tau_g)
      par_ml_g <- fit_w$coef
      if(fit_w$convergence != 0L) stop("Prevalence starting-value optimization failed to converge.")
      ini_par_g[j,] = spread[j] * (par_ml_g )
  
      prob_g_j = pnorm(drop(x1_g %*% ini_par_g[j,]))
      kappa_j = ini_kappa[j]
      
      # Prevalence state.
      if(dist == 'weibull' | dist == 'loglog') cur_dist_par_t = trans_par(x1_t, par = ini_par_t[j,])
      if(dist == 'lognormal') cur_dist_par_t = trans_par_ind_norm(x1 = x1_t, p = ini_par_t[j,1:p1_t], v= ini_par_t[j,(p1_t+1)])
      if(dist == 'gamma')  cur_dist_par_t = trans_par_gamma(x1 = x1_t, par = ini_par_t[j,])
      if(dist == 'gengamma')  cur_dist_par_t = trans_par_gengamma(x1 = x1_t, par = ini_par_t[j,])
      if(update_kappa) pobs = p_v_obs_rcpp(v_obs, kappa = kappa_j)
      
      omega = interval_probs_rcpp(pobs_vec = unlist(pobs), v_obs, v_obs_l, v_obs_r,
                                  cur_dist_par_t, dist)
      interval_sums          = omega$interval_sums
      interval_probabilities = omega$pobs_norm
      if(prev) g_j <- augment_g_collapsed_rcpp(interval_sums = interval_sums, v_obs = v_obs, kappa = kappa_j,
                                                 prob_g = prob_g_j, r = r, g_fixed = g_fixed )
      if(!prev) g_j <- rep(0, nrow(x1_t))
      ini_g[j,] = g_j
      
      # k
      k_j = sapply(v_obs, function(x) length(x)-1)
      incidence_rows <- which(g_j == 0L)
      prevalence_rows <- which(g_j == 1L)
      k_active_j <- sample_k_rcpp( interval_probabilities[incidence_rows] )
      phi_active_j <- look_up_mat_rcpp( v_obs[incidence_rows], k_active_j )
      k_j[incidence_rows] <- k_active_j
      ini_k[j,] <- k_j
      
      # Define incidence_data and evaluate posterior
      if (sampler == "slice_collapsed") {
        incidence_data <- list( L = phi_active_j[, 1L], R = phi_active_j[, 2L],rows = incidence_rows )
        
        pst <- log_pst(eta = ini_par_t[j,], log_ll = ll_aft_ic, log_prior_fun = log_prior_fun, 
                tau_t = tau_t, sig_prior = sig_prior, q_prior_sd = q_prior_sd, beta_prior = beta_prior,
                L = incidence_data$L, R = incidence_data$R, x = x1_t[incidence_rows,,drop=FALSE], 
                dist = dist)
        
      } else if (sampler %in% c("mh", "slice")) {
        if (length(incidence_rows) > 0L) {
          times[incidence_rows] <- r_trdist( par = cur_dist_par_t[incidence_rows, , drop = FALSE],
                                         a = phi_active_j[, 1L], b = phi_active_j[, 2L], dist = dist ) 
        }
        if (length(prevalence_rows) > 0L) { 
          times[prevalence_rows] <- rdist( n = length(prevalence_rows), par = cur_dist_par_t[prevalence_rows, , drop = FALSE], dist = dist ) 
        }
        times[times == 0] <- 1e-300
        incidence_data <- list(times = times, rows = incidence_rows)
        pst <- log_pst(eta = ini_par_t[j,], log_ll = ll_aft, log_prior_fun = log_prior_fun, 
                       tau_t = tau_t, sig_prior = sig_prior, q_prior_sd = q_prior_sd, beta_prior = beta_prior, 
                       y = times[incidence_rows], x = x1_t[incidence_rows,,drop=FALSE], 
                       dist = dist) 
        
      }
        ini_incidence_data[[j]] = incidence_data
        try_ini = !is.finite(pst)
        if(tries_left == 2) spread[j] = 1
        if(tries_left == 1) spread[j] = 0
    }
    if (!is.finite(pst)) stop('Cannot find starting values. You can try increasing or descreasing ini_spread. Check the input data closely.')
    rng_state[[j]] = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  list(ini_par_t = ini_par_t, ini_par_g = ini_par_g, ini_kappa = ini_kappa, 
       ini_incidence_data = ini_incidence_data,
       ini_g = ini_g, ini_k = ini_k, rng_state = rng_state)
}
