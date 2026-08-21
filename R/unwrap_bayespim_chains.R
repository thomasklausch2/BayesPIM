#' Assemble raw sampler output into named mcmc chains
#' @noRd
unwrap_bayespim_chains <- function(
    run,
    times_rescale_factor,
    p1_t,
    x_t,
    x_g,
    prev,
    has_q,
    update_kappa,
    sampler,
    saved_iterations = seq_len(nrow(as.matrix(run[[1L]]$par))),
    save_every = 1L,
    covariate_scaling = NULL,
    prev_run = NULL,
    previous_par = NULL
) {
  if (!length(run)) {
    stop("`run` must contain at least one MCMC chain.", call. = FALSE)
  }
  
  parameter_names <- if (prev) {
    if(!has_q){
    c(
      "beta_t_intercept",
      if (is.null(x_t)) NULL else paste0("beta_t_", colnames(x_t)),
      "sigma_t",
      "beta_g_intercept",
      if (is.null(x_g)) NULL else paste0("beta_g_", colnames(x_g))
    )} else{
      c(
        "beta_t_intercept",
        if (is.null(x_t)) NULL else paste0("beta_t_", colnames(x_t)),
        "sigma_t", "q_t",
        "beta_g_intercept",
        if (is.null(x_g)) NULL else paste0("beta_g_", colnames(x_g))
      )
    }
  } else {
    if(!has_q){
      c(
        "beta_t_intercept",
        if (is.null(x_t)) NULL else paste0("beta_t_", colnames(x_t)),
        "sigma_t"
      )
    } else{
      c(
        "beta_t_intercept",
        if (is.null(x_t)) NULL else paste0("beta_t_", colnames(x_t)),
        "sigma_t", "q_t"
      )
    }
    
  }
  
  # Initial values, with the incidence intercept returned to user time units.
  inis <- do.call(
    rbind,
    lapply(run, function(chain) chain$ini)
  )
  inis <- bayespim_transform_coefficients(
    par = inis,
    p1_t = p1_t,
    p1_g = if (is.null(x_g)) 1L else ncol(x_g) + 1L,
    prev = prev,
    has_q = has_q,
    covariate_scaling = covariate_scaling,
    direction = "internal_to_public"
  )
  inis[, 1L] <- inis[, 1L] + log(times_rescale_factor)
  
  # Transform and label each parameter chain.
  public_parameter_matrix <- function(par, include_kappa = FALSE) {
    par <- as.matrix(par)

    par <- bayespim_transform_coefficients(
      par = par,
      p1_t = p1_t,
      p1_g = if (is.null(x_g)) 1L else ncol(x_g) + 1L,
      prev = prev,
      has_q = has_q,
      covariate_scaling = covariate_scaling,
      direction = "internal_to_public"
    )

    # Return the incidence intercept and sigma to their public scales.
    par[, 1L] <- par[, 1L] + log(times_rescale_factor)
    par[, p1_t + 1L] <- exp(par[, p1_t + 1L])

    output_names <- c(parameter_names, if (include_kappa) "kappa")
    if (length(output_names) != ncol(par)) {
      stop(
        sprintf(
          paste0(
            "Parameter-name mismatch: generated %d names for ",
            "a parameter matrix with %d columns."
          ),
          length(output_names),
          ncol(par)
        ),
        call. = FALSE
      )
    }

    colnames(par) <- output_names
    par
  }

  if (length(saved_iterations) > 0L) {
    mcmc_par <- lapply(run, function(chain) {
      par <- public_parameter_matrix(chain$par)
      if (update_kappa) {
        par <- cbind(par, kappa = chain$kappa)
      }
      mcmc(
        par,
        start = saved_iterations[1L],
        end = saved_iterations[length(saved_iterations)],
        thin = save_every
      )
    })
    par <- mcmc.list(mcmc_par)
  } else {
    par <- NULL
  }

  terminal_par_internal <- do.call(
    rbind,
    lapply(run, function(chain) as.matrix(chain$terminal_par))
  )
  colnames(terminal_par_internal) <- c(
    parameter_names,
    if (update_kappa) "kappa"
  )
  terminal_par <- public_parameter_matrix(
    terminal_par_internal,
    include_kappa = update_kappa
  )
  
  # Acceptance indicators are matrices with draws in rows and chains
  # in columns. For slice samplers all elements are NULL, so this
  # remains NULL.
  ac <- do.call(
    cbind,
    lapply(run, function(chain) chain$ac)
  )
  
  # These are terminal subject-level states, concatenated by chain.
  times_draw <- if (sampler == "slice_collapsed") {
    NULL
  } else {
    unlist(
      lapply(run, function(chain) chain$times),
      use.names = FALSE
    )
  }
  
  k_draw <- unlist(
    lapply(run, function(chain) chain$k),
    use.names = FALSE
  )
  
  g_draw <- unlist(
    lapply(run, function(chain) chain$g_aug),
    use.names = FALSE
  )
  
  ac_cur <- NULL
  
  if (!is.null(prev_run)) {
    if (is.null(previous_par)) {
      stop(
        "`previous_par` must be supplied when updating `prev_run`.",
        call. = FALSE
      )
    }
    
    if (!is.null(par)) {
      par <- bind_mcmc_lists(previous_par, par)
    } else {
      par <- previous_par
    }
    
    # Preserve the acceptance indicators from this update separately
    # before joining them to the preceding run.
    ac_cur <- ac
    ac <- rbind(prev_run$ac, ac)
  }
  
  list(
    inis = inis,
    par = par,
    terminal_par = terminal_par,
    terminal_par_internal = terminal_par_internal,
    ac = ac,
    ac_cur = ac_cur,
    times_draw = times_draw,
    k_draw = k_draw,
    g_draw = g_draw
  )
}
