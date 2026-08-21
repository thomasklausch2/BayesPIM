#' Input validation for bayespim, internal use
#' @noRd
validate_bayespim_inputs <- function(
    # Group 1: Data inputs
    v_obs = NULL,
    x_t = NULL,
    x_g = NULL,
    r = NULL,

    # Group 2: Basic model settings
    dist = "weibull",
    kappa = NULL,
    update_kappa = FALSE,
    kappa_prior = NULL,

    # Group 3: Main MCMC sampler settings
    ndraws = 1000,
    warmup = NULL,
    warmup_updated = FALSE,
    prop_sd = NULL,
    slice_width = 1,
    chains = NULL,
    seed_chains = NULL,
    save_every = 1,
    update_till_converge = FALSE,
    maxit = Inf,
    max_rhat = 1.01,
    min_effss = NULL,
    sampler = "slice_collapsed",

    # Group 4: Prior specifications
    log_prior_fun = log_aft_prior,
    beta_prior = "norm",
    tau_t = 1,
    sig_prior = 1,
    tau_g = 1,
    fix_sigma = FALSE,
    q_prior_sd = 1,
    fix_q = FALSE,

    # Group 5: Updating previous run settings
    prev_run = NULL,
    ndraws_update = NULL,

    # Group 6: Additional model flags
    prev = TRUE,
    par_exp = FALSE,
    rescale_times = TRUE,
    standardize_covariates = TRUE,
    ini_spread = 0.5,

    # Group 7: Console output
    silent = FALSE,

    # Internal validation control
    stage = c("effective", "initial")
) {
  stage <- match.arg(stage)
  errors <- character(0)
  
  add_error <- function(x) {
    errors <<- c(errors, x)
  }
  
  stop_if_errors <- function() {
    if (length(errors) > 0L) {
      stop(
        paste0(
          "Invalid input to `bayespim()`:\n",
          paste0("  - ", errors, collapse = "\n")
        ),
        call. = FALSE
      )
    }
    invisible(TRUE)
  }
  
  fmt <- function(x, max_n = 6L) {
    x <- as.character(x)
    if (length(x) > max_n) {
      paste0(paste(x[seq_len(max_n)], collapse = ", "), ", ...")
    } else {
      paste(x, collapse = ", ")
    }
  }
  
  is_logical1 <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }
  
  is_num1 <- function(x, finite = TRUE, positive = FALSE, nonnegative = FALSE) {
    if (!(is.numeric(x) && length(x) == 1L && !is.na(x))) {
      return(FALSE)
    }
    if (finite && !is.finite(x)) {
      return(FALSE)
    }
    if (positive && !(x > 0)) {
      return(FALSE)
    }
    if (nonnegative && !(x >= 0)) {
      return(FALSE)
    }
    TRUE
  }
  
  is_count1 <- function(x, min = 1L) {
    is_num1(x, finite = TRUE) &&
      abs(x - round(x)) < sqrt(.Machine$double.eps) &&
      x >= min
  }

  is_count_or_inf1 <- function(x, min = 1L) {
    is_count1(x, min = min) ||
      (is.numeric(x) && length(x) == 1L && !is.na(x) &&
         is.infinite(x) && x > 0)
  }
  
  check_logical <- function(name, x) {
    if (!is_logical1(x)) {
      add_error(sprintf("`%s` must be a single TRUE or FALSE value.", name))
    }
  }
  
  check_count <- function(name, x, min = 1L) {
    if (!is_count1(x, min = min)) {
      add_error(sprintf(
        "`%s` must be a single whole number >= %s.",
        name, min
      ))
    }
  }
  
  check_positive_num <- function(name, x) {
    if (!is_num1(x, finite = TRUE, positive = TRUE)) {
      add_error(sprintf("`%s` must be a single positive finite number.", name))
    }
  }
  
  check_design_matrix <- function(x, name, n) {
    if (is.null(x)) {
      return(invisible(TRUE))
    }
    
    if (is.data.frame(x)) {
      non_numeric <- names(x)[!vapply(x, is.numeric, logical(1))]
      if (length(non_numeric) > 0L) {
        add_error(sprintf(
          "`%s` must contain only numeric columns; non-numeric column(s): %s.",
          name, fmt(non_numeric)
        ))
        return(invisible(FALSE))
      }
    }
    
    x_matrix <- tryCatch(
      as.matrix(x),
      error = function(e) e
    )
    
    if (inherits(x_matrix, "error")) {
      add_error(sprintf("`%s` could not be converted to a matrix.", name))
      return(invisible(FALSE))
    }
    
    if (!is.numeric(x_matrix)) {
      add_error(sprintf("`%s` must be numeric.", name))
      return(invisible(FALSE))
    }
    
    if (!is.na(n) && nrow(x_matrix) != n) {
      add_error(sprintf(
        "`%s` has %s row(s), but `v_obs` has length %s. These must match.",
        name, nrow(x_matrix), n
      ))
    }
    
    if (ncol(x_matrix) < 1L) {
      add_error(sprintf(
        "`%s` has zero columns. Use NULL for an intercept-only model.",
        name
      ))
    }
    
    if (anyNA(x_matrix)) {
      add_error(sprintf("`%s` contains missing values; missing covariates are not allowed.", name))
    }
    
    if (any(!is.finite(x_matrix))) {
      add_error(sprintf("`%s` must contain only finite numeric values.", name))
    }
    
    invisible(TRUE)
  }
  
  ## Validate previous run before bayespim dereferences it
  if (!is.null(prev_run)) {
    if (!is.list(prev_run)) {
      add_error("`prev_run` must be a previous `bayespim()` result object.")
      stop_if_errors()
    }

    # Every component `bayespim()` dereferences when continuing a run.
    required_prev <- c(
      "ini", "par", "terminal_par_internal", "times", "k", "g", "ac",
      "v_obs", "x_t", "x_g", "r",
      "dist", "kappa", "update_kappa", "kappa_prior",
      "ndraws", "warmup", "warmup_updated", "prop_sd", "slice_width",
      "seed_chains", "rng_state", "save_every", "total_iterations",
      "maxit", "max_rhat", "sampler",
      "log_prior_fun", "beta_prior", "priors", "fix_sigma", "fix_q",
      "prev", "par_exp", "rescale_times", "standardize_covariates",
      "covariate_scaling", "runtime"
    )

    missing_prev <- setdiff(required_prev, names(prev_run))
    if (length(missing_prev) > 0L) {
      add_error(sprintf(
        "`prev_run` is missing required component(s): %s.",
        fmt(paste0("$", missing_prev))
      ))
    }
    stop_if_errors()

    # From here on every component is present, so only its content is checked.
    if (!is.list(prev_run$priors)) {
      add_error("`prev_run$priors` must be a list.")
    } else {
      missing_priors <- setdiff(
        c("tau_t", "sig_prior", "tau_g", "q_prior_sd"),
        names(prev_run$priors)
      )
      if (length(missing_priors) > 0L) {
        add_error(sprintf(
          "`prev_run$priors` is missing required component(s): %s.",
          fmt(paste0("$", missing_priors))
        ))
      }
    }

    previous_par <- prev_run$par
    if (!(inherits(previous_par, "mcmc.list") || is.list(previous_par))) {
      add_error("`prev_run$par` must be an `mcmc.list` or list-like MCMC object.")
    } else if (length(previous_par) < 1L) {
      add_error("`prev_run$par` contains no chains.")
    }

    if (!is.list(prev_run$covariate_scaling) ||
        !all(c("x_t", "x_g") %in% names(prev_run$covariate_scaling))) {
      add_error("`prev_run$covariate_scaling` must be a list with `$x_t` and `$x_g` scaling records.")
    }
    stop_if_errors()

    n_chains_prev <- length(previous_par)
    n_states <- if (is.list(prev_run$v_obs)) {
      length(prev_run$v_obs) * n_chains_prev
    } else {
      NA_integer_
    }

    if (!is_count1(prev_run$total_iterations, min = 1L)) {
      add_error("`prev_run$total_iterations` must be a single positive whole number.")
    }

    if (!is_count1(prev_run$save_every, min = 1L)) {
      add_error("`prev_run$save_every` must be a single positive whole number.")
    }

    terminal_par <- prev_run$terminal_par_internal
    terminal_ok <- is.numeric(terminal_par) && is.matrix(terminal_par) &&
      nrow(terminal_par) == n_chains_prev && ncol(terminal_par) >= 2L &&
      !anyNA(terminal_par) && all(is.finite(terminal_par))
    if (!terminal_ok) {
      add_error("`prev_run$terminal_par_internal` must contain one finite internal-scale parameter row per chain.")
    }

    k_ok <- is.numeric(prev_run$k) && length(prev_run$k) >= 1L &&
      !anyNA(prev_run$k) && all(is.finite(prev_run$k)) &&
      all(prev_run$k == as.integer(prev_run$k)) && all(prev_run$k >= 1L) &&
      (is.na(n_states) || length(prev_run$k) == n_states)
    if (!k_ok) {
      add_error("`prev_run$k` must contain one positive whole-number interval index per subject and chain.")
    }

    g_ok <- is.numeric(prev_run$g) && !anyNA(prev_run$g) &&
      all(prev_run$g %in% c(0, 1)) &&
      (is.na(n_states) || length(prev_run$g) == n_states)
    if (!g_ok) {
      add_error("`prev_run$g` must contain one binary prevalence state per subject and chain.")
    }

    if (!identical(prev_run$sampler, "slice_collapsed")) {
      times_ok <- is.numeric(prev_run$times) && !anyNA(prev_run$times) &&
        all(is.finite(prev_run$times)) && all(prev_run$times > 0) &&
        (is.na(n_states) || length(prev_run$times) == n_states)
      if (!times_ok) {
        add_error("`prev_run$times` must contain one positive finite latent time per subject and chain for an exact-time sampler.")
      }
    }

    rng_state_ok <- is.list(prev_run$rng_state) &&
      length(prev_run$rng_state) == n_chains_prev &&
      all(vapply(
        prev_run$rng_state,
        function(x) {
          is.numeric(x) && length(x) >= 2L && !anyNA(x) &&
            all(is.finite(x)) && all(x == as.integer(x))
        },
        logical(1)
      ))
    if (!rng_state_ok) {
      add_error("`prev_run$rng_state` must contain one valid saved RNG state per chain.")
    }
  }

  if (!is.null(seed_chains)) {
    seeds_ok <- is.numeric(seed_chains) && length(seed_chains) >= 1L &&
      !anyNA(seed_chains) && all(is.finite(seed_chains)) &&
      all(seed_chains == as.integer(seed_chains)) &&
      all(seed_chains >= 0) && all(seed_chains <= .Machine$integer.max)

    if (!seeds_ok) {
      add_error(sprintf(
        "`seed_chains` must contain whole numbers between 0 and %s.",
        .Machine$integer.max
      ))
    } else if (anyDuplicated(seed_chains)) {
      add_error("`seed_chains` must contain a unique seed for each chain.")
    }

    expected_chains <- if (is.null(prev_run)) chains else length(prev_run$par)
    if (is_count1(expected_chains) && length(seed_chains) != expected_chains) {
      add_error(sprintf(
        "`seed_chains` must have length %s, one seed per chain.",
        expected_chains
      ))
    }
  }
  
  if (identical(stage, "initial")) {
    if (is.null(prev_run)) {
      if (is.null(v_obs)) {
        add_error("`v_obs` must be supplied unless `prev_run` is supplied.")
      }
      if (is.null(chains)) {
        add_error("`chains` must be supplied unless `prev_run` is supplied.")
      }
    }
    
    if (!is_count_or_inf1(maxit, min = 1L)) {
      add_error("`maxit` must be a single positive whole number or Inf.")
    }
    
    if (!is.null(ndraws_update)) {
      check_count("ndraws_update", ndraws_update, min = 2L)
    }

    if (is.null(prev_run) &&
        is_count1(ndraws, min = 2L) &&
        is_count1(warmup, min = 0L) &&
        warmup >= ndraws) {
      add_error("`warmup` must be smaller than `ndraws` for an initial run.")
    }
    
    stop_if_errors()

    if (is.null(prev_run) &&
        isTRUE(update_till_converge) &&
        is_count1(ndraws, min = 2L) &&
        is_count1(maxit, min = 1L) &&
        ndraws > maxit) {
      warning(
        "`ndraws` exceeds `maxit`; the initial run will still produce `ndraws` draws, but no automatic convergence update can reduce it to `maxit`.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }
  
  ## Logical scalars
  logical_args <- list(
    update_kappa = update_kappa,
    warmup_updated = warmup_updated,
    update_till_converge = update_till_converge,
    fix_sigma = fix_sigma,
    fix_q = fix_q,
    prev = prev,
    par_exp = par_exp,
    rescale_times = rescale_times,
    standardize_covariates = standardize_covariates,
    silent = silent
  )
  
  for (nm in names(logical_args)) {
    check_logical(nm, logical_args[[nm]])
  }
  
  ## Choices
  allowed_dist <- c("weibull", "lognormal", "loglog", "gamma", "gengamma")
  dist_ok <- is.character(dist) && length(dist) == 1L &&
    !is.na(dist) && dist %in% allowed_dist
  if (!dist_ok) {
    add_error(sprintf(
      "`dist` must be one of: %s.",
      fmt(allowed_dist)
    ))
  }
  
  allowed_beta_prior <- c("norm", "t")
  if (!(is.character(beta_prior) && length(beta_prior) == 1L && !is.na(beta_prior) &&
        beta_prior %in% allowed_beta_prior)) {
    add_error("`beta_prior` must be either 'norm' or 't'.")
  }
  
  allowed_sampler <- c("mh", "slice", "slice_collapsed")
  sampler_ok <- is.character(sampler) && length(sampler) == 1L &&
    !is.na(sampler) && sampler %in% allowed_sampler
  if (!sampler_ok) {
    add_error("`sampler` must be one of 'mh', 'slice', or 'slice_collapsed'.")
  }

  if (!is.function(log_prior_fun)) {
    add_error("`log_prior_fun` must be a function.")
  }

  if (dist_ok && sampler_ok &&
      identical(dist, "gengamma") &&
      !identical(sampler, "slice_collapsed")) {
    add_error(
      "`dist = \"gengamma\"` requires `sampler = \"slice_collapsed\"`."
    )
  }

  if (isTRUE(fix_q) &&
      !(identical(dist, "gengamma") &&
        identical(sampler, "slice_collapsed"))) {
    add_error(
      "`fix_q = TRUE` requires `dist = \"gengamma\"` and `sampler = \"slice_collapsed\"`."
    )
  }

  ## Counts and numeric tuning parameters
  check_count("ndraws", ndraws, min = 2L)
  if (!is.null(warmup)) {
    check_count("warmup", warmup, min = 0L)
  }
  check_count("chains", chains, min = 1L)
  check_count("save_every", save_every, min = 1L)
  if (!is.null(ndraws_update)) {
    check_count("ndraws_update", ndraws_update, min = 2L)
  }
  if (is.null(prev_run) &&
      is_count1(ndraws, min = 2L) &&
      is_count1(warmup, min = 0L) &&
      warmup >= ndraws) {
    add_error("`warmup` must be smaller than `ndraws` for an initial run.")
  }
  if (is.null(prev_run) &&
      is_count1(ndraws, min = 2L) &&
      is_count1(save_every, min = 1L) &&
      save_every > ndraws) {
    add_error("`save_every` must not exceed `ndraws` for an initial run.")
  }
  
  if (identical(sampler, "mh")) {
    if (is.null(prop_sd)) {
      add_error("`prop_sd` must be supplied for `sampler = \"mh\"`.")
    } else {
      check_positive_num("prop_sd", prop_sd)
    }
  } else if (!is.null(prop_sd)) {
    check_positive_num("prop_sd", prop_sd)
  }

  check_positive_num("slice_width", slice_width)

  if (!is_num1(ini_spread, finite = TRUE) || ini_spread < 0 || ini_spread > 1) {
    add_error("`ini_spread` must be a single finite number between 0 and 1.")
  }
  
  if (!is_count_or_inf1(maxit, min = 1L)) {
    add_error("`maxit` must be a single positive whole number or Inf.")
  }

  if (!is_num1(max_rhat, finite = TRUE) || max_rhat < 1) {
    add_error("`max_rhat` must be a single finite number >= 1.")
  }
  
  if (!is.null(min_effss)) {
    check_count("min_effss", min_effss, min = 1L)
  }
  
  if (isTRUE(update_till_converge)) {
    if (!is_count1(chains, min = 2L)) {
      add_error("`update_till_converge = TRUE` requires `chains >= 2` for Gelman-Rubin diagnostics.")
    }
    if (is.null(min_effss) || !is_count1(min_effss, min = 1L)) {
      add_error("`min_effss` must be supplied as a positive whole number when `update_till_converge = TRUE`.")
    }
  }
  
  check_positive_num("tau_t", tau_t)
  check_positive_num("sig_prior", sig_prior)
  check_positive_num("tau_g", tau_g)
  check_positive_num("q_prior_sd", q_prior_sd)
  
  ## Kappa
  if (is_logical1(update_kappa)) {
    if (!update_kappa) {
      if (is.null(kappa)) {
        add_error(
          "`kappa` must be specified when `update_kappa = FALSE` because test sensitivity is not estimated."
        )
      } else if (!is_num1(kappa, finite = TRUE) || !(kappa > 0 && kappa <= 1)) {
        add_error(
          "`kappa` must be a single number in (0, 1] when `update_kappa = FALSE`."
        )
      }
    } else if (is.null(prev_run)) {
      if (!is.null(kappa)) {
        warning(
          "`kappa` is ignored when `update_kappa = TRUE`.",
          call. = FALSE
        )
      }
      if (is.null(kappa_prior)) {
        warning(
          paste0(
            "`kappa_prior` was not specified. Kappa will be estimated using an ",
            "uninformative Uniform(0, 1) prior. An informative Beta prior can be ",
            "specified through `kappa_prior = c(mean, sd)`."
          ),
          call. = FALSE
        )
      }
    }
  }
  
  if (isTRUE(update_kappa) && !is.null(kappa_prior)) {
    if (!(is.numeric(kappa_prior) && length(kappa_prior) == 2L &&
          all(is.finite(kappa_prior)) && !anyNA(kappa_prior))) {
      add_error("`kappa_prior` must be NULL or a numeric vector of length 2: c(mean, sd).")
    } else {
      kappa_mean <- kappa_prior[1L]
      kappa_sd <- kappa_prior[2L]
      
      if (!(kappa_mean > 0 && kappa_mean < 1)) {
        add_error("`kappa_prior[1]`, the prior mean for kappa, must lie strictly between 0 and 1.")
      }
      
      if (!(kappa_sd > 0)) {
        add_error("`kappa_prior[2]`, the prior standard deviation for kappa, must be positive.")
      }
      
      if (kappa_mean > 0 && kappa_mean < 1 && kappa_sd > 0) {
        max_sd <- sqrt(kappa_mean * (1 - kappa_mean))
        if (kappa_sd >= max_sd) {
          add_error(sprintf(
            "`kappa_prior` is not feasible for a Beta prior: sd must be smaller than sqrt(mean * (1 - mean)) = %.4f.",
            max_sd
          ))
        }
      }
    }
  }
  
  ## v_obs
  n <- NA_integer_
  g_fixed <- NULL
  
  if (is.null(v_obs)) {
    add_error("`v_obs` must be supplied unless `prev_run` is supplied.")
  } else if (!is.list(v_obs)) {
    add_error("`v_obs` must be a list of numeric vectors, one vector per subject.")
  } else {
    n <- length(v_obs)
    
    if (n < 1L) {
      add_error("`v_obs` must contain at least one subject.")
    } else {
      bad_basic <- which(!vapply(
        v_obs,
        function(x) is.numeric(x) && length(x) >= 1L,
        logical(1)
      ))
      
      if (length(bad_basic) > 0L) {
        add_error(sprintf(
          "Each element of `v_obs` must be a numeric vector with at least one value; problem index/indices: %s.",
          fmt(bad_basic)
        ))
      }
      
      ok_idx <- setdiff(seq_along(v_obs), bad_basic)
      
      which_bad <- function(idx, pred) {
        if (length(idx) == 0L) {
          return(integer(0))
        }
        idx[vapply(v_obs[idx], pred, logical(1))]
      }
      
      bad_missing <- which_bad(ok_idx, function(x) anyNA(x))
      if (length(bad_missing) > 0L) {
        add_error(sprintf(
          "`v_obs` must not contain NA or NaN values; problem index/indices: %s.",
          fmt(bad_missing)
        ))
      }
      
      bad_first <- which_bad(ok_idx, function(x) x[1L] != 0)
      if (length(bad_first) > 0L) {
        add_error(sprintf(
          "Every vector in `v_obs` must start with 0; problem index/indices: %s.",
          fmt(bad_first)
        ))
      }
      
      bad_neg_inf <- which_bad(ok_idx, function(x) any(is.infinite(x) & x < 0))
      if (length(bad_neg_inf) > 0L) {
        add_error(sprintf(
          "`v_obs` may contain Inf only as a final right-censoring time; -Inf is not allowed. Problem index/indices: %s.",
          fmt(bad_neg_inf)
        ))
      }
      
      bad_inf_position <- which_bad(ok_idx, function(x) {
        y <- x[-length(x)]
        any(is.infinite(y))
      })
      if (length(bad_inf_position) > 0L) {
        add_error(sprintf(
          "`Inf` may only appear as the last value of a `v_obs` vector; problem index/indices: %s.",
          fmt(bad_inf_position)
        ))
      }
      
      bad_nonpositive_times <- which_bad(ok_idx, function(x) {
        y <- x[-1L]
        any(is.finite(y) & y <= 0)
      })
      if (length(bad_nonpositive_times) > 0L) {
        add_error(sprintf(
          "All screening/event times after the initial 0 in `v_obs` must be positive; problem index/indices: %s.",
          fmt(bad_nonpositive_times)
        ))
      }
      
      bad_order <- which_bad(ok_idx, function(x) {
        length(x) > 1L && any(diff(x) <= 0, na.rm = TRUE)
      })
      if (length(bad_order) > 0L) {
        add_error(sprintf(
          "Each `v_obs` vector must be strictly increasing; problem index/indices: %s.",
          fmt(bad_order)
        ))
      }
      
      g_fixed <- lengths(v_obs) == 1L
      
      if (is.null(prev_run) && all(g_fixed)) {
        add_error("All `v_obs` entries are `c(0)`, indicating baseline-positive/prevalent cases only. At least one non-prevalent screening history is needed for a new PIM run.")
      }
      
      if (identical(prev, FALSE) && any(g_fixed)) {
        add_error("`prev = FALSE` assumes no baseline prevalence, but one or more `v_obs` entries are `c(0)`, which encodes baseline prevalence.")
      }
    }
  }
  
  ## Design matrices
  check_design_matrix(x_t, "x_t", n)
  check_design_matrix(x_g, "x_g", n)
  
  ## r
  if (is.null(r)) {
    if (isTRUE(prev)) {
      add_error("`r` must be supplied when `prev = TRUE`.")
    }
  } else {
    r_type_ok <- is.numeric(r) || is.logical(r)
    
    if (!r_type_ok) {
      add_error("`r` must be a binary numeric, integer, or logical vector.")
    } else {
      if (!is.na(n) && length(r) != n) {
        add_error(sprintf(
          "`r` has length %s, but `v_obs` has length %s. These must match.",
          length(r), n
        ))
      }
      
      if (anyNA(r)) {
        add_error("`r` must not contain missing values.")
      }
      
      bad_r <- which(!(r %in% c(0, 1)))
      if (length(bad_r) > 0L) {
        add_error(sprintf(
          "`r` must contain only 0/1 or FALSE/TRUE values; problem index/indices: %s.",
          fmt(bad_r)
        ))
      }
      
      if (!is.null(g_fixed) && length(r) == length(g_fixed)) {
        bad_fixed_r <- which(g_fixed & as.numeric(r) != 1)
        if (length(bad_fixed_r) > 0L) {
          add_error(sprintf(
            "`v_obs[[i]] = c(0)` indicates a positive baseline test, so `r[i]` must be 1; problem index/indices: %s.",
            fmt(bad_fixed_r)
          ))
        }
      }
    }
  }
  
  stop_if_errors()
  invisible(TRUE)
}
