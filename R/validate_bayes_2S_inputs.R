#' Input validation for Bayes.2S, internal use
#' @noRd
.validate_bayes_2S_inputs <- function(
    Vobs = NULL,
    Z.X = NULL,
    Z.W = NULL,
    r = NULL,
    dist.X = "weibull",
    kappa = 0.5,
    update.kappa = FALSE,
    kappa.prior = NULL,
    ndraws = 1000,
    prop.sd.X = NULL,
    chains = NULL,
    thining = 1,
    parallel = TRUE,
    update.till.converge = FALSE,
    maxit = Inf,
    conv.crit = "upper",
    min_effss = NULL,
    beta.prior = "norm",
    beta.prior.X = 1,
    sig.prior.X = 1,
    tau.w = 1,
    fix.sigma.X = FALSE,
    prev.run = NULL,
    update.burnin = TRUE,
    ndraws.update = NULL,
    prev = TRUE,
    vanilla = FALSE,
    ndraws.naive = 1e4,
    naive.run.prop.sd.X = prop.sd.X,
    par.exp = FALSE,
    collapsed.g = TRUE,
    k.prior = 1,
    fix.k = FALSE,
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
          "Invalid input to `bayes.2S()`:\n",
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
  
  check_design_matrix <- function(Z, name, n) {
    if (is.null(Z)) {
      return(invisible(TRUE))
    }
    
    if (is.data.frame(Z)) {
      non_numeric <- names(Z)[!vapply(Z, is.numeric, logical(1))]
      if (length(non_numeric) > 0L) {
        add_error(sprintf(
          "`%s` must contain only numeric columns; non-numeric column(s): %s.",
          name, fmt(non_numeric)
        ))
        return(invisible(FALSE))
      }
    }
    
    Zmat <- tryCatch(
      as.matrix(Z),
      error = function(e) e
    )
    
    if (inherits(Zmat, "error")) {
      add_error(sprintf("`%s` could not be converted to a matrix.", name))
      return(invisible(FALSE))
    }
    
    if (!is.numeric(Zmat)) {
      add_error(sprintf("`%s` must be numeric.", name))
      return(invisible(FALSE))
    }
    
    if (!is.na(n) && nrow(Zmat) != n) {
      add_error(sprintf(
        "`%s` has %s row(s), but `Vobs` has length %s. These must match.",
        name, nrow(Zmat), n
      ))
    }
    
    if (ncol(Zmat) < 1L) {
      add_error(sprintf(
        "`%s` has zero columns. Use NULL for an intercept-only model.",
        name
      ))
    }
    
    if (anyNA(Zmat)) {
      add_error(sprintf("`%s` contains missing values; missing covariates are not allowed.", name))
    }
    
    if (any(!is.finite(Zmat))) {
      add_error(sprintf("`%s` must contain only finite numeric values.", name))
    }
    
    invisible(TRUE)
  }
  
  ## Validate previous run before bayes.2S dereferences it
  if (!is.null(prev.run)) {
    if (!is.list(prev.run)) {
      add_error("`prev.run` must be a previous `bayes.2S()` result object.")
    } else {
      required_prev <- c(
        "par.X.all", "Z.X", "Z.W", "r", "Vobs", "kappa",
        "update.kappa", "prev", "vanilla", "par.exp",
        "collapsed.g", "ndraws", "prop.sd.X", "priors",
        "thining", "X", "dist.X", "fix.sigma.X",
        "beta.prior", "runtime", "fix.k", "maxit"
      )
      
      missing_prev <- setdiff(required_prev, names(prev.run))
      if (length(missing_prev) > 0L) {
        add_error(sprintf(
          "`prev.run` is missing required component(s): %s.",
          fmt(paste0("$", missing_prev))
        ))
      }
      
      if ("priors" %in% names(prev.run)) {
        if (!is.list(prev.run$priors)) {
          add_error("`prev.run$priors` must be a list.")
        } else {
          missing_priors <- setdiff(
            c("beta.prior.X", "sig.prior.X"),
            names(prev.run$priors)
          )
          if (length(missing_priors) > 0L) {
            add_error(sprintf(
              "`prev.run$priors` is missing required component(s): %s.",
              fmt(paste0("$", missing_priors))
            ))
          }
        }
      }
      
      if ("par.X.all" %in% names(prev.run)) {
        if (!(inherits(prev.run$par.X.all, "mcmc.list") || is.list(prev.run$par.X.all))) {
          add_error("`prev.run$par.X.all` must be an `mcmc.list` or list-like MCMC object.")
        } else if (length(prev.run$par.X.all) < 1L) {
          add_error("`prev.run$par.X.all` contains no chains.")
        }
      }
    }
  }
  
  if (identical(stage, "initial")) {
    if (is.null(prev.run)) {
      if (is.null(Vobs)) {
        add_error("`Vobs` must be supplied unless `prev.run` is supplied.")
      }
      if (is.null(chains)) {
        add_error("`chains` must be supplied unless `prev.run` is supplied.")
      }
    }
    
    if (!is_num1(maxit, finite = FALSE, positive = TRUE)) {
      add_error("`maxit` must be a single positive number or Inf.")
    }
    
    if (!is.null(ndraws.update)) {
      check_count("ndraws.update", ndraws.update, min = 1L)
    }
    
    stop_if_errors()
    return(invisible(TRUE))
  }
  
  ## Logical scalars
  logical_args <- list(
    update.kappa = update.kappa,
    parallel = parallel,
    update.till.converge = update.till.converge,
    fix.sigma.X = fix.sigma.X,
    update.burnin = update.burnin,
    prev = prev,
    vanilla = vanilla,
    par.exp = par.exp,
    collapsed.g = collapsed.g,
    fix.k = fix.k
  )
  
  for (nm in names(logical_args)) {
    check_logical(nm, logical_args[[nm]])
  }
  
  ## Choices
  allowed_dist <- c("weibull", "lognormal", "loglog", "gengamma")
  if (!(is.character(dist.X) && length(dist.X) == 1L && !is.na(dist.X) && dist.X %in% allowed_dist)) {
    add_error(sprintf(
      "`dist.X` must be one of: %s.",
      fmt(allowed_dist)
    ))
  }
  
  allowed_beta_prior <- c("norm", "t")
  if (!(is.character(beta.prior) && length(beta.prior) == 1L && !is.na(beta.prior) &&
        beta.prior %in% allowed_beta_prior)) {
    add_error("`beta.prior` must be either 'norm' or 't'.")
  }
  
  if (!(is.character(conv.crit) && length(conv.crit) == 1L && !is.na(conv.crit) &&
        conv.crit %in% c("point", "upper"))) {
    add_error("`conv.crit` must be either 'point' or 'upper'.")
  }
  
  ## Counts and numeric tuning parameters
  check_count("ndraws", ndraws, min = 2L)
  check_count("chains", chains, min = 1L)
  check_count("thining", thining, min = 1L)
  check_count("ndraws.naive", ndraws.naive, min = 2L)
  
  if (!is.null(ndraws.update)) {
    check_count("ndraws.update", ndraws.update, min = 1L)
  }
  
  if (is.null(prop.sd.X)) {
    add_error("`prop.sd.X` must be supplied and must be a single positive finite number.")
  } else {
    check_positive_num("prop.sd.X", prop.sd.X)
  }
  
  if (is.null(naive.run.prop.sd.X)) {
    add_error("`naive.run.prop.sd.X` must be a single positive finite number. It defaults to `prop.sd.X`, so this usually means `prop.sd.X` was not supplied.")
  } else {
    check_positive_num("naive.run.prop.sd.X", naive.run.prop.sd.X)
  }
  
  if (!is_num1(maxit, finite = FALSE, positive = TRUE)) {
    add_error("`maxit` must be a single positive number or Inf.")
  }
  
  if (!is.null(min_effss)) {
    check_count("min_effss", min_effss, min = 1L)
  }
  
  if (isTRUE(update.till.converge)) {
    if (!is_count1(chains, min = 2L)) {
      add_error("`update.till.converge = TRUE` requires `chains >= 2` for Gelman-Rubin diagnostics.")
    }
    if (is.null(min_effss) || !is_count1(min_effss, min = 1L)) {
      add_error("`min_effss` must be supplied as a positive whole number when `update.till.converge = TRUE`.")
    }
  }
  
  check_positive_num("beta.prior.X", beta.prior.X)
  check_positive_num("sig.prior.X", sig.prior.X)
  check_positive_num("tau.w", tau.w)
  check_positive_num("k.prior", k.prior)
  
  ## Kappa
  if (!is_num1(kappa, finite = TRUE)) {
    add_error("`kappa` must be a single finite number.")
  } else if (is_logical1(update.kappa)) {
    if (update.kappa) {
      if (!(kappa > 0 && kappa < 1)) {
        add_error("When `update.kappa = TRUE`, `kappa` is a starting value and must lie strictly between 0 and 1.")
      }
    } else {
      if (!(kappa > 0 && kappa <= 1)) {
        add_error("When `update.kappa = FALSE`, `kappa` must lie in the interval (0, 1].")
      }
    }
  }
  
  if (isTRUE(update.kappa) && !is.null(kappa.prior)) {
    if (!(is.numeric(kappa.prior) && length(kappa.prior) == 2L &&
          all(is.finite(kappa.prior)) && !anyNA(kappa.prior))) {
      add_error("`kappa.prior` must be NULL or a numeric vector of length 2: c(mean, sd).")
    } else {
      kappa_mean <- kappa.prior[1L]
      kappa_sd <- kappa.prior[2L]
      
      if (!(kappa_mean > 0 && kappa_mean < 1)) {
        add_error("`kappa.prior[1]`, the prior mean for kappa, must lie strictly between 0 and 1.")
      }
      
      if (!(kappa_sd > 0)) {
        add_error("`kappa.prior[2]`, the prior standard deviation for kappa, must be positive.")
      }
      
      if (kappa_mean > 0 && kappa_mean < 1 && kappa_sd > 0) {
        max_sd <- sqrt(kappa_mean * (1 - kappa_mean))
        if (kappa_sd >= max_sd) {
          add_error(sprintf(
            "`kappa.prior` is not feasible for a Beta prior: sd must be smaller than sqrt(mean * (1 - mean)) = %.4f.",
            max_sd
          ))
        }
      }
    }
  }
  
  if (isTRUE(vanilla) && isTRUE(update.kappa)) {
    add_error("`vanilla = TRUE` is incompatible with `update.kappa = TRUE`; the vanilla model does not update kappa.")
  }
  
  ## Vobs
  n <- NA_integer_
  g.fixed <- NULL
  
  if (is.null(Vobs)) {
    add_error("`Vobs` must be supplied unless `prev.run` is supplied.")
  } else if (!is.list(Vobs)) {
    add_error("`Vobs` must be a list of numeric vectors, one vector per subject.")
  } else {
    n <- length(Vobs)
    
    if (n < 1L) {
      add_error("`Vobs` must contain at least one subject.")
    } else {
      bad_basic <- which(!vapply(
        Vobs,
        function(x) is.numeric(x) && length(x) >= 1L,
        logical(1)
      ))
      
      if (length(bad_basic) > 0L) {
        add_error(sprintf(
          "Each element of `Vobs` must be a numeric vector with at least one value; problem index/indices: %s.",
          fmt(bad_basic)
        ))
      }
      
      ok_idx <- setdiff(seq_along(Vobs), bad_basic)
      
      which_bad <- function(idx, pred) {
        if (length(idx) == 0L) {
          return(integer(0))
        }
        idx[vapply(Vobs[idx], pred, logical(1))]
      }
      
      bad_missing <- which_bad(ok_idx, function(x) anyNA(x))
      if (length(bad_missing) > 0L) {
        add_error(sprintf(
          "`Vobs` must not contain NA or NaN values; problem index/indices: %s.",
          fmt(bad_missing)
        ))
      }
      
      bad_first <- which_bad(ok_idx, function(x) x[1L] != 0)
      if (length(bad_first) > 0L) {
        add_error(sprintf(
          "Every vector in `Vobs` must start with 0; problem index/indices: %s.",
          fmt(bad_first)
        ))
      }
      
      bad_neg_inf <- which_bad(ok_idx, function(x) any(is.infinite(x) & x < 0))
      if (length(bad_neg_inf) > 0L) {
        add_error(sprintf(
          "`Vobs` may contain Inf only as a final right-censoring time; -Inf is not allowed. Problem index/indices: %s.",
          fmt(bad_neg_inf)
        ))
      }
      
      bad_inf_position <- which_bad(ok_idx, function(x) {
        y <- x[-length(x)]
        any(is.infinite(y))
      })
      if (length(bad_inf_position) > 0L) {
        add_error(sprintf(
          "`Inf` may only appear as the last value of a `Vobs` vector; problem index/indices: %s.",
          fmt(bad_inf_position)
        ))
      }
      
      bad_nonpositive_times <- which_bad(ok_idx, function(x) {
        y <- x[-1L]
        any(is.finite(y) & y <= 0)
      })
      if (length(bad_nonpositive_times) > 0L) {
        add_error(sprintf(
          "All screening/event times after the initial 0 in `Vobs` must be positive; problem index/indices: %s.",
          fmt(bad_nonpositive_times)
        ))
      }
      
      bad_order <- which_bad(ok_idx, function(x) {
        length(x) > 1L && any(diff(x) <= 0, na.rm = TRUE)
      })
      if (length(bad_order) > 0L) {
        add_error(sprintf(
          "Each `Vobs` vector must be strictly increasing; problem index/indices: %s.",
          fmt(bad_order)
        ))
      }
      
      g.fixed <- lengths(Vobs) == 1L
      
      if (!isTRUE(vanilla) && is.null(prev.run) && all(g.fixed)) {
        add_error("All `Vobs` entries are `c(0)`, indicating baseline-positive/prevalent cases only. At least one non-prevalent screening history is needed for a new PIM run.")
      }
      
      if (isTRUE(vanilla) && any(g.fixed)) {
        add_error("`Vobs` entries equal to `c(0)` encode baseline prevalence and are incompatible with `vanilla = TRUE`.")
      }
      
      if (!isTRUE(vanilla) && identical(prev, FALSE) && any(g.fixed)) {
        add_error("`prev = FALSE` assumes no baseline prevalence, but one or more `Vobs` entries are `c(0)`, which encodes baseline prevalence.")
      }
    }
  }
  
  ## Design matrices
  check_design_matrix(Z.X, "Z.X", n)
  check_design_matrix(Z.W, "Z.W", n)
  
  ## r
  if (is.null(r)) {
    if (!isTRUE(vanilla) && isTRUE(prev)) {
      add_error("`r` must be supplied when `prev = TRUE` and `vanilla = FALSE`.")
    }
  } else {
    r_type_ok <- is.numeric(r) || is.logical(r)
    
    if (!r_type_ok) {
      add_error("`r` must be a binary numeric, integer, or logical vector.")
    } else {
      if (!is.na(n) && length(r) != n) {
        add_error(sprintf(
          "`r` has length %s, but `Vobs` has length %s. These must match.",
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
      
      if (!is.null(g.fixed) && length(r) == length(g.fixed)) {
        bad_fixed_r <- which(g.fixed & as.numeric(r) != 1)
        if (length(bad_fixed_r) > 0L) {
          add_error(sprintf(
            "`Vobs[[i]] = c(0)` indicates a positive baseline test, so `r[i]` must be 1; problem index/indices: %s.",
            fmt(bad_fixed_r)
          ))
        }
      }
    }
  }
  
  stop_if_errors()
  invisible(TRUE)
}