library(BayesPIM)

set.seed(771)
storage_data <- gen_data(
  n = 24,
  p = 1,
  beta_t = 0.1,
  beta_g = 0.1,
  theta = 0.15,
  mu_t = 1,
  sigma_t = 0.4,
  v_min = 0.5,
  v_max = 1,
  mean_rc = 4,
  prob_r = 1
)

base_args <- list(
  v_obs = storage_data$v_obs,
  x_t = storage_data$x,
  x_g = storage_data$x,
  r = storage_data$r,
  kappa = 0.8,
  ndraws = 8,
  warmup = 3,
  chains = 1,
  seed_chains = 817,
  min_effss = 1,
  silent = TRUE
)

# Storage selection must not change the Markov chain trajectory. Check this for
# every incidence sampler by comparing a sparse fit with the corresponding rows
# from a fit that stores every draw.
for (sampler in c("slice_collapsed", "slice", "mh")) {
  sampler_args <- c(
    base_args,
    list(sampler = sampler),
    if (sampler == "mh") list(prop_sd = 0.01) else list()
  )
  full_fit <- do.call(bayespim, sampler_args)
  sparse_fit <- do.call(
    bayespim,
    c(sampler_args, list(save_every = 2))
  )

  stopifnot(
    nrow(as.matrix(full_fit$par[[1L]])) == 8L,
    nrow(as.matrix(sparse_fit$par[[1L]])) == 4L,
    object.size(sparse_fit$par) < object.size(full_fit$par),
    identical(
      as.matrix(sparse_fit$par[[1L]]),
      as.matrix(full_fit$par[[1L]])[c(2L, 4L, 6L, 8L), , drop = FALSE]
    ),
    identical(sparse_fit$terminal_par, full_fit$terminal_par),
    identical(sparse_fit$rng_state, full_fit$rng_state),
    sparse_fit$total_iterations == 8L,
    sparse_fit$save_every == 2L,
    coda::thin(sparse_fit$par) == 2
  )
}

# Estimated sensitivity is part of both the stored and terminal parameter
# state, so storage thinning must preserve its trajectory as well.
kappa_args <- base_args
kappa_args$kappa <- NULL
kappa_args$update_kappa <- TRUE
kappa_args$kappa_prior <- c(0.8, 0.05)
kappa_args$seed_chains <- 818
full_kappa_fit <- do.call(bayespim, kappa_args)
sparse_kappa_fit <- do.call(
  bayespim,
  c(kappa_args, list(save_every = 2))
)
stopifnot(
  identical(
    as.matrix(sparse_kappa_fit$par[[1L]]),
    as.matrix(full_kappa_fit$par[[1L]])[c(2L, 4L, 6L, 8L), , drop = FALSE]
  ),
  identical(sparse_kappa_fit$terminal_par, full_kappa_fit$terminal_par),
  identical(colnames(sparse_kappa_fit$terminal_par), colnames(sparse_kappa_fit$par[[1L]]))
)

# A continuation must start at the true terminal state even if the most recent
# generated iteration was not stored. Global storage indices must remain
# regular across update boundaries, including an update that stores no rows.
initial_args <- base_args
initial_args$ndraws <- 5
initial_args$save_every <- 5
initial_fit <- do.call(bayespim, initial_args)
no_new_stored_draw <- bayespim(
  prev_run = initial_fit,
  ndraws_update = 2,
  min_effss = 1,
  silent = TRUE
)
continued_fit <- bayespim(
  prev_run = no_new_stored_draw,
  ndraws_update = 3,
  min_effss = 1,
  silent = TRUE
)

uninterrupted_args <- base_args
uninterrupted_args$ndraws <- 10
uninterrupted_args$save_every <- 5
uninterrupted_fit <- do.call(bayespim, uninterrupted_args)

stopifnot(
  nrow(as.matrix(no_new_stored_draw$par[[1L]])) == 1L,
  no_new_stored_draw$total_iterations == 7L,
  identical(
    as.matrix(continued_fit$par[[1L]]),
    as.matrix(uninterrupted_fit$par[[1L]])
  ),
  identical(continued_fit$terminal_par, uninterrupted_fit$terminal_par),
  identical(continued_fit$rng_state, uninterrupted_fit$rng_state),
  continued_fit$total_iterations == 10L,
  identical(as.numeric(stats::time(continued_fit$par[[1L]])), c(5, 10))
)

# Analysis helpers use every stored post-warm-up draw and never thin again.
analysis_chains <- BayesPIM:::bayespim_analysis_chains(sparse_fit)
stopifnot(
  identical(as.numeric(stats::time(analysis_chains[[1L]])), c(4, 6, 8)),
  sparse_fit$convergence$n_iter == 8L,
  sparse_fit$convergence$n_iter_diagnostic == 3L
)

sparse_summary <- NULL
summary_output <- capture.output(sparse_summary <- summary(sparse_fit))
stopifnot(
  sparse_summary$draws$total_generated == 8L,
  sparse_summary$draws$total_stored == 4L,
  sparse_summary$draws$warmup_per_chain == 3L,
  sparse_summary$draws$stored_warmup == 1L,
  sparse_summary$draws$stored_warmup_per_chain == 1L,
  any(grepl(
    "Warm-up cutoff: 3 generated iterations per chain",
    summary_output,
    fixed = TRUE
  )),
  any(grepl(
    "Stored warm-up draws omitted: 1 (1 per chain)",
    summary_output,
    fixed = TRUE
  )),
  sparse_summary$draws$retained_stored == 3L,
  sparse_summary$draws$retained_stored_per_chain == 3L,
  any(grepl(
    "Posterior draws used: 3 (3 per chain)",
    summary_output,
    fixed = TRUE
  ))
)

# `warmup` is always a raw generated-iteration cutoff, whether or not it is
# divisible by `save_every`. With raw stored indices 10, 20, ..., 10000, both
# cutoffs below omit 50 stored draws per chain and retain raw iteration 510
# onward; neither cutoff is multiplied by 10.
indexed_chain <- coda::mcmc(
  matrix(seq_len(1000L), ncol = 1L),
  start = 10L,
  end = 10000L,
  thin = 10L
)
indexed_chains <- coda::mcmc.list(
  replicate(4L, indexed_chain, simplify = FALSE)
)
indexed_fit <- structure(
  list(
    par = indexed_chains,
    total_iterations = 10000L,
    save_every = 10L
  ),
  class = c("bayespim", "list")
)
for (warmup_cutoff in c(500L, 505L)) {
  indexed_retained <- BayesPIM:::bayespim_analysis_chains(
    indexed_fit,
    warmup = warmup_cutoff
  )
  indexed_counts <- BayesPIM:::bayespim_summary_draw_counts(
    object = indexed_fit,
    retained_chains = indexed_retained,
    warmup = warmup_cutoff
  )
  stopifnot(
    stats::start(indexed_retained[[1L]]) == 510L,
    indexed_counts$warmup_per_chain == warmup_cutoff,
    indexed_counts$stored_warmup == 200L,
    indexed_counts$stored_warmup_per_chain == 50L,
    indexed_counts$retained_stored == 3800L,
    indexed_counts$retained_stored_per_chain == 950L
  )
}

trimmed_sparse <- trim_mcmc(sparse_fit$par, burnin = 1, thinning = 2)
stopifnot(
  identical(as.numeric(stats::time(trimmed_sparse[[1L]])), c(4, 8)),
  coda::thin(trimmed_sparse) == 4
)

plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
plot_return <- plot(sparse_fit)
grDevices::dev.off()
unlink(plot_file)
stopifnot(identical(plot_return, sparse_fit))

too_many_ic_draws <- try(get_ic(sparse_fit, samples = 4, cores = 1), silent = TRUE)
too_many_ppcif_draws <- try(
  ppCIF(sparse_fit, pst_samples = 4, quant = c(0, 1)),
  silent = TRUE
)
stopifnot(
  inherits(too_many_ic_draws, "try-error"),
  inherits(too_many_ppcif_draws, "try-error"),
  grepl(
    "cannot exceed the number of stored post-warm-up posterior draws (3)",
    as.character(too_many_ppcif_draws),
    fixed = TRUE
  )
)

invalid_save_every <- try(
  do.call(bayespim, c(base_args, list(save_every = 9))),
  silent = TRUE
)
stopifnot(
  inherits(invalid_save_every, "try-error"),
  grepl(
    "`save_every` must not exceed `ndraws` for an initial run.",
    as.character(invalid_save_every),
    fixed = TRUE
  )
)
