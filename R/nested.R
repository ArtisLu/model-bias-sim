# Nested (double) simulation, paper Section 2.7.2: every outer dataset is treated like the study
# data, i.e. its own refitted prediction model, inverse calibration, inner simulated studies, bias
# correction and non-parametric CI (and, with no effect, its own rejection limits).

#' Seeds for the inner simulated studies and the truth of one outer replicate
#' @param seed_aux The outer replicate's auxiliary seed (`draw_replicate_rng()`).
#' @param n_inner Number of inner replicates.
#' @return list(inner = seeds, truth = one seed).
nested_seeds <- function(seed_aux, n_inner) {
  s <- make_crn_seeds(seed_aux, n_inner + 1L)
  list(inner = s[seq_len(n_inner)], truth = s[n_inner + 1L])
}

#' Empty nested result row for a failed outer replicate
#' @param rep,seed Replicate index and seed; @param res Outer pipeline result.
#' @param status,message Failure status and message; @param cols Columns to fill with NA.
#' @return One-row data.table.
nested_fail_row <- function(rep, seed, res, status, message, cols) {
  out <- data.table::data.table(rep = rep, seed = seed, estimate = res$estimate, se = res$se)
  for (cl in cols) out[, (cl) := NA_real_]
  out[, `:=`(status = status, message = message)]
  out[]
}

nested_cols <- c("theta_hat_r", "theta_r", "n_eval", "at_bound", "mean_inner", "mean_inner_mcse",
                 "n_inner_ok", "truth_r", "truth_r_mcse", "bias_r", "est_adj", "q_lo", "q_hi",
                 "ci_lo_np", "ci_hi_np", "ci_lo_np_refl", "ci_hi_np_refl", "seconds", "seconds_cal")

#' One outer replicate at the calibrated effect (paper 2.7.2, nested steps 2-8)
#' @param data Base data; @param dgm Fitted DGM; @param calibration From `calibrate_observed()`.
#' @param rep,seed Outer replicate index and seed; @param spec Dataset entry of `cfg$datasets`.
#' @param cfg Full config (uses `nested`, `calibration$bound`); @param level Nominal level.
#' @param keep_inner Keep the inner estimates at theta_r in list column `inner_est` (pilot only).
#' @return One-row data.table (see `nested_cols`); `seconds_cal` is the time of the inner
#'   simulations (calibration steps), `seconds` the whole outer replicate.
nested_outer_one <- function(data, dgm, calibration, rep, seed, spec, cfg, level = 0.95,
                             keep_inner = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  nc <- cfg$nested
  rng <- draw_replicate_rng(seed, nrow(data))
  sim <- simulate_outcomes(data, dgm, calibration$theta_cal, rng$u_event, rng$u_cens)
  res <- safely_run(run_pipeline, sim, spec$pipeline)
  if (is.na(res$estimate)) {
    return(nested_fail_row(rep, seed, res, "error", res$message, nested_cols))
  }
  out <- tryCatch({
    # Refit the prediction model to the outer dataset, as for the study data
    dgm_r <- fit_dgm(sim, spec$dgm, spec$censoring)
    sd <- nested_seeds(rng$seed_aux, nc$inner)
    scale <- if (is.null(spec$tol_scale)) 1 else spec$tol_scale
    cc <- list(max_iter = nc$max_iter, tol = nc$tol * scale, bound = cfg$calibration$bound)
    t_cal <- proc.time()[["elapsed"]]
    cal <- calibrate_mean(dgm_r, data, res$estimate, sd$inner, spec$pipeline, cc,
                          start = calibration$theta_cal, slope = calibration$slope)
    seconds_cal <- proc.time()[["elapsed"]] - t_cal
    tr <- compute_truth(dgm_r, data, cal$theta_cal, spec$pipeline, cfg, n_rep = nc$truth_n_rep,
                        n_target = nc$truth_n_target, seed = sd$truth)
    est_in <- cal$estimates[!is.na(cal$estimates)]
    bias_r <- mean(est_in) - tr$truth
    est_adj <- res$estimate - bias_r
    q <- centred_quantiles(est_in, level)
    data.table::data.table(
      rep = rep, seed = seed, estimate = res$estimate, se = res$se,
      theta_hat_r = dgm_r$theta_hat, theta_r = cal$theta_cal, n_eval = cal$n_eval,
      at_bound = as.numeric(cal$at_bound), mean_inner = mean(est_in),
      mean_inner_mcse = stats::sd(est_in) / sqrt(length(est_in)), n_inner_ok = length(est_in),
      truth_r = tr$truth, truth_r_mcse = tr$truth_mcse, bias_r = bias_r, est_adj = est_adj,
      q_lo = q[["q_lo"]], q_hi = q[["q_hi"]],
      ci_lo_np = est_adj + q[["q_lo"]], ci_hi_np = est_adj + q[["q_hi"]],
      ci_lo_np_refl = est_adj - q[["q_hi"]], ci_hi_np_refl = est_adj - q[["q_lo"]],
      seconds = NA_real_, seconds_cal = seconds_cal, status = res$status, message = res$message
    )
  }, error = function(e) {
    nested_fail_row(rep, seed, res, "error_nested", conditionMessage(e), nested_cols)
  })
  if (keep_inner && !is.na(out$est_adj)) out[, inner_est := list(list(est_in))]
  out[, seconds := proc.time()[["elapsed"]] - t0]
  out[]
}

#' One outer replicate with no effect: refitted model, inner null studies, own rejection limits
#' @param data Base data; @param dgm Fitted DGM; @param rep,seed Outer replicate index and seed.
#' @param spec Dataset entry of `cfg$datasets`; @param cfg Full config; @param level Nominal level.
#' @return One-row data.table: estimate, inner null mean and limits, reject, reject_fav (rejection
#'   in the direction favouring arm 1).
nested_null_one <- function(data, dgm, rep, seed, spec, cfg, level = 0.95) {
  t0 <- proc.time()[["elapsed"]]
  cols <- c("mean_inner", "n_inner_ok", "q_lo", "q_hi", "reject", "reject_fav", "seconds")
  low_fav <- endpoint_info(spec$pipeline)$favours_low
  rng <- draw_replicate_rng(seed, nrow(data))
  sim <- simulate_outcomes(data, dgm, 0, rng$u_event, rng$u_cens)
  res <- safely_run(run_pipeline, sim, spec$pipeline)
  if (is.na(res$estimate)) return(nested_fail_row(rep, seed, res, "error", res$message, cols))
  out <- tryCatch({
    dgm_r <- fit_dgm(sim, spec$dgm, spec$censoring)
    # No calibration under the null: conditional and marginal effects are both 0 (paper 2.7.3)
    seeds_in <- nested_seeds(rng$seed_aux, cfg$nested$inner)$inner
    est_in <- replicate_estimates(data, dgm_r, 0, seeds_in, spec$pipeline)
    est_in <- est_in[!is.na(est_in)]
    q <- centred_quantiles(est_in, level)
    d <- res$estimate - mean(est_in)
    data.table::data.table(
      rep = rep, seed = seed, estimate = res$estimate, se = res$se,
      mean_inner = mean(est_in), n_inner_ok = length(est_in),
      q_lo = q[["q_lo"]], q_hi = q[["q_hi"]],
      reject = as.numeric(d < q[["q_lo"]] | d > q[["q_hi"]]),
      reject_fav = as.numeric(if (low_fav) d < q[["q_lo"]] else d > q[["q_hi"]]),
      seconds = NA_real_, status = res$status, message = res$message
    )
  }, error = function(e) nested_fail_row(rep, seed, res, "error_nested", conditionMessage(e), cols))
  out[, seconds := proc.time()[["elapsed"]] - t0]
  out[]
}

#' Outer replicates of one batch (calibrated effect or no effect)
#' @param data Base data; @param dgm Fitted DGM; @param calibration From `calibrate_observed()`
#'   (NULL: no-effect scenario). @param seeds Outer seeds (the validation seeds).
#' @param batch_id Batch index; @param spec Dataset entry; @param cfg Full config.
#' @param keep_inner Keep the inner estimates (pilot only, see `nested_outer_one()`).
#' @return data.table, one row per outer replicate.
run_nested_batch <- function(data, dgm, calibration, seeds, batch_id, spec, cfg,
                             keep_inner = FALSE) {
  checkmate::assert_count(batch_id, positive = TRUE)
  nc <- cfg$nested
  checkmate::assert_true(nc$outer <= length(seeds))
  reps <- (batch_id - 1L) * nc$batch_size + seq_len(nc$batch_size)
  reps <- reps[reps <= nc$outer]
  level <- spec$pipeline$conf_level
  out <- lapply(reps, function(r) {
    if (is.null(calibration)) {
      nested_null_one(data, dgm, r, seeds[r], spec, cfg, level)
    } else {
      nested_outer_one(data, dgm, calibration, r, seeds[r], spec, cfg, level, keep_inner)
    }
  })
  data.table::rbindlist(out, fill = TRUE)
}

#' Summarise the nested simulation and derive the nested CI of the observed analysis
#' @param nested_obs,nested_null Outer replicates at the calibrated effect and with no effect.
#' @param truth_table From `make_truth_table()`; @param adj From `adjust_simple()`.
#' @param level Nominal level.
#' @param shrink Shrink the outer deviations of the adjusted estimate by sqrt(V / (V + sigma_w^2 / R))
#'   before taking the quantiles of the nested CI (paper Appendix); both CIs are reported.
#' @return list(summary = one-row data.table, long = rows as `perf_rows()`).
summarise_nested <- function(nested_obs, nested_null, truth_table, adj, level = 0.95,
                             shrink = FALSE) {
  tt <- truth_table[scenario == "observed"]
  o <- nested_obs[!is.na(est_adj)]
  n0 <- nested_null[!is.na(reject)]
  long <- rbind(
    perf_rows("Nested adjustment", "Non-parametric", o$est_adj, o$ci_lo_np, o$ci_hi_np, tt$truth,
              rej_null = n0$reject == 1, low_null = n0$reject_fav == 1,
              bias_extra_mcse = tt$truth_mcse),
    perf_rows("Nested adjustment", "Non-parametric (reflected)", o$est_adj, o$ci_lo_np_refl,
              o$ci_hi_np_refl, tt$truth, bias_extra_mcse = tt$truth_mcse)
  )
  # Same outer datasets with the correction of the simple adjustment (paired comparison)
  fx <- apply_adjustment(o$estimate, o$se, adj, level)
  fx0 <- apply_adjustment(n0$estimate, n0$se, adj, level)
  # Nested CI of the observed analysis: spread of the bias-adjusted estimator across outer
  # datasets, which includes the variation of the estimated bias correction
  q <- centred_quantiles(o$est_adj, level)
  comp <- nested_variance_components(o$est_adj, o$mean_inner_mcse^2 * o$n_inner_ok,
                                     stats::median(o$n_inner_ok))
  q_s <- centred_quantiles(shrink_deviations(o$est_adj, comp$shrink), level)
  q_use <- if (shrink && is.finite(comp$shrink)) q_s else q
  within <- sqrt(mean(o$mean_inner_mcse^2 + data.table::fifelse(is.na(o$truth_r_mcse), 0,
                                                                o$truth_r_mcse^2)))
  summary <- data.table::data.table(
    n_outer = nrow(nested_obs), n_outer_ok = nrow(o), n_fail = nrow(nested_obs) - nrow(o),
    n_null = nrow(nested_null), n_null_ok = nrow(n0), n_null_fail = nrow(nested_null) - nrow(n0),
    n_inner = max(c(o$n_inner_ok, n0$n_inner_ok)), n_at_bound = sum(o$at_bound == 1),
    mean_eval = mean(o$n_eval),
    seconds = sum(nested_obs$seconds, na.rm = TRUE) + sum(nested_null$seconds, na.rm = TRUE),
    theta_r_mean = mean(o$theta_r), theta_r_sd = stats::sd(o$theta_r),
    bias_r_mean = mean(o$bias_r), bias_r_sd = stats::sd(o$bias_r), bias_r_within_mcse = within,
    var_obs = comp$var_obs, sigma_w2 = comp$sigma_w2, var_V = comp$V, k_ratio = comp$k,
    shrink_factor = comp$shrink, shrunk = shrink && is.finite(comp$shrink),
    q_lo = q_use[["q_lo"]], q_hi = q_use[["q_hi"]],
    ci_lo_np = adj$est_adj + q_use[["q_lo"]], ci_hi_np = adj$est_adj + q_use[["q_hi"]],
    ci_lo_np_unshrunk = adj$est_adj + q[["q_lo"]], ci_hi_np_unshrunk = adj$est_adj + q[["q_hi"]],
    ci_lo_np_refl = adj$est_adj - q_use[["q_hi"]], ci_hi_np_refl = adj$est_adj - q_use[["q_lo"]],
    coverage_fixed = mean(fx$ci_lo_np <= tt$truth & tt$truth <= fx$ci_hi_np),
    type1_fixed = mean(fx0$rej_np)
  )
  list(summary = summary, long = long)
}
