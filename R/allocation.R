# Split of the nested simulation into outer replicates (B) and inner runs (R): ANOVA decomposition,
# optimal R from the pilot, shrinkage of the outer deviations (paper Appendix).

#' ANOVA decomposition of the spread of an outer-level estimate
#' Var_obs = V + sigma_w^2 / R: V is the variance of the estimator without inner noise, sigma_w^2
#' the variance of single inner estimates (pooled within outer replicates).
#' @param x Outer-level estimates (one per outer replicate), e.g. the bias-adjusted estimates.
#' @param sigma_w2_i Within-replicate variance of single inner estimates, one per outer replicate.
#' @param R Inner runs behind each outer estimate.
#' @return list(var_obs, sigma_w2, R, V, k = sigma_w2 / V, shrink = sqrt(V / var_obs)).
nested_variance_components <- function(x, sigma_w2_i, R) {
  x <- x[is.finite(x)]
  var_obs <- stats::var(x)
  sigma_w2 <- mean(sigma_w2_i[is.finite(sigma_w2_i)])
  v <- var_obs - sigma_w2 / R
  list(
    var_obs = var_obs, sigma_w2 = sigma_w2, R = R, V = v,
    k = if (is.finite(v) && v > 0) sigma_w2 / v else NA_real_,
    # Inner noise inflates the spread; a non-positive V means R is too small to separate them
    shrink = if (is.finite(v) && v > 0) sqrt(v / var_obs) else NA_real_
  )
}

#' Shrink deviations from the mean by a factor (percentile intervals, paper Appendix)
#' @param x Estimates; @param factor Shrink factor in (0, 1]; NA leaves x unchanged.
#' @return Numeric vector.
shrink_deviations <- function(x, factor) {
  if (!is.finite(factor)) return(x)
  mean(x) + factor * (x - mean(x))
}

#' Optimal number of inner runs for the variance target (paper Appendix)
#' R* = (k / 2) * (1 + sqrt(1 + 8 c_o / (k c_i))), approximately sqrt(2 k c_o / c_i) for c_o >> c_i.
#' @param k sigma_w^2 / V; @param c_o Cost of an outer replicate without inner runs.
#' @param c_i Cost of one inner run (per unit of R, incl. all calibration steps).
#' @return Numeric R* (not rounded).
optimal_inner <- function(k, c_o, c_i) {
  checkmate::assert_number(k, lower = 0)
  checkmate::assert_number(c_o, lower = 0)
  checkmate::assert_number(c_i, lower = 1e-12)
  (k / 2) * (1 + sqrt(1 + 8 * c_o / (k * c_i)))
}

#' Outer replicates needed for a relative precision of the empirical SE
#' MCSE(SE) = SE / sqrt(2 (B - 1)), so B = 1 + 1 / (2 rel^2).
#' @param rel Target relative MCSE (0.02 gives about 1,250).
#' @return Integer B.
outer_for_precision <- function(rel) {
  checkmate::assert_number(rel, lower = 1e-4)
  as.integer(ceiling(1 + 1 / (2 * rel^2)))
}

#' Allocation for percentile limits: R proportional to C^(1/3), B to C^(2/3) (Gordy and Juneja
#' 2010; the exact statement should be checked before citing)
#' @param budget Computing budget; @param c_o,c_i Costs; @param ref list(R, B) of a reference design
#'   (the variance-target design) whose budget is `ref$budget`.
#' @return list(R, B, cost).
scale_percentile_allocation <- function(budget, c_o, c_i, ref) {
  R <- max(2, round(ref$R * (budget / ref$budget)^(1 / 3)))
  B <- floor(budget / (c_o + R * c_i))
  list(R = R, B = B, cost = B * (c_o + R * c_i))
}

#' Time of one inner run: simulate a study at theta and apply the pipeline
#' @param data Base data; @param dgm Fitted DGM; @param theta Generating value.
#' @param pipeline Pipeline spec; @param n Number of studies timed.
#' @return Seconds per study (elapsed).
time_inner_run <- function(data, dgm, theta, pipeline, n = 20) {
  seeds <- make_crn_seeds(1L, n)
  t0 <- proc.time()[["elapsed"]]
  invisible(replicate_estimates(data, dgm, theta, seeds, pipeline))
  (proc.time()[["elapsed"]] - t0) / n
}

#' ANOVA check on sub-samples of the inner runs: V = Var_obs(r) - sigma_w^2 / r should not depend on r
#' @param pilot Pilot outer replicates with list column `inner_est` (estimates at theta_r).
#' @param r_grid Numbers of inner runs to use (first r of the R stored).
#' @return data.table: r, var_obs, sigma_w2, V (of the bias-adjusted estimate).
anova_by_inner_runs <- function(pilot, r_grid) {
  p <- pilot[!is.na(est_adj) & lengths(inner_est) > 0]
  R <- min(lengths(p$inner_est))
  data.table::rbindlist(lapply(r_grid[r_grid <= R], function(r) {
    m_r <- vapply(p$inner_est, function(e) mean(e[seq_len(r)]), numeric(1))
    v_r <- vapply(p$inner_est, function(e) stats::var(e[seq_len(r)]), numeric(1))
    # Bias-adjusted estimate recomputed with r inner runs: estimate - (mean_r - truth_r)
    x <- p$estimate - (m_r - p$truth_r)
    cp <- nested_variance_components(x, v_r, r)
    data.table::data.table(r = r, var_obs = cp$var_obs, sigma_w2 = cp$sigma_w2, V = cp$V)
  }))
}

#' Width of the per-outer-dataset limits from the inner studies, by number of inner runs
#' The centred 2.5% and 97.5% quantiles of r inner estimates are the limits of the outer CI; with few
#' inner runs they are too narrow (the range of r draws covers a new draw with probability (r-1)/(r+1)).
#' @param pilot Pilot outer replicates with list column `inner_est`; @param r_grid Inner runs to use.
#' @param level Nominal level.
#' @return data.table: r, width (mean over outer replicates), width_rel (relative to the largest r).
limit_width_by_inner_runs <- function(pilot, r_grid, level = 0.95) {
  p <- pilot[!is.na(est_adj) & lengths(inner_est) > 0]
  R <- min(lengths(p$inner_est))
  out <- data.table::rbindlist(lapply(r_grid[r_grid <= R], function(r) {
    w <- vapply(p$inner_est, function(e) {
      q <- centred_quantiles(e[seq_len(r)], level)
      q[["q_hi"]] - q[["q_lo"]]
    }, numeric(1))
    data.table::data.table(r = r, width = mean(w))
  }))
  out[, width_rel := width / width[which.max(r)]]
  out[]
}

#' Recommended split from the pilot study (paper Appendix, practical recipe steps 1-3)
#' R* minimises the variance of V-hat per unit cost. The limits of the outer CIs and tests are
#' quantiles of the inner studies and are too narrow with few of them, so the recommended R is at
#' least `pc$min_inner_limits` (if set).
#' @param pilot Pilot outer replicates at theta* (as `nested_obs`); @param pc Pilot settings
#'   (`cfg$pilot`); @param current The production design list(outer, inner) for comparison.
#' @param c_study Directly timed seconds per inner study (cross-check of the cost split).
#' @return One-row data.table: variance components, costs, R*, B, cost of the recommended and of
#'   the current design (CPU seconds) and the relative variance of V-hat of both at equal cost.
recommend_allocation <- function(pilot, pc, current, c_study = NA_real_) {
  p <- pilot[!is.na(est_adj)]
  R0 <- stats::median(p$n_inner_ok)
  cp <- nested_variance_components(p$est_adj, p$mean_inner_mcse^2 * p$n_inner_ok, R0)
  n_eval <- mean(p$n_eval)
  # Cost split measured in the pilot: calibration steps are the inner runs, the rest is the outer cost
  c_i <- mean(p$seconds_cal) / R0
  c_o <- max(mean(p$seconds - p$seconds_cal), 1e-6)
  k_use <- if (is.finite(cp$k)) cp$k else 1
  r_star <- optimal_inner(k_use, c_o, c_i)
  r_var <- max(2L, as.integer(round(r_star)))
  r_rec <- if (is.null(pc$min_inner_limits)) r_var else max(r_var, as.integer(pc$min_inner_limits))
  b_prec <- outer_for_precision(pc$precision_rel_se)
  b_rec <- max(b_prec, pc$min_outer_percentile)
  cost <- function(B, R) B * (c_o + R * c_i)
  # Variance of V-hat (units of V^2) for B outer and R inner runs: B Var = 2 (1 + k/R)^2 + 2 k^2 / (R^2 (R - 1))
  vv <- function(B, R) (2 * (1 + k_use / R)^2 + 2 * k_use^2 / (R^2 * (R - 1))) / B
  data.table::data.table(
    n_outer_pilot = nrow(p), R_pilot = R0, var_obs = cp$var_obs, sigma_w2 = cp$sigma_w2, V = cp$V,
    k = cp$k, shrink_factor = cp$shrink, mean_eval = n_eval,
    c_o_sec = c_o, c_i_sec = c_i, c_study_direct_sec = c_study, cost_ratio = c_o / c_i,
    R_star = r_star, R_variance = r_var, R_rec = r_rec, B_precision = b_prec, B_rec = b_rec,
    cost_rec_cpu_h = cost(b_rec, r_rec) / 3600,
    cost_rstar_cpu_h = cost(b_rec, r_var) / 3600,
    relvar_rstar = vv(b_rec, r_var), relvar_rec = vv(b_rec, r_rec),
    cost_current_cpu_h = cost(current$outer, current$inner) / 3600,
    # Same cost as the current design: B that the recommended R would allow, and its precision gain
    B_at_current_cost = floor(cost(current$outer, current$inner) / (c_o + r_rec * c_i)),
    relvar_current = vv(current$outer, current$inner),
    relvar_rec_same_cost = vv(floor(cost(current$outer, current$inner) / (c_o + r_rec * c_i)), r_rec)
  )
}

#' Pilot study: outer replicates at theta* with few outer and inner runs, plus timing
#' @param data Base data; @param dgm Fitted DGM; @param calibration From `calibrate_observed()`.
#' @param seeds Outer seeds; @param spec Dataset entry; @param cfg Full config.
#' @param batch_id Batch index (batches of `cfg$pilot$batch_size` outer replicates).
#' @return data.table of outer replicates with the inner estimates in `inner_est`.
run_pilot_batch <- function(data, dgm, calibration, seeds, batch_id, spec, cfg) {
  cfg_p <- cfg
  cfg_p$nested$outer <- cfg$pilot$outer
  cfg_p$nested$inner <- cfg$pilot$inner
  cfg_p$nested$batch_size <- cfg$pilot$batch_size
  run_nested_batch(data, dgm, calibration, seeds, batch_id, spec, cfg_p, keep_inner = TRUE)
}
