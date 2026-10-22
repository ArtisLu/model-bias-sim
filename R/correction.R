# QBA adjustment of the observed analysis (paper Section 2.7): bias correction by inverse
# calibration, non-parametric (primary) and parametric CI and test, their check in an independent
# set of simulated studies, and long-format performance rows shared by all analyses.

#' Quantiles of the centred Monte Carlo errors (paper 2.7.2, step 2a)
#' @param x Simulated estimates; @param level Nominal level.
#' @return Named numeric: q_lo, q_hi of x - mean(x).
centred_quantiles <- function(x, level = 0.95) {
  x <- x[!is.na(x)]
  a <- (1 - level) / 2
  q <- stats::quantile(x - mean(x), c(a, 1 - a), names = FALSE, type = 8)
  c(q_lo = q[1], q_hi = q[2])
}

#' Simple QBA adjustment of the observed analysis (paper 2.7.1-2.7.3)
#' @param obs_result Observed pipeline result.
#' @param sim_results Main replicates at the calibrated theta and at theta = 0.
#' @param truth_table From `make_truth_table()`; @param dataset Dataset key.
#' @param level Nominal level.
#' @return One-row data.table on the log-HR scale.
adjust_simple <- function(obs_result, sim_results, truth_table, dataset, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  tt <- truth_table[scenario == "observed"]
  o <- sim_results[theta == tt$theta & !is.na(estimate)]
  n0 <- sim_results[theta == 0 & !is.na(estimate) & !is.na(se)]
  m <- mean(o$estimate)
  # Bias at the calibrated value; since m(theta_cal) = observed, psi(theta_cal) = observed - bias
  bias <- m - tt$truth
  est_adj <- obs_result$estimate - bias
  q <- centred_quantiles(o$estimate, level)
  emp_se <- stats::sd(o$estimate)
  k <- emp_se / sqrt(mean(o$se^2, na.rm = TRUE))
  # No effect: rejection limits from the centred null distribution (paper 2.7.3, steps 6-7)
  m0 <- mean(n0$estimate)
  e0 <- n0$estimate - m0
  q0 <- centred_quantiles(n0$estimate, level)
  d_obs <- obs_result$estimate - m0
  z0 <- abs(e0 / n0$se)
  z_obs <- abs(d_obs / obs_result$se)
  z_crit <- stats::quantile(z0, level, names = FALSE, type = 8)
  data.table::data.table(
    dataset = dataset, theta_cal = tt$theta,
    est_obs = obs_result$estimate, se_obs = obs_result$se,
    ci_lo_obs = obs_result$ci_lo, ci_hi_obs = obs_result$ci_hi, p_obs = obs_result$p_value,
    n_sim = nrow(o), mean_est = m, mean_est_mcse = emp_se / sqrt(nrow(o)),
    calib_resid = obs_result$estimate - m,
    truth = tt$truth, truth_mcse = tt$truth_mcse, bias = bias,
    # Truth is itself simulated: its MC error adds to that of the mean estimate
    bias_mcse = sqrt(emp_se^2 / nrow(o) + tt$truth_mcse^2),
    est_adj = est_adj,
    q_lo = q[["q_lo"]], q_hi = q[["q_hi"]],
    # Paper 2.7.2, step 5a: add the centred quantiles to the adjusted estimate
    ci_lo_np = est_adj + q[["q_lo"]], ci_hi_np = est_adj + q[["q_hi"]],
    # Reflected (basic bootstrap) form: matches the paper's form when the errors are symmetric
    ci_lo_np_refl = est_adj - q[["q_hi"]], ci_hi_np_refl = est_adj - q[["q_lo"]],
    emp_se = emp_se, mod_se = sqrt(mean(o$se^2, na.rm = TRUE)), se_inflation = k,
    ci_lo_par = est_adj - z * k * obs_result$se, ci_hi_par = est_adj + z * k * obs_result$se,
    n_null = nrow(n0), null_mean = m0, null_q_lo = q0[["q_lo"]], null_q_hi = q0[["q_hi"]],
    reject_np = d_obs < q0[["q_lo"]] | d_obs > q0[["q_hi"]],
    # Two-sided p from the centred null distribution (equal tails, as the rejection limits)
    p_np = min(1, 2 * min(mean(e0 <= d_obs), mean(e0 >= d_obs))),
    z_crit_par = z_crit, z_obs_centred = z_obs, reject_par = z_obs > z_crit,
    p_par = mean(z0 >= z_obs)
  )
}

#' Apply the adjustment learned by `adjust_simple()` to other estimates (validation, nested)
#' @param estimate,se Replicate estimates and model-based SEs; @param adj From `adjust_simple()`.
#' @param level Nominal level.
#' @return data.table: adjusted estimate, CIs (Wald, non-parametric, reflected, parametric) and
#'   rejections of H0 (overall and favouring arm 1) under each rule.
apply_adjustment <- function(estimate, se, adj, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  ea <- estimate - adj$bias
  d0 <- estimate - adj$null_mean
  data.table::data.table(
    est_adj = ea,
    ci_lo_wald = estimate - z * se, ci_hi_wald = estimate + z * se,
    ci_lo_np = ea + adj$q_lo, ci_hi_np = ea + adj$q_hi,
    ci_lo_np_refl = ea - adj$q_hi, ci_hi_np_refl = ea - adj$q_lo,
    ci_lo_par = ea - z * adj$se_inflation * se, ci_hi_par = ea + z * adj$se_inflation * se,
    rej_wald = abs(estimate / se) > z, low_wald = estimate / se < -z,
    rej_np = d0 < adj$null_q_lo | d0 > adj$null_q_hi, low_np = d0 < adj$null_q_lo,
    rej_par = abs(d0 / se) > adj$z_crit_par, low_par = d0 / se < -adj$z_crit_par
  )
}

#' Performance rows in long format (one analysis and CI/test variant)
#' @param analysis,method Labels. @param est Estimates (adjusted where applicable).
#' @param lo,hi CI limits (NULL: no coverage rows); @param truth True value.
#' @param rej_alt Rejections at the calibrated effect (power); NULL to skip.
#' @param rej_null,low_null Rejections and rejections favouring arm 1 with no effect.
#' @param bias_extra_mcse MC error of the reference the bias is measured against.
#' @return data.table: analysis, method, measure, value, mcse, ref, n.
perf_rows <- function(analysis, method, est, lo = NULL, hi = NULL, truth, rej_alt = NULL,
                      rej_null = NULL, low_null = NULL, bias_extra_mcse = 0) {
  prop <- function(measure, x, ref) {
    m <- mean(x)
    data.table::data.table(measure = measure, value = m, mcse = sqrt(m * (1 - m) / length(x)),
                           ref = ref, n = length(x))
  }
  n <- length(est)
  s <- stats::sd(est)
  out <- list(
    data.table::data.table(measure = "Mean estimate", value = mean(est), mcse = s / sqrt(n),
                           ref = truth, n = n),
    data.table::data.table(measure = "Bias (log HR)", value = mean(est) - truth,
                           mcse = sqrt(s^2 / n + bias_extra_mcse^2), ref = 0, n = n),
    data.table::data.table(measure = "Empirical SE", value = s, mcse = s / sqrt(2 * (n - 1)),
                           ref = NA_real_, n = n),
    data.table::data.table(measure = "RMSE", value = sqrt((mean(est) - truth)^2 + s^2),
                           mcse = NA_real_, ref = NA_real_, n = n)
  )
  if (!is.null(lo)) {
    out <- c(out, list(
      prop("95% CI coverage", lo <= truth & truth <= hi, 0.95),
      prop("CI entirely below truth", hi < truth, 0.025),
      prop("CI entirely above truth", lo > truth, 0.025)
    ))
  }
  if (!is.null(rej_alt)) out <- c(out, list(prop("Power", rej_alt, NA_real_)))
  if (!is.null(rej_null)) {
    out <- c(out, list(prop("Type I error", rej_null, 0.05),
                       prop("Rejections favouring arm 1", low_null, 0.025)))
  }
  cbind(data.table::data.table(analysis = analysis, method = method), data.table::rbindlist(out))
}

#' Performance of the published CCW analysis in long format (main replicates)
#' @param sim_results Main replicates; @param truth_table From `make_truth_table()`.
#' @param level Nominal level.
#' @return Rows as `perf_rows()`.
ccw_performance_long <- function(sim_results, truth_table, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  tt <- truth_table[scenario == "observed"]
  o <- sim_results[theta == tt$theta & !is.na(estimate)]
  n0 <- sim_results[theta == 0 & !is.na(estimate)]
  perf_rows("CCW", "Wald (robust SE)", o$estimate, o$estimate - z * o$se, o$estimate + z * o$se,
            tt$truth, rej_alt = abs(o$estimate / o$se) > z, rej_null = abs(n0$estimate / n0$se) > z,
            low_null = n0$estimate / n0$se < -z, bias_extra_mcse = tt$truth_mcse)
}

#' Check the simple adjustment in an independent set of simulated studies (paper 2.7.2 step 5c,
#' 2.7.3 step 9)
#' @param adj From `adjust_simple()`; @param val_results Validation replicates (both scenarios).
#' @param truth_table From `make_truth_table()`; @param level Nominal level.
#' @return Rows as `perf_rows()`, analysis "Simple adjustment", one block per CI/test variant.
validate_adjustment <- function(adj, val_results, truth_table, level = 0.95) {
  d <- merge(val_results[!is.na(estimate) & !is.na(se)],
             truth_table[, .(dataset, theta, scenario, truth)], by = c("dataset", "theta"))
  a <- cbind(d, apply_adjustment(d$estimate, d$se, adj, level))
  o <- a[scenario == "observed"]
  n0 <- a[scenario == "null"]
  tr <- o$truth[1]
  # The correction was estimated from the main replicates: its MC error adds to the bias MCSE
  extra <- adj$mean_est_mcse
  rbind(
    perf_rows("Simple adjustment", "Non-parametric", o$est_adj, o$ci_lo_np, o$ci_hi_np, tr,
              o$rej_np, n0$rej_np, n0$low_np, extra),
    perf_rows("Simple adjustment", "Non-parametric (reflected)", o$est_adj, o$ci_lo_np_refl,
              o$ci_hi_np_refl, tr, bias_extra_mcse = extra),
    perf_rows("Simple adjustment", "Parametric", o$est_adj, o$ci_lo_par, o$ci_hi_par, tr,
              o$rej_par, n0$rej_par, n0$low_par, extra),
    perf_rows("CCW", "Wald (robust SE), validation set", o$estimate, o$ci_lo_wald, o$ci_hi_wald,
              tr, o$rej_wald, n0$rej_wald, n0$low_wald)
  )
}

#' Performance of all analyses (CCW, simple and nested adjustment) in one long table
#' @param sim_results Main replicates; @param truth_table From `make_truth_table()`.
#' @param validation From `validate_adjustment()`; @param nested From `summarise_nested()`.
#' @param level Nominal level.
#' @return Rows as `perf_rows()`.
compare_performance <- function(sim_results, truth_table, validation, nested, level = 0.95) {
  rbind(ccw_performance_long(sim_results, truth_table, level), validation, nested$long)
}

#' Quantiles of the studentised statistic T = (estimate - truth) / se (parametric bootstrap-t)
#' @param estimate,se Replicate estimates and SEs; @param truth True value in that scenario.
#' @param level Nominal level.
#' @return Named numeric: q_lo, q_hi.
boot_t_quantiles <- function(estimate, se, truth, level = 0.95) {
  ok <- !is.na(estimate) & !is.na(se) & se > 0
  t <- (estimate[ok] - truth) / se[ok]
  a <- (1 - level) / 2
  q <- stats::quantile(t, c(a, 1 - a), names = FALSE, type = 8)
  c(q_lo = q[1], q_hi = q[2])
}

#' Bootstrap-t CI: invert T using its simulated quantiles
#' @param estimate,se Estimate(s) and SE(s) to interval; @param q From `boot_t_quantiles()`.
#' @return list(lo, hi).
boot_t_ci <- function(estimate, se, q) {
  # T already contains the bias, so the interval is built around the raw estimate
  list(lo = estimate - q[["q_hi"]] * se, hi = estimate - q[["q_lo"]] * se)
}

#' Bootstrap-t and confidence-belt CIs for the observed analysis (optional CI comparison)
#' @param obs_result Observed pipeline result; @param sim_results Main replicates.
#' @param truth_table From `make_truth_table()`; @param belt From `fit_belt()`.
#' @param level Nominal level.
#' @return One-row data.table: ci_lo/hi for bootstrap-t and belt, T quantiles.
extra_cis <- function(obs_result, sim_results, truth_table, belt, level = 0.95) {
  tt <- truth_table[scenario == "observed"]
  o <- sim_results[theta == tt$theta]
  q <- boot_t_quantiles(o$estimate, o$se, tt$truth, level)
  ci_t <- boot_t_ci(obs_result$estimate, obs_result$se, q)
  ci_b <- belt_invert(obs_result$estimate, obs_result$se, belt$curve)
  data.table::data.table(t_q_lo = q[["q_lo"]], t_q_hi = q[["q_hi"]],
                         ci_lo_boot_t = ci_t$lo, ci_hi_boot_t = ci_t$hi,
                         ci_lo_belt = ci_b$lo, ci_hi_belt = ci_b$hi)
}
