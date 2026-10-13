# Bias-corrected estimate, widened CI and simulation-calibrated p-value (paper Step 7).

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

#' Correct the observed analysis using the observed-effect and null scenarios
#' @param obs_result Observed pipeline result; @param performance From `summarise_performance()`.
#' @param sim_results Replicate table; @param dataset Dataset key; @param level Nominal level.
#' @param belt From `fit_belt()` (optional): adds the confidence-belt CI.
#' @return One-row data.table on the log-HR scale.
correct_estimate <- function(obs_result, performance, sim_results, dataset, level = 0.95,
                             belt = NULL) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  po <- performance[scenario == "observed"]
  # Paper 2.7: truth = observed estimate, so subtracting the bias is the grid-search answer
  # to first order (bias is locally constant in the truth)
  est_c <- obs_result$estimate - po$bias
  # Widen the CI by empirical/model SE so that coverage returns to nominal
  infl <- max(1, po$emp_se / po$mod_se)
  null <- sim_results[theta == 0 & !is.na(estimate)]
  z_null <- abs(null$estimate / null$se)
  z_obs <- abs(obs_result$estimate / obs_result$se)
  # Bootstrap-t from the observed-effect scenario (truth = observed estimate): asymmetric CI
  # that also captures the estimate-SE correlation the symmetric widening misses
  obs_sim <- sim_results[theta == po$theta]
  q <- boot_t_quantiles(obs_sim$estimate, obs_sim$se, po$truth, level)
  ci_t <- boot_t_ci(obs_result$estimate, obs_result$se, q)
  ci_b <- if (is.null(belt)) list(lo = NA_real_, hi = NA_real_) else
    belt_invert(obs_result$estimate, obs_result$se, belt$curve)
  data.table::data.table(
    dataset = dataset,
    est_obs = obs_result$estimate, ci_lo_obs = obs_result$ci_lo, ci_hi_obs = obs_result$ci_hi,
    p_obs = obs_result$p_value,
    est_corr = est_c,
    se_inflation = infl,
    ci_lo_corr = est_c - z * infl * obs_result$se,
    ci_hi_corr = est_c + z * infl * obs_result$se,
    t_q_lo = q[["q_lo"]], t_q_hi = q[["q_hi"]],
    ci_lo_boot_t = ci_t$lo, ci_hi_boot_t = ci_t$hi,
    ci_lo_belt = ci_b$lo, ci_hi_belt = ci_b$hi,
    # Critical value that gives a 5% rejection rate when treatment has no effect
    z_crit_calibrated = stats::quantile(z_null, level, names = FALSE),
    p_calibrated = mean(z_null >= z_obs)
  )
}
