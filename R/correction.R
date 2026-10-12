# Bias-corrected estimate, widened CI and simulation-calibrated p-value (paper Step 7).

#' Correct the observed analysis using the observed-effect and null scenarios
#' @param obs_result Observed pipeline result; @param performance From `summarise_performance()`.
#' @param sim_results Replicate table; @param dataset Dataset key; @param level Nominal level.
#' @return One-row data.table on the log-HR scale.
correct_estimate <- function(obs_result, performance, sim_results, dataset, level = 0.95) {
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
  data.table::data.table(
    dataset = dataset,
    est_obs = obs_result$estimate, ci_lo_obs = obs_result$ci_lo, ci_hi_obs = obs_result$ci_hi,
    p_obs = obs_result$p_value,
    est_corr = est_c,
    se_inflation = infl,
    ci_lo_corr = est_c - z * infl * obs_result$se,
    ci_hi_corr = est_c + z * infl * obs_result$se,
    # Critical value that gives a 5% rejection rate when treatment has no effect
    z_crit_calibrated = stats::quantile(z_null, level, names = FALSE),
    p_calibrated = mean(z_null >= z_obs)
  )
}
