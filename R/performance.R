# Performance measures with Monte Carlo SE (paper Step 6; ADEMP P).

#' Label the two scenarios and attach the truth used for each
#' @param truth Truth rows (theta*, 0); @param dataset Dataset key.
#' @return data.table: dataset, scenario, theta, truth, truth_mcse, psi_large_sample.
make_truth_table <- function(truth, dataset) {
  tt <- data.table::copy(truth)
  tt[, `:=`(
    dataset = dataset,
    scenario = data.table::fifelse(theta == 0, "null", "observed"),
    psi_large_sample = truth
  )]
  # With no effect every estimand is exactly 0; psi(0) is kept as a diagnostic only
  tt[scenario == "null", `:=`(truth = 0, truth_mcse = 0)]
  tt[]
}

#' Summarise replicate results against the truth
#' @param sim_results Stacked replicate table from `run_batch()`.
#' @param truth_table From `make_truth_table()`; @param level Nominal confidence level.
#' @param favours_low TRUE if a smaller estimate favours arm 1 (log HR); FALSE for a survival
#'   difference.
#' @return data.table, one row per scenario. Rejection rates (`reject_0`) are reported for the
#'   no-effect scenario only (type I error); no power is computed.
summarise_performance <- function(sim_results, truth_table, level = 0.95, favours_low = TRUE) {
  checkmate::assert_names(names(sim_results), must.include = c("theta", "estimate", "se"))
  d <- merge(sim_results, truth_table, by = c("dataset", "theta"))
  z <- stats::qnorm(1 - (1 - level) / 2)
  # Failures are counted; estimates from failed fits are NA and excluded (ADEMP P4)
  d[, .(
    dataset = dataset[1], truth = truth[1], truth_mcse = truth_mcse[1],
    psi_large_sample = psi_large_sample[1],
    B = .N, n_fail = sum(status == "error"), n_warn = sum(status == "warning"),
    mean_est = mean(estimate, na.rm = TRUE),
    bias = mean(estimate - truth, na.rm = TRUE),
    bias_mcse = stats::sd(estimate, na.rm = TRUE) / sqrt(sum(!is.na(estimate))),
    emp_se = stats::sd(estimate, na.rm = TRUE),
    mod_se = sqrt(mean(se^2, na.rm = TRUE)),
    mse = mean((estimate - truth)^2, na.rm = TRUE),
    mse_mcse = stats::sd((estimate - truth)^2, na.rm = TRUE) / sqrt(sum(!is.na(estimate))),
    coverage = mean(abs(estimate - truth) <= z * se, na.rm = TRUE),
    reject_0 = if (scenario[1] == "null") mean(abs(estimate) > z * se, na.rm = TRUE) else NA_real_,
    reject_0_fav = if (scenario[1] != "null") NA_real_ else if (favours_low) {
      mean(estimate / se < -z, na.rm = TRUE)
    } else {
      mean(estimate / se > z, na.rm = TRUE)
    },
    mean_outcomes = mean(n_outcomes)
  ), by = .(scenario, theta)][, `:=`(
    # Truth is itself simulated (observed scenario): its MC error adds to the bias MC error
    bias_mcse_total = sqrt(bias_mcse^2 + truth_mcse^2),
    emp_se_mcse = emp_se / sqrt(2 * (B - n_fail - 1)),
    rel_err_modse = mod_se / emp_se - 1,
    rmse = sqrt(mse),
    rmse_mcse = mse_mcse / (2 * sqrt(mse)),
    std_bias = bias / emp_se,
    coverage_mcse = sqrt(coverage * (1 - coverage) / (B - n_fail)),
    reject_0_mcse = sqrt(reject_0 * (1 - reject_0) / (B - n_fail))
  )][order(-scenario)]
}

#' Cross-check of headline measures with rsimsum
#' @param sim_results Replicate table; @param truth_table From `make_truth_table()`.
#' @return Tidy rsimsum summary (bias, empse, modelse, cover, mse) by scenario.
rsimsum_check <- function(sim_results, truth_table) {
  d <- as.data.frame(merge(sim_results, truth_table, by = c("dataset", "theta")))
  d <- d[!is.na(d$estimate), ]
  out <- lapply(split(d, d$scenario), function(x) {
    s <- rsimsum::simsum(x, estvarname = "estimate", true = x$truth[1], se = "se")
    tb <- rsimsum::tidy(summary(s), stats = c("bias", "empse", "modelse", "cover", "mse"))
    cbind(dataset = x$dataset[1], scenario = x$scenario[1], tb)
  })
  data.table::rbindlist(out)
}
