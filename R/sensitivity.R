# Sensitivity of the QBA adjustment to the prediction model (paper Section 3.8): the same
# pipeline, seeds and simple (non-nested) adjustment under each outcome model, compared side by
# side, plus the conservative combination across models.

#' Fit statistics of the prediction models
#' @param dgms Named list of fitted DGMs (names end in the dataset key).
#' @param labels Named character vector: dataset key -> display label.
#' @return data.table, one row per model: dataset, label, model, n_par, loglik, aic, theta_hat,
#'   theta_se, shape (text), lr_stat and lr_p (1-df likelihood-ratio test against the Weibull).
sens_model_fit <- function(dgms, labels) {
  f2 <- function(x) formatC(x, format = "f", digits = 2)
  data.table::rbindlist(lapply(names(labels), function(ds) {
    m <- dgms[[grep(paste0("(^|_)", ds, "$"), names(dgms))[1]]]
    model <- dgm_model_label(m)
    # Exponential: restricted Weibull; generalized gamma: the Weibull is its Q = 1 member
    lr <- switch(model,
      Exponential = 2 * (m$loglik_weibull - m$loglik),
      `Generalized gamma` = 2 * (m$loglik - m$loglik_weibull),
      NA_real_
    )
    data.table::data.table(
      dataset = ds, label = labels[[ds]], model = model,
      n_par = length(m$se), loglik = m$loglik, aic = 2 * length(m$se) - 2 * m$loglik,
      theta_hat = m$theta_hat, theta_se = m$se[["theta"]],
      shape = switch(model,
        Exponential = "shape fixed at 1",
        `Generalized gamma` = paste0("sigma = ", f2(m$sigma), ", Q = ", f2(m$Q)),
        paste0("shape = ", f2(m$rho))
      ),
      lr_stat = lr,
      lr_p = if (is.na(lr)) NA_real_ else stats::pchisq(max(lr, 0), 1, lower.tail = FALSE)
    )
  }))
}

#' Monte Carlo SE of an empirical quantile from the binomial ranks around it
#' @param x Simulated values; @param p Probability.
#' @return Half the distance between the order statistics one binomial SD either side of n * p.
quantile_mcse <- function(x, p) {
  x <- sort(x[!is.na(x)])
  n <- length(x)
  half <- sqrt(n * p * (1 - p))
  lo <- max(1L, floor(n * p - half))
  hi <- min(n, ceiling(n * p + half))
  (x[hi] - x[lo]) / 2
}

#' Performance and adjustment under each prediction model, side by side
#' @param adjustment_all Stacked `adjust_simple()` rows; @param performance_all Stacked
#'   `summarise_performance()` rows; @param perf_long Stacked `compare_performance()` rows with a
#'   column `target` whose name ends in the dataset key; @param sim_all Stacked main replicates.
#' @param primary Dataset key of the primary model (differences are taken against it).
#' @param level Nominal level.
#' @return data.table, one row per model. `d_bias` is the bias minus that of the primary model;
#'   its MCSE pairs the replicates (common random numbers) and adds both truth MCSEs.
sens_compare <- function(adjustment_all, performance_all, perf_long, sim_all, primary,
                         level = 0.95) {
  checkmate::assert_choice(primary, adjustment_all$dataset)
  a <- (1 - level) / 2
  po <- performance_all[scenario == "observed"]
  pn <- performance_all[scenario == "null"]
  est_of <- function(ds, scen) {
    th <- performance_all[dataset == ds & scenario == scen, theta]
    sim_all[dataset == ds & theta == th & !is.na(estimate), .(rep, estimate)]
  }
  ref <- est_of(primary, "observed")
  adj_perf <- function(ds, what) {
    r <- perf_long[grepl(paste0("(^|_)", ds, "$"), target) & analysis == "Simple adjustment" &
                     method == "Non-parametric" & measure == what]
    if (nrow(r)) r[1, .(value, mcse)] else data.table::data.table(value = NA_real_, mcse = NA_real_)
  }
  out <- data.table::rbindlist(lapply(adjustment_all$dataset, function(ds) {
    ad <- adjustment_all[dataset == ds]
    o <- po[dataset == ds]
    n0 <- pn[dataset == ds]
    e <- est_of(ds, "observed")
    e0 <- est_of(ds, "null")
    pair <- merge(e, ref, by = "rep", suffixes = c("", "_ref"))
    d_mean_mcse <- stats::sd(pair$estimate - pair$estimate_ref) / sqrt(nrow(pair))
    cov <- adj_perf(ds, "95% CI coverage")
    t1 <- adj_perf(ds, "Type I error")
    data.table::data.table(
      dataset = ds, theta_cal = ad$theta_cal, truth = ad$truth, truth_mcse = ad$truth_mcse,
      n_sim = ad$n_sim, n_fail = o$n_fail + n0$n_fail,
      bias = ad$bias, bias_mcse = ad$bias_mcse,
      d_mean_mcse = if (ds == primary) 0 else d_mean_mcse,
      emp_se = ad$emp_se, emp_se_mcse = o$emp_se_mcse, mod_se = ad$mod_se,
      se_inflation = ad$se_inflation,
      coverage_wald = o$coverage, coverage_wald_mcse = o$coverage_mcse,
      type1_wald = n0$reject_0, type1_wald_mcse = n0$reject_0_mcse,
      q_lo = ad$q_lo, q_lo_mcse = quantile_mcse(e$estimate, a),
      q_hi = ad$q_hi, q_hi_mcse = quantile_mcse(e$estimate, 1 - a),
      null_q_lo = ad$null_q_lo, null_q_lo_mcse = quantile_mcse(e0$estimate, a),
      null_q_hi = ad$null_q_hi, null_q_hi_mcse = quantile_mcse(e0$estimate, 1 - a),
      z_crit_par = ad$z_crit_par,
      est_obs = ad$est_obs, est_adj = ad$est_adj,
      ci_lo_np = ad$ci_lo_np, ci_hi_np = ad$ci_hi_np, p_np = ad$p_np, reject_np = ad$reject_np,
      ci_lo_par = ad$ci_lo_par, ci_hi_par = ad$ci_hi_par, p_par = ad$p_par,
      coverage_adj = cov$value, coverage_adj_mcse = cov$mcse,
      type1_adj = t1$value, type1_adj_mcse = t1$mcse
    )
  }))
  pr <- out[dataset == primary]
  out[, d_bias := bias - pr$bias]
  out[, d_bias_mcse := data.table::fifelse(
    dataset == primary, 0, sqrt(d_mean_mcse^2 + truth_mcse^2 + pr$truth_mcse^2)
  )]
  out[, d_mean_mcse := NULL]
  out[]
}

#' Conservative adjustment across prediction models (paper Section 3.8)
#' Adjusted estimate closest to no effect, the widest CI limits and rejection limits on each
#' side, and the largest p-value; H0 is rejected only if it is rejected under every model.
#' @param adjustment_all Stacked `adjust_simple()` rows; @param keep Dataset keys to combine
#'   (models with a clearly worse fit should be left out; default all).
#' @return One-row data.table: models, est_adj, ci_lo_np, ci_hi_np, null_q_lo, null_q_hi,
#'   z_crit_par, p_np, reject_np.
sens_conservative <- function(adjustment_all, keep = adjustment_all$dataset) {
  checkmate::assert_subset(keep, adjustment_all$dataset)
  ad <- adjustment_all[dataset %in% keep]
  data.table::data.table(
    models = paste(ad$dataset, collapse = ", "),
    est_adj = ad$est_adj[which.min(abs(ad$est_adj))],
    ci_lo_np = min(ad$ci_lo_np), ci_hi_np = max(ad$ci_hi_np),
    null_q_lo = min(ad$null_q_lo), null_q_hi = max(ad$null_q_hi),
    z_crit_par = max(ad$z_crit_par),
    p_np = max(ad$p_np), reject_np = all(ad$reject_np)
  )
}

#' Figure: Kaplan-Meier curve of the outcome against the marginal survival of each model
#' @param data Base data; @param dgms Named list of fitted DGMs; @param labels Dataset labels.
#' @return ggplot.
plot_sens_fit <- function(data, dgms, labels) {
  times <- seq(0, max(data$time), length.out = 200)
  km <- survival::survfit(survival::Surv(data$time, data$status) ~ 1)
  d <- data.table::rbindlist(lapply(names(labels), function(ds) {
    m <- dgms[[grep(paste0("(^|_)", ds, "$"), names(dgms))[1]]]
    data.table::data.table(time = times, surv = dgm_marginal_surv(m, data, times),
                           model = dgm_model_label(m))
  }))
  ggplot2::ggplot(d, ggplot2::aes(time, surv)) +
    ggplot2::geom_step(data = data.table::data.table(time = c(0, km$time), surv = c(1, km$surv)),
                       colour = viz_cols[["ink"]], linewidth = 0.9) +
    ggplot2::geom_line(ggplot2::aes(colour = model, linetype = model), linewidth = 0.8) +
    ggplot2::labs(x = "Days since time zero", y = "Proportion alive (black: Kaplan-Meier)",
                  colour = NULL, linetype = NULL) +
    theme_qba()
}

#' Figure: centred sampling distributions of the estimate under each model, both scenarios
#' @param sim_all Stacked main replicates; @param performance_all Stacked performance rows.
#' @param labels Dataset labels.
#' @return ggplot.
plot_sens_distribution <- function(sim_all, performance_all, labels) {
  d <- merge(sim_all[!is.na(estimate) & dataset %in% names(labels)],
             performance_all[, .(dataset, theta, scenario, mean_est)], by = c("dataset", "theta"))
  d[, `:=`(
    centred = estimate - mean_est,
    model = factor(labels[dataset], levels = labels),
    panel = data.table::fifelse(scenario == "null", "No effect", "Calibrated effect")
  )]
  ggplot2::ggplot(d, ggplot2::aes(centred, colour = model, linetype = model)) +
    ggplot2::geom_density(linewidth = 0.7, key_glyph = "path") +
    ggplot2::facet_wrap(~panel) +
    ggplot2::labs(x = "Estimate minus its mean (log HR)", y = "Density", colour = NULL,
                  linetype = NULL) +
    theme_qba()
}
