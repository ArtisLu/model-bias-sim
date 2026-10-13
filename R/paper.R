# Publication outputs for the SynPUF case study (paper Section 4): baseline table in the
# style of the Maringe Table S2, results tables, performance after QBA adjustment, figures
# in the style of the earlier simulation papers (grey strips, framed panels, MCSE bars).

#' Standardised mean difference (%) for a binary indicator or a continuous variable
#' @param x Values; @param g Logical group (TRUE = first column).
#' @return Scalar SMD in percent.
smd_pct <- function(x, g) {
  m1 <- mean(x[g]); m0 <- mean(x[!g])
  v <- if (all(x %in% 0:1)) {
    (m1 * (1 - m1) + m0 * (1 - m0)) / 2
  } else {
    (stats::var(x[g]) + stats::var(x[!g])) / 2
  }
  100 * (m1 - m0) / sqrt(v)
}

#' Read the raw base data file (read-only source)
#' @param file Path; @param ds Dataset key.
#' @return data.table as stored in the file.
read_raw_data <- function(file, ds) {
  switch(ds,
    maringe = data.table::fread(file),
    synpuf = data.table::as.data.table(haven::read_sas(file))
  )
}

#' Table rows in the Table S2 layout (helpers shared by both datasets)
#' @param g Logical group indicator (TRUE = first column).
#' @return List of row-building functions.
s2_rows <- function(g) {
  f1 <- function(x) formatC(x, format = "f", digits = 1)
  row <- function(label, indent, n1, p1, n0, p0, smd) {
    data.table::data.table(label = label, indent = indent, n1 = n1, p1 = p1, n0 = n0, p0 = p0, smd = smd)
  }
  cat_row <- function(label, x, smd = TRUE) {
    row(label, 1L, format(sum(x[g])), f1(100 * mean(x[g])), format(sum(x[!g])), f1(100 * mean(x[!g])),
        if (smd) f1(smd_pct(as.numeric(x), g)) else "")
  }
  list(
    total = function() row("Total", 0L, format(sum(g)), "100.0", format(sum(!g)), "100.0", ""),
    head = function(label) row(label, 0L, "", "", "", "", ""),
    cat = cat_row,
    # Binary characteristics get one SMD, on the first level (as in Table S2)
    bin = function(label_a, label_b, x) rbind(cat_row(label_a, x), cat_row(label_b, !x, smd = FALSE)),
    qhead = function() row("", 0L, "Median", "Q1 Q3", "Median", "Q1 Q3", ""),
    q = function(label, x, digits = 1, smd = TRUE) {
      q <- function(v) formatC(stats::quantile(v, c(0.5, 0.25, 0.75), names = FALSE, na.rm = TRUE),
                               format = "f", digits = digits)
      a <- q(x[g])
      b <- if (any(!g & !is.na(x))) q(x[!g]) else c("", "", "")
      row(label, 0L, a[1], paste(a[2], a[3]), b[1], paste(b[2], b[3]),
          if (smd) f1(smd_pct(x, g)) else "")
    }
  )
}

#' Baseline characteristics by observed treatment (Table S2 layout)
#' @param raw Raw data from `read_raw_data()`; @param ds Dataset key; @param spec Dataset config.
#' @return data.table: label, indent, n1, p1, n0, p0, smd (text), attribute `groups`.
make_baseline_table <- function(raw, ds, spec) {
  d <- data.table::copy(raw)
  if (ds == "maringe") {
    g <- !is.na(d$timetosurgery)
    r <- s2_rows(g)
    tb <- rbind(
      r$total(),
      r$head("Sex"), r$bin("Men", "Women", d$sex == 1),
      r$head("Performance status"),
      r$cat("0", d$perf == 0), r$cat("1", d$perf == 1), r$cat("2", d$perf == 2),
      r$head("Charlson's comorbidity index"),
      r$cat("0", d$charlson == 0), r$cat("1", d$charlson == 1), r$cat("2", d$charlson == 2),
      r$head("TNM stage"), r$bin("I", "II", d$stage == 1),
      r$head("Route to diagnosis"), r$bin("Emergency", "Non-emergency", d$emergency == 1),
      r$qhead(),
      r$q("Age at diagnosis (years)", d$age),
      r$q("Deprivation score", d$deprivation, digits = 2),
      r$q("Time to surgery (days)", data.table::fifelse(g, d$timetosurgery, NA_real_), 0, FALSE)
    )
    attr(tb, "groups") <- c("Received surgery", "Did not receive surgery")
  } else {
    h <- spec$horizon
    # Initiation counts only within the follow-up and before the outcome or loss
    g <- !is.na(d$RX_start) & d$RX_start <= pmin(d$followup, h)
    r <- s2_rows(g)
    tb <- rbind(
      r$total(),
      r$head("Sex"), r$bin("Men", "Women", d$sex == "1"),
      r$head("Race"),
      r$cat("White", d$race == "1"), r$cat("Black", d$race == "2"),
      r$cat("Other", d$race == "3"), r$cat("Hispanic", d$race == "5"),
      r$head("End-stage renal disease"), r$bin("Yes", "No", d$renal == "Y"),
      r$head("Outcome within 180 days"),
      r$cat("Readmission or death", d$outcome == 1 & d$followup <= h, smd = FALSE),
      r$cat("Lost to follow-up", d$outcome == 0 & d$followup < h, smd = FALSE),
      r$qhead(),
      r$q("Age at MI (years)", d$baseline_age, 0),
      r$q("Time to treatment start (days)", data.table::fifelse(g, d$RX_start, NA_real_), 0, FALSE)
    )
    attr(tb, "groups") <- c("Started treatment", "Did not start treatment")
  }
  tb
}

#' Weibull outcome and censoring model parameters (supplementary table)
#' @param dgm Fitted DGM with `cens` element; @param spec Dataset config (outcome label).
#' @return data.table: model, term, hr, ci (text), n_events.
make_dgm_table <- function(dgm, spec) {
  # Display labels and the unit each HR refers to (coefficients are per unit of the model term)
  lab <- c(
    age_c = "Age (per 5 years)", sex = "Female", perf1 = "Performance status 1 vs 0",
    perf2 = "Performance status 2 vs 0", stage = "TNM stage II vs I",
    deprivation = "Deprivation score (per 0.1)", charlson1 = "Charlson index 1 vs 0",
    charlson2 = "Charlson index 2 vs 0", emergency = "Emergency presentation",
    renal = "End-stage renal disease", nonwhite = "Non-white race",
    theta = "On treatment (time-varying)"
  )
  # Age is centred and scaled in data.R; the unit is recorded in the config
  lab["age_c"] <- paste0("Age (per ", spec$dgm$age_unit, " years)")
  scale <- c(deprivation = 0.1)
  one <- function(m, model) {
    terms <- c(names(m$gamma), "theta")
    k <- ifelse(terms %in% names(scale), scale[terms], 1)
    b <- k * c(m$gamma, theta = m$theta_hat)
    se <- k * m$se[terms]
    lb <- ifelse(terms %in% names(lab), lab[terms], terms)
    f2 <- function(x) formatC(x, format = "f", digits = 2)
    f4 <- function(x) formatC(x, format = "f", digits = 4)
    rbind(
      data.table::data.table(model = model, term = "Weibull shape (rho)", hr = f2(m$rho), ci = ""),
      data.table::data.table(
        model = model, term = unname(lb), hr = f4(exp(b)),
        ci = paste0(f4(exp(b - 1.96 * se)), " to ", f4(exp(b + 1.96 * se)))
      )
    )[, n_events := m$n_events]
  }
  out <- one(dgm, paste0("Outcome: ", tolower(spec$outcome_label)))
  if (identical(dgm$cens$model, "weibull")) out <- rbind(out, one(dgm$cens, "Censoring: loss to follow-up"))
  out
}

#' Performance of the published CCW analysis and of the QBA-adjusted analysis
#' @param sim_results Replicates (one dataset, both scenarios).
#' @param truth_table From `make_truth_table()`; @param level Nominal level.
#' @return data.table: analysis, scenario, measure, value, mcse, ref.
evaluate_adjustment <- function(sim_results, truth_table, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  d <- merge(sim_results[!is.na(estimate)], truth_table[, .(dataset, theta, scenario, truth)],
             by = c("dataset", "theta"))
  d[, fold := rep %% 2L]
  # Two-fold cross-fitting: corrections learned on the other half of the replicates, so the
  # adjusted procedure is not judged on the replicates that defined it
  corr <- data.table::rbindlist(lapply(0:1, function(k) {
    o <- d[fold != k & scenario == "observed"]
    n0 <- d[fold != k & scenario == "null"]
    data.table::data.table(
      fold = k, bias = mean(o$estimate - o$truth),
      infl = max(1, stats::sd(o$estimate) / sqrt(mean(o$se^2))),
      z_crit = stats::quantile(abs(n0$estimate / n0$se), level, names = FALSE)
    )
  }))
  d <- merge(d, corr, by = "fold")
  # The adjusted procedure is the same whatever the true effect, so it is applied in both scenarios
  d[, `:=`(est_adj = estimate - bias, se_adj = infl * se)]
  p_row <- function(analysis, scenario, measure, x, ref) {
    n <- length(x)
    m <- mean(x)
    mc <- if (all(x %in% 0:1)) sqrt(m * (1 - m) / n) else stats::sd(x) / sqrt(n)
    data.table::data.table(analysis = analysis, scenario = scenario, measure = measure,
                           value = m, mcse = mc, ref = ref, n = n)
  }
  o <- d[scenario == "observed"]
  n0 <- d[scenario == "null"]
  rbind(
    p_row("CCW", "observed", "Bias (log HR)", o$estimate - o$truth, 0),
    p_row("CCW", "observed", "95% CI coverage", as.integer(abs(o$estimate - o$truth) <= z * o$se), level),
    p_row("CCW", "null", "Type I error", as.integer(abs(n0$estimate / n0$se) > z), 1 - level),
    p_row("CCW", "observed", "Power", as.integer(abs(o$estimate / o$se) > z), NA_real_),
    p_row("CCW + QBA adjustment", "observed", "Bias (log HR)", o$est_adj - o$truth, 0),
    p_row("CCW + QBA adjustment", "observed", "95% CI coverage",
          as.integer(abs(o$est_adj - o$truth) <= z * o$se_adj), level),
    p_row("CCW + QBA adjustment", "null", "Type I error",
          as.integer(abs(n0$estimate / n0$se) > n0$z_crit), 1 - level),
    p_row("CCW + QBA adjustment", "observed", "Power",
          as.integer(abs(o$estimate / o$se) > o$z_crit), NA_real_)
  )
}

#' Observed analysis, reference model, calibration and performance in one list of tables
#' @param obs Observed pipeline result; @param dgm Fitted DGM; @param cal Calibration.
#' @param truth_table,performance,correction,perf_adj Pipeline targets for one dataset.
#' @param spec Dataset config (`arm1` names the strategy coded 1).
#' @return list(analysis = data.table, performance = data.table, dgm = data.table).
make_results_table <- function(obs, dgm, cal, truth_table, performance, correction, perf_adj, spec) {
  f2 <- function(x) formatC(x, format = "f", digits = 2)
  f3 <- function(x) formatC(x, format = "f", digits = 3)
  f4 <- function(x) formatC(x, format = "f", digits = 4)
  pc <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
  # HRs and their CIs with 4 decimals, type I error with 3 (percent)
  hr_ci <- function(e, lo, hi) paste0(f4(exp(e)), " (", f4(exp(lo)), " to ", f4(exp(hi)), ")")
  po <- performance[scenario == "observed"]
  pn <- performance[scenario == "null"]
  tn <- truth_table[scenario == "null"]
  cr <- correction
  analysis <- data.table::data.table(
    item = c(
      "Observed CCW HR (95% CI)", "Observed p-value",
      "Weibull outcome model: on-treatment HR (exp(theta-hat))",
      "Calibrated generating HR, exp(theta*)",
      "True CCW HR in the simulated world, exp(psi(theta*)) (MCSE of psi)",
      "Large-sample CCW HR with no treatment effect, exp(psi(0))",
      "QBA-adjusted HR (95% CI)", "CI inflation factor (empirical / model SE)",
      "Calibrated critical value |z| for a 5% test", "Simulation-calibrated p-value"
    ),
    value = c(
      hr_ci(obs$estimate, obs$ci_lo, obs$ci_hi),
      formatC(obs$p_value, format = "f", digits = 3),
      f4(exp(dgm$theta_hat)), f4(exp(cal$theta_star)),
      paste0(f4(exp(po$truth)), " (", f3(po$truth_mcse), ")"),
      f4(exp(tn$psi_large_sample)),
      hr_ci(cr$est_corr, cr$ci_lo_corr, cr$ci_hi_corr),
      f2(cr$se_inflation), f2(cr$z_crit_calibrated),
      formatC(cr$p_calibrated, format = "f", digits = 3)
    )
  )
  # Unadjusted CCW performance with full MCSE (incl. truth MC error) plus cross-fitted adjusted
  pa <- perf_adj[analysis == "CCW + QBA adjustment"]
  getp <- function(m) pa[measure == m]
  fmt <- function(v, mc, pct, d = 1) {
    if (pct) paste0(pc(v, d), " (", pc(mc, d), ")") else paste0(f3(v), " (", f3(mc), ")")
  }
  performance_tab <- data.table::data.table(
    measure = c(
      "Mean estimate, log HR (true value)", "Bias, log HR", "Relative bias, % of empirical SE",
      "Empirical SE", "Model-based (robust) SE", "Relative error in model SE, %", "RMSE",
      "95% CI coverage, %", "Power, % (true HR = observed)", "Type I error, % (no effect)",
      paste0("Rejections favouring ", spec$arm1, ", % (no effect)"), "Failed replicates, n"
    ),
    ccw = c(
      paste0(f3(po$mean_est), " (", f3(po$truth), ")"),
      fmt(po$bias, po$bias_mcse_total, FALSE),
      pc(po$std_bias),
      fmt(po$emp_se, po$emp_se_mcse, FALSE), f3(po$mod_se), pc(po$rel_err_modse), f3(po$rmse),
      fmt(po$coverage, po$coverage_mcse, TRUE),
      fmt(po$reject_0, po$reject_0_mcse, TRUE),
      fmt(pn$reject_0, pn$reject_0_mcse, TRUE, 3),
      pc(pn$reject_0_lower, 3),
      format(po$n_fail + pn$n_fail)
    ),
    adjusted = c(
      # With two equal folds the cross-fitted bias correction averages to exactly zero
      "", "0 (by construction)", "", "", "", "", "",
      fmt(getp("95% CI coverage")$value, getp("95% CI coverage")$mcse, TRUE),
      fmt(getp("Power")$value, getp("Power")$mcse, TRUE),
      fmt(getp("Type I error")$value, getp("Type I error")$mcse, TRUE, 3), "", ""
    )
  )
  list(analysis = analysis, performance = performance_tab, dgm = make_dgm_table(dgm, spec))
}

#' Theme of the earlier simulation papers: framed panels, grey strips, legend on the right
theme_paper <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", linewidth = 0.8, fill = NA),
      panel.grid.minor = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.4),
      strip.background = ggplot2::element_rect(fill = "grey78", colour = "black", linewidth = 0.8),
      strip.text = ggplot2::element_text(face = "bold", size = base_size + 2, colour = "black"),
      axis.text = ggplot2::element_text(colour = "grey30"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.6, "lines")
    )
}

paper_cols <- c("CCW" = "black", "CCW + QBA adjustment" = "red")

#' Figure: performance of CCW before and after QBA adjustment (points +- 1.96 MCSE)
#' @param perf_adj From `evaluate_adjustment()`.
#' @return ggplot.
plot_paper_performance <- function(perf_adj) {
  d <- data.table::copy(perf_adj)
  d[, measure := factor(measure, levels = c("Bias (log HR)", "95% CI coverage", "Type I error", "Power"))]
  d[, analysis := factor(analysis, levels = names(paper_cols))]
  refs <- unique(d[!is.na(ref), .(measure, ref)])
  ggplot2::ggplot(d, ggplot2::aes(analysis, value, colour = analysis)) +
    ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = ref), linetype = "dashed") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = value - 1.96 * mcse, ymax = value + 1.96 * mcse),
                            linewidth = 1.1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~measure, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = paper_cols, name = "Method:") +
    ggplot2::scale_x_discrete(labels = c("CCW", "CCW +\nadjustment")) +
    ggplot2::labs(x = NULL, y = "Estimate with 95% Monte Carlo interval") +
    theme_paper()
}

#' Figure: distribution of simulated estimates (true HR = observed) and of null z statistics
#' @param sim_results Replicates; @param performance Performance table (one dataset).
#' @param z_crit Calibrated critical value (from the correction step).
#' @return ggplot.
plot_paper_distribution <- function(sim_results, performance, z_crit) {
  po <- performance[scenario == "observed"]
  a <- sim_results[theta == po$theta & !is.na(estimate), .(x = estimate, panel = "True HR = observed estimate")]
  b <- sim_results[theta == 0 & !is.na(estimate), .(x = estimate / se, panel = "No treatment effect")]
  d <- rbind(a, b)
  d[, panel := factor(panel, levels = c("True HR = observed estimate", "No treatment effect"))]
  lines <- rbind(
    data.table::data.table(panel = "True HR = observed estimate", x = c(po$truth, po$mean_est),
                           what = c("True value", "Mean of estimates")),
    data.table::data.table(panel = "No treatment effect", x = c(-1.96, 1.96, -z_crit, z_crit),
                           what = rep(c("Nominal critical value", "Calibrated critical value"), each = 2))
  )
  lines[, panel := factor(panel, levels = levels(d$panel))]
  ggplot2::ggplot(d, ggplot2::aes(x)) +
    ggplot2::geom_histogram(bins = 60, fill = "grey70", colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = x, colour = what, linetype = what),
                        linewidth = 0.9) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::scale_colour_manual(
      values = c("True value" = "black", "Mean of estimates" = "red",
                 "Nominal critical value" = "black", "Calibrated critical value" = "red"),
      name = "Reference:"
    ) +
    ggplot2::scale_linetype_manual(
      values = c("True value" = "solid", "Mean of estimates" = "dashed",
                 "Nominal critical value" = "solid", "Calibrated critical value" = "dashed"),
      name = "Reference:"
    ) +
    ggplot2::labs(x = "Estimated log HR (left) / z statistic (right)", y = "Number of simulated studies") +
    theme_paper()
}

#' Marginal survival implied by a Weibull model with time-varying treatment
#' @param m Fitted Weibull list; @param data Base data; @param times Time grid.
#' @return Numeric vector: mean over patients of S_i(t).
weibull_marginal_surv <- function(m, data, times) {
  X <- dgm_design(data, m$covariates)
  eh <- exp(m$alpha + drop(X %*% m$gamma))
  vapply(times, function(t) {
    pre <- pmin(data$s, t)
    H <- eh * (pre^m$rho + (data$s < t) * exp(m$theta_hat) * (t^m$rho - pre^m$rho))
    mean(exp(-H))
  }, numeric(1))
}

#' Figure: fit of the Weibull outcome and censoring models (Kaplan-Meier vs model)
#' @param data Base data; @param dgm Fitted DGM; @param spec Dataset config (outcome label).
#' @return ggplot.
plot_paper_fit <- function(data, dgm, spec) {
  times <- seq(0, max(data$time), length.out = 200)
  km <- function(status, panel) {
    f <- survival::survfit(survival::Surv(data$time, status) ~ 1)
    data.table::data.table(time = c(0, f$time), surv = c(1, f$surv), what = "Kaplan-Meier", panel = panel)
  }
  wb <- function(m, panel) {
    data.table::data.table(time = times, surv = weibull_marginal_surv(m, data, times),
                           what = "Weibull model", panel = panel)
  }
  p1 <- spec$outcome_label
  p2 <- "Loss to follow-up"
  d <- rbind(km(data$status, p1), wb(dgm, p1))
  if (identical(dgm$cens$model, "weibull")) d <- rbind(d, km(data$cens_event, p2), wb(dgm$cens, p2))
  d[, panel := factor(panel, levels = c(p1, p2))]
  ggplot2::ggplot(d, ggplot2::aes(time, surv, colour = what, linetype = what)) +
    ggplot2::geom_step(data = d[what == "Kaplan-Meier"], linewidth = 0.9) +
    ggplot2::geom_line(data = d[what == "Weibull model"], linewidth = 0.9) +
    ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c("Kaplan-Meier" = "black", "Weibull model" = "red"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("Kaplan-Meier" = "solid", "Weibull model" = "dashed"), name = NULL) +
    ggplot2::labs(x = "Days since time zero", y = "Proportion free of the event") +
    theme_paper() +
    ggplot2::theme(legend.position = "bottom")
}

#' Save all paper figures as PNG (600 dpi) and EPS
#' @param performance,perf_adj,sim_results,data,dgm Pipeline targets (one dataset).
#' @param z_crit Calibrated critical value; @param spec Dataset config; @param dir Output folder.
#' @return Character vector of file paths.
save_paper_figures <- function(performance, perf_adj, sim_results, data, dgm, z_crit, spec, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  figs <- list(
    fig_performance = list(plot_paper_performance(perf_adj), 13, 4.2),
    fig_distribution = list(plot_paper_distribution(sim_results, performance, z_crit), 12, 4.6),
    fig_model_fit = list(
      plot_paper_fit(data, dgm, spec), if (identical(dgm$cens$model, "weibull")) 10 else 6, 4.6
    )
  )
  unlist(lapply(names(figs), function(nm) {
    f <- figs[[nm]]
    png <- file.path(dir, paste0(nm, ".png"))
    eps <- file.path(dir, paste0(nm, ".eps"))
    ggplot2::ggsave(png, f[[1]], width = f[[2]], height = f[[3]], dpi = 600, bg = "white")
    ggplot2::ggsave(eps, f[[1]], width = f[[2]], height = f[[3]], device = grDevices::cairo_ps)
    c(png, eps)
  }))
}

#' Covariate balance between the strategy arms at the end of the grace period (Figure S2 of
#' Maringe et al. 2020): standardised differences before and after censoring weights
#' @param data Base data; @param spec Pipeline spec.
#' @param covariates Covariates to show (factors are shown per level).
#' @return data.table: term, weighting, smd (%), n0, n1 (clones at risk in each arm).
ccw_balance <- function(data, spec, covariates) {
  long <- run_ccw_long(data, spec)
  # Balance is assessed where artificial censoring ends: grace period or longest window
  t_star <- if (spec$design == "grace_vs_none") spec$grace else max(unlist(spec$windows))
  at_risk <- long[start <= t_star & stop > t_star, .(id, arm, w)]
  X <- stats::model.matrix(stats::reformulate(covariates), data = data)[, -1, drop = FALSE]
  # Factors get every level (the reference level too), as in a Table 1
  for (v in covariates[vapply(data[, ..covariates], is.factor, logical(1))]) {
    lv <- levels(data[[v]])
    X <- cbind(X, stats::setNames(as.data.frame(as.integer(data[[v]] == lv[1])), paste0(v, lv[1])))
  }
  X <- as.matrix(X)
  d <- merge(at_risk, data.table::data.table(id = data$id, X), by = "id")
  smd <- function(x, a, w) {
    m1 <- stats::weighted.mean(x[a == 1], w[a == 1]); m0 <- stats::weighted.mean(x[a == 0], w[a == 0])
    v1 <- sum(w[a == 1] * (x[a == 1] - m1)^2) / sum(w[a == 1])
    v0 <- sum(w[a == 0] * (x[a == 0] - m0)^2) / sum(w[a == 0])
    # Pooled SD from the unweighted samples, so both rows share one denominator
    s <- sqrt((stats::var(x[a == 1]) + stats::var(x[a == 0])) / 2)
    100 * (m1 - m0) / s
  }
  out <- data.table::rbindlist(lapply(colnames(X), function(v) {
    x <- d[[v]]
    data.table::data.table(
      term = v, weighting = c("Unweighted", "Weighted"),
      smd = c(smd(x, d$arm, rep(1, nrow(d))), smd(x, d$arm, d$w))
    )
  }))
  out[, `:=`(n0 = sum(d$arm == 0), n1 = sum(d$arm == 1), t_star = t_star)]
  out[]
}

#' Long CCW data (clones, intervals, weights) exactly as used by the pipeline
#' @param data Patient-level data; @param spec Pipeline spec.
#' @return Stacked long data (id, arm, start, stop, ev, w).
run_ccw_long <- function(data, spec) {
  switch(spec$design,
    grace_vs_none = ccw_grace_vs_none(data, spec),
    window_vs_window = ccw_window_vs_window(data, spec)
  )
}

#' Figure: standardised differences between arms before and after weighting
#' @param bal From `ccw_balance()`; @param spec Dataset config (for labels).
#' @return ggplot.
plot_balance <- function(bal, spec) {
  lab <- c(
    age_c = "Age", sex = "Women",
    perf0 = "Performance status 0", perf1 = "Performance status 1", perf2 = "Performance status 2",
    stage = "TNM stage II", deprivation = "Deprivation score",
    charlson0 = "Charlson index 0", charlson1 = "Charlson index 1", charlson2 = "Charlson index 2",
    emergency = "Emergency presentation", renal = "End-stage renal disease", nonwhite = "Non-white race"
  )
  d <- data.table::copy(bal)
  d[, term := ifelse(term %in% names(lab), lab[term], term)]
  lev <- unique(d$term)
  d[, term := factor(term, levels = rev(lev[order(match(lev, lab))]))]
  d[, weighting := factor(weighting, levels = c("Unweighted", "Weighted"))]
  ggplot2::ggplot(d, ggplot2::aes(smd, term, colour = weighting, shape = weighting)) +
    ggplot2::geom_vline(xintercept = c(-10, 10), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = c(Unweighted = "black", Weighted = "red"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(Unweighted = 1, Weighted = 16), name = NULL) +
    ggplot2::labs(
      x = "Standardised difference, surgery vs no-surgery arm (%)", y = NULL,
      subtitle = paste0("Clones at risk at day ", formatC(d$t_star[1], format = "f", digits = 0),
                        ": ", d$n1[1], " surgery, ", d$n0[1], " no surgery")
    ) +
    theme_paper() +
    ggplot2::theme(legend.position = "bottom")
}

#' Save the balance figure (PNG 600 dpi and EPS)
#' @param bal From `ccw_balance()`; @param spec Dataset config; @param dir Output folder.
#' @return File paths.
save_balance_figure <- function(bal, spec, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  p <- plot_balance(bal, spec)
  png <- file.path(dir, "fig_balance.png")
  eps <- file.path(dir, "fig_balance.eps")
  ggplot2::ggsave(png, p, width = 7, height = 5, dpi = 600, bg = "white")
  ggplot2::ggsave(eps, p, width = 7, height = 5, device = grDevices::cairo_ps)
  c(png, eps)
}

#' Render the Quarto report for one results store into the profile's output folder
#' @param qmd Report source; @param store Results store of this run; @param out_dir Output folder.
#' @param deps Targets the report reads (passed only so the report reruns when they change).
#' @return Path of the HTML report.
render_report <- function(qmd, store, out_dir, deps) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  src <- file.path(out_dir, basename(qmd))
  # Other profiles render a temporary copy, so the paper report in reports/ is never touched
  copied <- normalizePath(src, winslash = "/", mustWork = FALSE) != normalizePath(qmd, winslash = "/")
  if (copied) {
    file.copy(qmd, src, overwrite = TRUE)
    on.exit(unlink(src), add = TRUE)
  }
  # Quarto starts its own R session in the report folder, where renv is not activated: hand it
  # this session's libraries so it finds the same packages
  libs <- paste(.libPaths(), collapse = .Platform$path.sep)
  withr::with_envvar(c(R_LIBS = libs, R_LIBS_USER = libs), {
    quarto::quarto_render(
      src, execute_params = list(store = normalizePath(store, winslash = "/")), quiet = FALSE
    )
  })
  sub("\\.qmd$", ".html", src)
}
