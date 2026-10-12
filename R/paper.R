# Publication outputs for the SynPUF case study (paper Section 4): baseline table in the
# style of the Maringe Table S2, results tables, performance after QBA adjustment, figures
# in the style of the earlier simulation papers (grey strips, framed panels, MCSE bars).

#' Standardised mean difference (%) for a binary indicator or a continuous variable
#' @param x Values; @param g Logical group (TRUE = first column).
#' @return Scalar SMD in percent.
smd_pct <- function(x, g) {
  m1 <- mean(x[g]); m0 <- mean(x[!g])
  v <- if (all(x %in% 0:1)) (m1 * (1 - m1) + m0 * (1 - m0)) / 2 else (stats::var(x[g]) + stats::var(x[!g])) / 2
  100 * (m1 - m0) / sqrt(v)
}

#' Baseline characteristics by observed treatment initiation (Table S2 layout)
#' @param raw Raw SynPUF cohort (`baseccwcohort.sas7bdat`).
#' @param spec `cfg$datasets$synpuf` (uses `horizon`).
#' @return data.table: label, indent, n1, p1, n0, p0, smd (text columns, ready to print).
make_baseline_table <- function(raw, spec) {
  d <- data.table::as.data.table(raw)
  h <- spec$horizon
  # Initiation counts only within the 180-day follow-up and before the outcome or loss
  fu <- pmin(d$followup, h)
  g <- !is.na(d$RX_start) & d$RX_start <= fu
  n1 <- sum(g); n0 <- sum(!g)
  f1 <- function(x) formatC(x, format = "f", digits = 1)
  cat_row <- function(label, x) {
    data.table::data.table(
      label = label, indent = 1L,
      n1 = format(sum(x[g])), p1 = f1(100 * mean(x[g])),
      n0 = format(sum(x[!g])), p0 = f1(100 * mean(x[!g])),
      smd = f1(smd_pct(as.numeric(x), g))
    )
  }
  # Binary variables get one SMD (on the first level), as in Table S2
  bin_rows <- function(label_a, label_b, x) {
    r <- rbind(cat_row(label_a, x), cat_row(label_b, !x))
    r[2, smd := ""]
    r
  }
  head_row <- function(label) {
    data.table::data.table(label = label, indent = 0L, n1 = "", p1 = "", n0 = "", p0 = "", smd = "")
  }
  q_row <- function(label, x, smd = TRUE) {
    q <- function(v) formatC(stats::quantile(v, c(0.5, 0.25, 0.75), names = FALSE), format = "f", digits = 0)
    a <- q(x[g]); b <- if (any(!g & !is.na(x))) q(x[!g]) else c("", "", "")
    data.table::data.table(
      label = label, indent = 0L, n1 = a[1], p1 = paste(a[2], a[3]),
      n0 = b[1], p0 = paste(b[2], b[3]),
      smd = if (smd) f1(smd_pct(x, g)) else ""
    )
  }
  race <- d$race
  rbind(
    data.table::data.table(label = "Total", indent = 0L, n1 = format(n1), p1 = "100.0",
                           n0 = format(n0), p0 = "100.0", smd = ""),
    head_row("Sex"),
    bin_rows("Men", "Women", d$sex == "1"),
    head_row("Race"),
    cat_row("White", race == "1"), cat_row("Black", race == "2"),
    cat_row("Other", race == "3"), cat_row("Hispanic", race == "5"),
    head_row("End-stage renal disease"),
    bin_rows("Yes", "No", d$renal == "Y"),
    head_row("Outcome within 180 days"),
    cat_row("Readmission or death", d$outcome == 1 & d$followup <= h)[, smd := ""],
    cat_row("Lost to follow-up", d$outcome == 0 & d$followup < h)[, smd := ""],
    data.table::data.table(label = "", indent = 0L, n1 = "Median", p1 = "Q1 Q3",
                           n0 = "Median", p0 = "Q1 Q3", smd = ""),
    q_row("Age at MI (years)", d$baseline_age),
    q_row("Time to treatment start (days)", data.table::fifelse(g, d$RX_start, NA_real_), smd = FALSE)
  )
}

#' Weibull outcome and censoring model parameters (supplementary table)
#' @param dgm Fitted DGM with `cens` element.
#' @return data.table: model, term, hr, ci (text), n_events.
make_dgm_table <- function(dgm) {
  lab <- c(age_c = "Age (per 10 years)", sex = "Female", renal = "End-stage renal disease",
           nonwhite = "Non-white race", theta = "On treatment (time-varying)")
  one <- function(m, model) {
    terms <- c(names(m$gamma), "theta")
    b <- c(m$gamma, theta = m$theta_hat)
    se <- m$se[terms]
    f2 <- function(x) formatC(x, format = "f", digits = 2)
    rbind(
      data.table::data.table(model = model, term = "Weibull shape (rho)", hr = f2(m$rho), ci = ""),
      data.table::data.table(
        model = model, term = unname(lab[terms]), hr = f2(exp(b)),
        ci = paste0(f2(exp(b - 1.96 * se)), " to ", f2(exp(b + 1.96 * se)))
      )
    )[, n_events := m$n_events]
  }
  out <- one(dgm, "Outcome: readmission or death")
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
#' @return list(analysis = data.table, performance = data.table, dgm = data.table).
make_results_table <- function(obs, dgm, cal, truth_table, performance, correction, perf_adj) {
  f2 <- function(x) formatC(x, format = "f", digits = 2)
  f3 <- function(x) formatC(x, format = "f", digits = 3)
  pc <- function(x) formatC(100 * x, format = "f", digits = 1)
  hr_ci <- function(e, lo, hi) paste0(f2(exp(e)), " (", f2(exp(lo)), " to ", f2(exp(hi)), ")")
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
      "Large-sample CCW HR with no treatment effect, exp(psi(0)) (MCSE of psi)",
      "QBA-adjusted HR (95% CI)", "CI inflation factor (empirical / model SE)",
      "Calibrated critical value |z| for a 5% test", "Simulation-calibrated p-value"
    ),
    value = c(
      hr_ci(obs$estimate, obs$ci_lo, obs$ci_hi),
      formatC(obs$p_value, format = "g", digits = 2),
      f2(exp(dgm$theta_hat)), f2(exp(cal$theta_star)),
      paste0(f2(exp(po$truth)), " (", f3(po$truth_mcse), ")"),
      paste0(f2(exp(tn$psi_large_sample)), " (", f3(pn$truth_mcse), ")"),
      hr_ci(cr$est_corr, cr$ci_lo_corr, cr$ci_hi_corr),
      f2(cr$se_inflation), f2(cr$z_crit_calibrated),
      formatC(cr$p_calibrated, format = "g", digits = 2)
    )
  )
  # Unadjusted CCW performance with full MCSE (incl. truth MC error) plus cross-fitted adjusted
  pa <- perf_adj[analysis == "CCW + QBA adjustment"]
  getp <- function(m) pa[measure == m]
  fmt <- function(v, mc, pct) if (pct) paste0(pc(v), " (", pc(mc), ")") else paste0(f3(v), " (", f3(mc), ")")
  performance_tab <- data.table::data.table(
    measure = c(
      "Mean estimate, log HR (true value)", "Bias, log HR", "Relative bias, % of empirical SE",
      "Empirical SE", "Model-based (robust) SE", "Relative error in model SE, %", "RMSE",
      "95% CI coverage, %", "Power, % (true HR = observed)", "Type I error, % (no effect)",
      "Rejections favouring the longer window, % (no effect)", "Failed replicates, n"
    ),
    ccw = c(
      paste0(f3(po$mean_est), " (", f3(po$truth), ")"),
      fmt(po$bias, po$bias_mcse_total, FALSE),
      pc(po$std_bias),
      fmt(po$emp_se, po$emp_se_mcse, FALSE), f3(po$mod_se), pc(po$rel_err_modse), f3(po$rmse),
      fmt(po$coverage, po$coverage_mcse, TRUE),
      fmt(po$reject_0, po$reject_0_mcse, TRUE),
      fmt(pn$reject_0, pn$reject_0_mcse, TRUE),
      pc(pn$reject_0_lower),
      format(po$n_fail + pn$n_fail)
    ),
    adjusted = c(
      "", fmt(getp("Bias (log HR)")$value, getp("Bias (log HR)")$mcse, FALSE), "", "", "", "", "",
      fmt(getp("95% CI coverage")$value, getp("95% CI coverage")$mcse, TRUE),
      fmt(getp("Power")$value, getp("Power")$mcse, TRUE),
      fmt(getp("Type I error")$value, getp("Type I error")$mcse, TRUE), "", ""
    )
  )
  list(analysis = analysis, performance = performance_tab, dgm = make_dgm_table(dgm))
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
    ggplot2::labs(x = NULL, y = "Estimate (± 1.96 MCSE)") +
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
#' @param data Base data; @param dgm Fitted DGM.
#' @return ggplot.
plot_paper_fit <- function(data, dgm) {
  times <- seq(0, max(data$time), length.out = 200)
  km <- function(status, panel) {
    f <- survival::survfit(survival::Surv(data$time, status) ~ 1)
    data.table::data.table(time = c(0, f$time), surv = c(1, f$surv), what = "Kaplan-Meier", panel = panel)
  }
  wb <- function(m, panel) {
    data.table::data.table(time = times, surv = weibull_marginal_surv(m, data, times),
                           what = "Weibull model", panel = panel)
  }
  p1 <- "Readmission or death"
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
    ggplot2::labs(x = "Days since MI discharge", y = "Proportion free of the event") +
    theme_paper() +
    ggplot2::theme(legend.position = "bottom")
}

#' Save all paper figures as PNG (600 dpi) and EPS
#' @param performance,perf_adj,sim_results,data,dgm Pipeline targets (one dataset).
#' @param dir Output folder; @param z_crit Calibrated critical value.
#' @return Character vector of file paths.
save_paper_figures <- function(performance, perf_adj, sim_results, data, dgm, dir, z_crit) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  figs <- list(
    fig_performance = list(plot_paper_performance(perf_adj), 13, 4.2),
    fig_distribution = list(plot_paper_distribution(sim_results, performance, z_crit), 12, 4.6),
    fig_model_fit = list(plot_paper_fit(data, dgm), 10, 4.6)
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
