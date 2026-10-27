# Publication outputs for the case study (paper Section 4): baseline table in the
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

#' Outcome and censoring model parameters (supplementary table; Weibull or generalized gamma)
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
    shape <- if (identical(m$family, "gengamma")) {
      data.table::data.table(
        model = model, term = c("Generalized gamma scale (sigma)", "Generalized gamma shape (Q)"),
        hr = c(f2(m$sigma), f2(m$Q)),
        ci = c("", paste0(f2(m$Q - 1.96 * m$se[["Q"]]), " to ", f2(m$Q + 1.96 * m$se[["Q"]])))
      )
    } else {
      data.table::data.table(model = model, term = "Weibull shape (rho)", hr = f2(m$rho), ci = "")
    }
    rbind(
      shape,
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

#' Coverage of 95% CIs: Wald (robust SE), QBA widening, bootstrap-t and (optional) confidence belt
#' @param sim_results Replicates (one dataset, both scenarios).
#' @param truth_table From `make_truth_table()`; @param belt From `fit_belt()` or NULL.
#' @param level Nominal level.
#' @return data.table: method, scenario (where coverage is checked), source (scenario the
#'   CI was learned from), check, coverage, mcse, miss_below, miss_above, n.
compare_ci_methods <- function(sim_results, truth_table, belt = NULL, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  d <- merge(sim_results[!is.na(estimate)], truth_table[, .(dataset, theta, scenario, truth)],
             by = c("dataset", "theta"))
  d[, fold := rep %% 2L]
  summ <- function(lo, hi, truth) {
    hit <- lo <= truth & truth <= hi
    m <- mean(hit)
    data.table::data.table(
      coverage = m, mcse = sqrt(m * (1 - m) / length(hit)),
      # Misses split by side: a symmetric CI should miss 2.5% on each
      miss_below = mean(hi < truth), miss_above = mean(lo > truth), n = length(hit)
    )
  }
  # Learned CI: fitted on one half of the replicates of `src`, checked on the other half of `tg`
  learned <- function(method, src, tg) {
    res <- data.table::rbindlist(lapply(0:1, function(k) {
      s <- d[fold != k & scenario == src]
      x <- d[fold == k & scenario == tg]
      if (method == "QBA widening") {
        b <- mean(s$estimate - s$truth)
        k_se <- max(1, stats::sd(s$estimate) / sqrt(mean(s$se^2)))
        data.table::data.table(lo = x$estimate - b - z * k_se * x$se,
                               hi = x$estimate - b + z * k_se * x$se, truth = x$truth)
      } else {
        ci <- boot_t_ci(x$estimate, x$se, boot_t_quantiles(s$estimate, s$se, s$truth[1], level))
        data.table::data.table(lo = ci$lo, hi = ci$hi, truth = x$truth)
      }
    }))
    # Same scenario is near-circular; transfer tests whether the CI holds at other true effects
    cbind(data.table::data.table(method = method, scenario = tg, source = src,
                                 check = if (src == tg) "cross-fitted" else "transfer"),
          summ(res$lo, res$hi, res$truth))
  }
  wald <- lapply(c("observed", "null"), function(tg) {
    x <- d[scenario == tg]
    cbind(data.table::data.table(method = "Wald (robust SE)", scenario = tg, source = NA_character_,
                                 check = "none"),
          summ(x$estimate - z * x$se, x$estimate + z * x$se, x$truth))
  })
  grid <- data.table::CJ(method = c("QBA widening", "Bootstrap-t"), src = c("observed", "null"),
                         tg = c("observed", "null"), sorted = FALSE)
  out <- rbind(data.table::rbindlist(wald),
               data.table::rbindlist(Map(learned, grid$method, grid$src, grid$tg)))
  if (is.null(belt)) return(out)
  rbind(out, belt_coverage(sim_results, truth_table, belt), fill = TRUE)
}

#' Table of CI coverage by method (estimate (MCSE), misses by side)
#' @param ci From `compare_ci_methods()`.
#' @return data.table ready for printing.
make_ci_table <- function(ci) {
  pc <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
  lab <- c(observed = "True HR = observed", null = "No effect")
  learned_from <- data.table::fifelse(
    is.na(ci$source), "Not learned",
    data.table::fifelse(
      ci$source == "grid", "Grid without this effect (cross-fitted)",
      paste0(lab[ci$source], data.table::fifelse(ci$check == "cross-fitted", " (cross-fitted)", " (transfer)"))
    )
  )
  data.table::data.table(
    method = ci$method, checked_in = unname(lab[ci$scenario]), learned_from = learned_from,
    coverage = paste0(pc(ci$coverage), " (", pc(ci$mcse), ")"),
    below = pc(ci$miss_below), above = pc(ci$miss_above)
  )[order(match(method, names(ci_cols)), checked_in)]
}

ci_cols <- c("Wald (robust SE)" = "black", "QBA widening" = "red", "Bootstrap-t" = "#1f6fb4",
             "Confidence belt" = "#2e8b57")

#' Figure: coverage and one-sided misses of the three CIs (points +- 1.96 MCSE)
#' @param ci From `compare_ci_methods()`.
#' @return ggplot.
plot_ci_methods <- function(ci) {
  d <- data.table::copy(ci)
  # Belt rows carry bootstrap MCSEs that include the grid's Monte Carlo error; others binomial
  for (v in c("mcse_below", "mcse_above")) if (!v %in% names(d)) d[, (v) := NA_real_]
  mcse_side <- function(p, n, boot) data.table::fcoalesce(boot, sqrt(p * (1 - p) / n))
  lv <- c("Not learned", "Same scenario (cross-fitted)", "Other scenario (transfer)",
          "Grid without this effect (cross-fitted)")
  d[, learned := factor(data.table::fcase(
    check == "none", lv[1], source %in% "grid", lv[4], check == "transfer", lv[3], default = lv[2]
  ), levels = lv)]
  d[, scenario := factor(scenario, c("observed", "null"),
                         c("Checked at true HR = observed", "Checked at no effect"))]
  long <- data.table::rbindlist(list(
    d[, .(method, scenario, learned, measure = "Coverage", value = coverage, mcse, ref = 0.95)],
    d[, .(method, scenario, learned, measure = "CI entirely above truth", value = miss_above,
          mcse = mcse_side(miss_above, n, mcse_above), ref = 0.025)],
    d[, .(method, scenario, learned, measure = "CI entirely below truth", value = miss_below,
          mcse = mcse_side(miss_below, n, mcse_below), ref = 0.025)]
  ))
  long[, measure := factor(measure, c("Coverage", "CI entirely above truth", "CI entirely below truth"))]
  long[, method := factor(method, intersect(names(ci_cols), unique(method)))]
  refs <- unique(long[, .(measure, ref)])
  pd <- ggplot2::position_dodge(width = 0.6)
  ggplot2::ggplot(long, ggplot2::aes(method, 100 * value, colour = method, shape = learned)) +
    ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = 100 * ref), linetype = "dashed") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = 100 * (value - 1.96 * mcse),
                                         ymax = 100 * (value + 1.96 * mcse), group = learned),
                            position = pd, linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(group = learned), position = pd, size = 2.6) +
    ggplot2::facet_grid(measure ~ scenario, scales = "free_y") +
    ggplot2::scale_colour_manual(values = ci_cols, guide = "none") +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 18), name = "CI learned from:") +
    ggplot2::scale_x_discrete(labels = function(x) sub(" ", "\n", sub(" \\(robust SE\\)", "", x))) +
    ggplot2::labs(x = NULL, y = "% of replicates (95% Monte Carlo interval)") +
    theme_paper(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

#' Save the CI comparison figure (PNG and EPS, like the other paper figures)
#' @param ci From `compare_ci_methods()`; @param dir Output folder.
#' @return File paths.
save_ci_figure <- function(ci, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  p <- plot_ci_methods(ci)
  png <- file.path(dir, "fig_ci_methods.png")
  eps <- file.path(dir, "fig_ci_methods.eps")
  ggplot2::ggsave(png, p, width = 11, height = 7.5, dpi = 600, bg = "white")
  ggplot2::ggsave(eps, p, width = 11, height = 7.5, device = grDevices::cairo_ps)
  c(png, eps)
}

#' Display scale of an endpoint: log HR is shown as HR; a survival difference in percentage points
#' @param spec Dataset config. @return list(hr, k, unit, est, ci) with formatters.
display_scale <- function(spec) {
  ep <- endpoint_info(spec$pipeline)
  f <- function(x, d) formatC(x, format = "f", digits = d)
  if (ep$endpoint == "hr") {
    list(
      hr = TRUE, k = 1, unit = "log HR", name = "HR",
      est = function(e) f(exp(e), 4),
      ci = function(e, lo, hi) paste0(f(exp(e), 4), " (", f(exp(lo), 4), " to ", f(exp(hi), 4), ")")
    )
  } else {
    list(
      hr = FALSE, k = 100, unit = "percentage points", name = "survival difference",
      est = function(e) f(100 * e, 2),
      ci = function(e, lo, hi) paste0(f(100 * e, 2), " (", f(100 * lo, 2), " to ", f(100 * hi, 2), ")")
    )
  }
}

#' Observed analysis, calibration, adjustments and performance in one list of tables
#' @param obs Observed pipeline result; @param dgm Fitted DGM.
#' @param cal From `calibrate_observed()`.
#' @param truth_table,performance,adj,perf,nested Pipeline targets for one dataset
#'   (`perf` from `compare_performance()`, `nested` from `summarise_nested()`).
#' @param spec Dataset config (`arm1` names the strategy coded 1).
#' @return list(analysis, performance, ci_variants, dgm), each a data.table. Estimates, biases, SEs
#'   and RMSE are shown on the display scale of the endpoint (HR, or percentage points for the
#'   survival difference); MSE = RMSE^2 in the squared unit of the log HR or of percentage points.
make_results_table <- function(obs, dgm, cal, truth_table, performance, adj, perf, nested, spec) {
  ds <- display_scale(spec)
  f2 <- function(x) formatC(x, format = "f", digits = 2)
  f3 <- function(x) formatC(x, format = "f", digits = 3)
  f4 <- function(x) formatC(x, format = "f", digits = 4)
  pc <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
  # Log-scale quantities (log HR) keep their natural scale; survival differences are in points
  sc <- function(x, d = 3) formatC(ds$k * x, format = "f", digits = d)
  est_mc <- function(v, mc) paste0(sc(v), " (", sc(mc), ")")
  po <- performance[scenario == "observed"]
  pn <- performance[scenario == "null"]
  tn <- truth_table[scenario == "null"]
  ns <- nested$summary
  nm <- ds$name
  lim <- function(e) if (ds$hr) f4(exp(e)) else f2(100 * e)
  analysis <- data.table::data.table(
    item = c(
      paste0("Observed CCW ", nm, " (95% CI)"), "Observed p-value",
      "Prediction model: on-treatment HR (exp(theta-hat))",
      "Calibrated generating HR, exp(theta*), with mean CCW estimate = observed",
      paste0("Mean CCW estimate at theta* over the simulated studies (MCSE, ", ds$unit, ")"),
      paste0("True CCW value at theta*, psi(theta*) (MCSE of psi, ", ds$unit, ")"),
      paste0("Bias at theta*, ", ds$unit, " (MCSE)"),
      paste0("Large-sample CCW value with no treatment effect, psi(0)"),
      paste0("QBA-adjusted ", nm, ", simple adjustment (non-parametric 95% CI)"),
      paste0("QBA-adjusted ", nm, ", nested simulation (non-parametric 95% CI)"),
      paste0("QBA-adjusted ", nm, ", parametric 95% CI (CI inflation factor)"),
      paste0("Non-parametric 5% test: H0 rejected if CCW estimate is below / above"),
      "Simulation-calibrated p-value (non-parametric)",
      "Nested: SD of the estimated bias correction across outer datasets (within MC error)"
    ),
    value = c(
      ds$ci(obs$estimate, obs$ci_lo, obs$ci_hi),
      formatC(obs$p_value, format = "f", digits = 3),
      f4(exp(dgm$theta_hat)), f4(exp(cal$theta_cal)),
      paste0(lim(adj$mean_est), " (", f3(ds$k * adj$mean_est_mcse), ")"),
      paste0(lim(adj$truth), " (", f3(ds$k * adj$truth_mcse), ")"),
      est_mc(adj$bias, adj$bias_mcse),
      lim(tn$psi_large_sample),
      ds$ci(adj$est_adj, adj$ci_lo_np, adj$ci_hi_np),
      ds$ci(adj$est_adj, ns$ci_lo_np, ns$ci_hi_np),
      paste0(ds$ci(adj$est_adj, adj$ci_lo_par, adj$ci_hi_par), "; k = ", f2(adj$se_inflation)),
      paste0(lim(adj$null_mean + adj$null_q_lo), " / ", lim(adj$null_mean + adj$null_q_hi)),
      formatC(adj$p_np, format = "f", digits = 3),
      paste0(sc(ns$bias_r_sd), " (", sc(ns$bias_r_within_mcse), ")")
    )
  )
  cols <- list(
    ccw = c("CCW", "Wald (robust SE)"),
    simple = c("Simple adjustment", "Non-parametric"),
    nested = c("Nested adjustment", "Non-parametric")
  )
  get <- function(col, m) {
    r <- perf[analysis == cols[[col]][1] & method == cols[[col]][2] & measure == m]
    if (nrow(r) == 0) NULL else r
  }
  cell <- function(m, pct, d = 1) {
    vapply(names(cols), function(col) {
      r <- get(col, m)
      if (is.null(r)) return("")
      if (pct) paste0(pc(r$value, d), " (", pc(r$mcse, d), ")") else est_mc(r$value, r$mcse)
    }, character(1))
  }
  plain <- function(fun) vapply(names(cols), function(col) fun(col), character(1))
  rel_bias <- plain(function(col) {
    pc(get(col, "Bias")$value / get(col, "Empirical SE")$value)
  })
  misses <- plain(function(col) {
    paste0(pc(get(col, "CI entirely below truth")$value), " / ",
           pc(get(col, "CI entirely above truth")$value))
  })
  favour <- plain(function(col) {
    r <- get(col, "Rejections favouring arm 1")
    if (is.null(r)) "" else pc(r$value, 3)
  })
  # RMSE and MSE with MCSE on the display scale (MSE scales with k^2)
  rmse_mc <- plain(function(col) {
    r <- get(col, "RMSE")
    paste0(sc(r$value), " (", sc(r$mcse), ")")
  })
  mse_mc <- plain(function(col) {
    r <- get(col, "MSE")
    paste0(formatC(ds$k^2 * r$value, format = "f", digits = if (ds$hr) 4 else 2), " (",
           formatC(ds$k^2 * r$mcse, format = "f", digits = if (ds$hr) 4 else 2), ")")
  })
  n_sims <- c(
    ccw = paste0(format(po$B, big.mark = ","), " per scenario (", po$n_fail + pn$n_fail,
                 " failed)"),
    simple = paste0(format(get("simple", "95% CI coverage")$n, big.mark = ","), " per scenario"),
    nested = paste0(ns$n_outer, " outer x ", ns$n_inner, " inner (", ns$n_fail + ns$n_null_fail,
                    " failed)")
  )
  rows <- list(
    c(paste0("Mean estimate, ", ds$unit, " (true value)"),
      plain(function(col) {
        r <- get(col, "Mean estimate")
        paste0(sc(r$value), " (", sc(r$ref), ")")
      })),
    c(paste0("Bias, ", ds$unit, " (MCSE)"), cell("Bias", FALSE)),
    c("Relative bias, % of empirical SE", rel_bias),
    c("Empirical SE (MCSE)", cell("Empirical SE", FALSE)),
    c("Model-based (robust) SE", sc(po$mod_se, 3), "", ""),
    c("Relative error in model SE, %", pc(po$rel_err_modse), "", ""),
    c("RMSE (MCSE)", rmse_mc),
    c(paste0("MSE = RMSE^2 (MCSE), ", ds$unit, "^2"), mse_mc),
    c("95% CI coverage, % (MCSE)", cell("95% CI coverage", TRUE)),
    c("CI entirely below / above the true value, %", misses),
    c("Type I error, % (no effect)", cell("Type I error", TRUE, 3)),
    c(paste0("Rejections favouring ", spec$arm1, ", % (no effect)"), favour),
    c("Simulated studies", n_sims)
  )
  performance_tab <- data.table::as.data.table(do.call(rbind, rows))
  data.table::setnames(performance_tab, c("measure", "ccw", "simple", "nested"))
  list(analysis = analysis, performance = performance_tab,
       ci_variants = make_ci_variant_table(perf), dgm = make_dgm_table(dgm, spec))
}

#' Coverage of the CI variants (non-parametric as in the paper, reflected, parametric, Wald)
#' @param perf From `compare_performance()`.
#' @return data.table: analysis, method, coverage (MCSE), below, above (%).
make_ci_variant_table <- function(perf) {
  pc <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
  w <- data.table::dcast(perf[grepl("CI", measure)], analysis + method ~ measure,
                         value.var = c("value", "mcse"))
  ord <- c("Wald (robust SE)", "Wald (robust SE), validation set", "Non-parametric",
           "Non-parametric (reflected)", "Parametric")
  an <- c("CCW", "Simple adjustment", "Nested adjustment")
  w <- w[order(match(analysis, an), match(method, ord))]
  w[, .(
    analysis, method,
    coverage = paste0(pc(`value_95% CI coverage`), " (", pc(`mcse_95% CI coverage`), ")"),
    below = pc(`value_CI entirely below truth`), above = pc(`value_CI entirely above truth`)
  )]
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

paper_cols <- c("CCW" = "black", "Simple adjustment" = "red", "Nested adjustment" = "#1f6fb4")

#' Rows of the three analyses compared in the paper (CCW, simple and nested non-parametric)
#' @param perf From `compare_performance()`.
#' @return Subset of `perf`.
main_analyses <- function(perf) {
  perf[(analysis == "CCW" & method == "Wald (robust SE)") |
         (analysis %in% c("Simple adjustment", "Nested adjustment") & method == "Non-parametric")]
}

#' Figure: bias, RMSE, MSE, coverage and type I error of CCW, simple and nested adjustment
#' @param perf From `compare_performance()`.
#' @return ggplot (points +- 1.96 MCSE).
plot_paper_performance <- function(perf) {
  ms <- c("Bias", "RMSE", "MSE", "95% CI coverage", "Type I error")
  d <- main_analyses(perf)[measure %in% ms]
  d[, measure := factor(measure, levels = ms)]
  d[, analysis := factor(analysis, levels = names(paper_cols))]
  refs <- unique(d[!is.na(ref), .(measure, ref)])
  ggplot2::ggplot(d, ggplot2::aes(analysis, value, colour = analysis)) +
    ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = ref), linetype = "dashed") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = value - 1.96 * mcse, ymax = value + 1.96 * mcse),
                            linewidth = 1.1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~measure, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = paper_cols, name = "Method:") +
    ggplot2::scale_x_discrete(labels = function(x) sub(" adjustment", "", x), drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Estimate with 95% Monte Carlo interval") +
    theme_paper()
}

#' Figure: simulated estimates at theta* and with no effect, with true values and test limits
#' @param sim_results Replicates; @param adj From `adjust_simple()`.
#' @return ggplot.
plot_paper_distribution <- function(sim_results, adj, xlab = "Estimated log HR") {
  # Unicode escape keeps the source ASCII; the theta panel is the calibrated scenario
  p1 <- "Generating effect \u03b8*"
  p2 <- "No treatment effect"
  d <- rbind(
    sim_results[theta == adj$theta_cal & !is.na(estimate), .(x = estimate, panel = p1)],
    sim_results[theta == 0 & !is.na(estimate), .(x = estimate, panel = p2)]
  )
  d[, panel := factor(panel, levels = c(p1, p2))]
  lines <- data.table::data.table(
    panel = c(p1, p1, p2, p2, p2),
    x = c(adj$truth, adj$mean_est, 0, adj$null_mean + adj$null_q_lo, adj$null_mean + adj$null_q_hi),
    what = c("True value", "Mean of estimates", "True value", rep("Rejection limits (5% test)", 2))
  )
  lines[, panel := factor(panel, levels = levels(d$panel))]
  vals <- c("True value" = "black", "Mean of estimates" = "red",
            "Rejection limits (5% test)" = "#1f6fb4")
  ggplot2::ggplot(d, ggplot2::aes(x)) +
    ggplot2::geom_histogram(bins = 60, fill = "grey70", colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = x, colour = what, linetype = what),
                        linewidth = 0.9) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::scale_colour_manual(values = vals, name = "Reference:") +
    ggplot2::scale_linetype_manual(
      values = c("True value" = "solid", "Mean of estimates" = "dashed",
                 "Rejection limits (5% test)" = "dotted"),
      name = "Reference:"
    ) +
    ggplot2::labs(x = xlab, y = "Number of simulated studies") +
    theme_paper()
}

#' Figure: nested simulation, estimated bias correction and adjusted estimate per outer dataset
#' @param nested_obs Outer replicates at theta*; @param adj From `adjust_simple()`.
#' @param xlab Axis label of the estimate scale.
#' @return ggplot.
plot_nested <- function(nested_obs, adj, xlab = "log HR") {
  p1 <- "Estimated bias correction"
  p2 <- "Bias-adjusted estimate"
  ref <- "Simple adjustment (left) / true value (right)"
  avg <- "Mean over outer datasets"
  o <- nested_obs[!is.na(est_adj)]
  d <- rbind(o[, .(x = bias_r, panel = p1)], o[, .(x = est_adj, panel = p2)])
  d[, panel := factor(panel, levels = c(p1, p2))]
  lines <- data.table::data.table(
    panel = c(p1, p1, p2, p2),
    x = c(adj$bias, mean(o$bias_r), adj$truth, mean(o$est_adj)),
    what = c(ref, avg, ref, avg)
  )
  lines[, panel := factor(panel, levels = levels(d$panel))]
  ggplot2::ggplot(d, ggplot2::aes(x)) +
    ggplot2::geom_histogram(bins = 40, fill = "grey70", colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = x, colour = what, linetype = what),
                        linewidth = 0.9) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::scale_colour_manual(values = stats::setNames(c("black", "red"), c(ref, avg)),
                                 name = "Reference:") +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "dashed"), c(ref, avg)),
                                   name = "Reference:") +
    ggplot2::labs(x = xlab, y = "Number of outer datasets") +
    theme_paper()
}

#' Marginal survival implied by the fitted PH model (Weibull or generalized gamma baseline)
#' @param m Fitted model list; @param data Base data; @param times Time grid.
#' @return Numeric vector: mean over patients of S_i(t).
dgm_marginal_surv <- function(m, data, times) {
  X <- dgm_design(data, m$covariates)
  eh <- exp(m$alpha + drop(X %*% m$gamma))
  h0 <- if (identical(m$family, "gengamma")) {
    function(t) gg_cumhaz0(t / m$tau, m$sigma, m$Q)
  } else {
    function(t) t^m$rho
  }
  vapply(times, function(t) {
    pre <- pmin(data$s, t)
    H <- eh * (h0(pre) + (data$s < t) * exp(m$theta_hat) * (h0(t) - h0(pre)))
    mean(exp(-H))
  }, numeric(1))
}

#' @rdname dgm_marginal_surv
weibull_marginal_surv <- dgm_marginal_surv

#' Figure: fit of the outcome and censoring models (Kaplan-Meier vs model)
#' @param data Base data; @param dgm Fitted DGM; @param spec Dataset config (outcome label).
#' @return ggplot.
plot_paper_fit <- function(data, dgm, spec) {
  times <- seq(0, max(data$time), length.out = 200)
  lab_m <- if (identical(dgm$family, "gengamma")) "Generalized gamma model" else "Weibull model"
  km <- function(status, panel) {
    f <- survival::survfit(survival::Surv(data$time, status) ~ 1)
    data.table::data.table(time = c(0, f$time), surv = c(1, f$surv), what = "Kaplan-Meier", panel = panel)
  }
  wb <- function(m, panel) {
    data.table::data.table(time = times, surv = dgm_marginal_surv(m, data, times),
                           what = lab_m, panel = panel)
  }
  p1 <- spec$outcome_label
  p2 <- "Loss to follow-up"
  d <- rbind(km(data$status, p1), wb(dgm, p1))
  if (identical(dgm$cens$model, "weibull")) d <- rbind(d, km(data$cens_event, p2), wb(dgm$cens, p2))
  d[, panel := factor(panel, levels = c(p1, p2))]
  ggplot2::ggplot(d, ggplot2::aes(time, surv, colour = what, linetype = what)) +
    ggplot2::geom_step(data = d[what == "Kaplan-Meier"], linewidth = 0.9) +
    ggplot2::geom_line(data = d[what == lab_m], linewidth = 0.9) +
    ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::scale_colour_manual(values = stats::setNames(c("black", "red"), c("Kaplan-Meier", lab_m)), name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "dashed"), c("Kaplan-Meier", lab_m)), name = NULL) +
    ggplot2::labs(x = "Days since time zero", y = "Proportion free of the event") +
    theme_paper() +
    ggplot2::theme(legend.position = "bottom")
}

#' Save all paper figures as PNG (600 dpi) and EPS
#' @param perf From `compare_performance()`; @param sim_results Main replicates.
#' @param nested_obs Outer replicates at theta*; @param data Base data; @param dgm Fitted DGM.
#' @param adj From `adjust_simple()`; @param spec Dataset config; @param dir Output folder.
#' @return Character vector of file paths.
save_paper_figures <- function(perf, sim_results, nested_obs, data, dgm, adj, spec, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  xl <- endpoint_info(spec$pipeline)
  xl <- if (xl$endpoint == "hr") "log HR" else "Difference in 12-month survival"
  figs <- list(
    fig_performance = list(plot_paper_performance(perf), 15, 4.4),
    fig_distribution = list(plot_paper_distribution(sim_results, adj, paste("Estimated", xl)),
                            12, 4.6),
    fig_nested = list(plot_nested(nested_obs, adj, xl), 12, 4.6),
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
  # Render in a temporary folder and copy only the HTML back: in a synced project folder
  # (Dropbox) the sync client locks Quarto's intermediate files and the cleanup fails
  tmp <- tempfile("qba_report_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  src <- file.path(tmp, basename(qmd))
  file.copy(qmd, src)
  # Quarto starts its own R session, where renv is not activated: hand it this session's
  # libraries so it finds the same packages
  libs <- paste(.libPaths(), collapse = .Platform$path.sep)
  log <- file.path(out_dir, "render_error.log")
  unlink(log)
  params <- list(store = normalizePath(store, winslash = "/"),
                 root = normalizePath(".", winslash = "/"))
  withr::with_envvar(c(R_LIBS = libs, R_LIBS_USER = libs), tryCatch(
    quarto::quarto_render(src, execute_params = params, quiet = FALSE),
    error = function(e) {
      # targets shortens error messages, which cuts off the R error inside the Quarto output
      msg <- c(conditionMessage(e), e$stderr)
      p <- e$parent
      while (!is.null(p)) {
        msg <- c(msg, "", conditionMessage(p), p$stderr, p$stdout)
        p <- p$parent
      }
      writeLines(as.character(unlist(msg)), log)
      stop("Quarto rendering failed; full output in ", log, call. = FALSE)
    }
  ))
  html <- sub("\\.qmd$", ".html", basename(qmd))
  out <- file.path(out_dir, html)
  if (!file.copy(file.path(tmp, html), out, overwrite = TRUE)) {
    cli::cli_abort("Could not write {out} (is it open or locked?).")
  }
  out
}
