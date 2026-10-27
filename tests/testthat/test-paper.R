test_that("baseline tables have Table S2 columns and counts that add up", {
  for (x in list(list("maringe", spec_m), list("synpuf", spec_s))) {
    raw <- read_raw_data(test_path("..", "..", x[[2]]$file), x[[1]])
    tb <- make_baseline_table(raw, x[[1]], x[[2]])
    expect_named(tb, c("label", "indent", "n1", "p1", "n0", "p0", "smd"))
    tot <- tb[label == "Total"]
    expect_equal(as.integer(tot$n1) + as.integer(tot$n0), nrow(raw))
    men <- tb[label == "Men"]; women <- tb[label == "Women"]
    expect_equal(as.integer(men$n1) + as.integer(women$n1), as.integer(tot$n1))
  }
})

test_that("Maringe table reproduces the surgery split of the data", {
  raw <- read_raw_data(test_path("..", "..", spec_m$file), "maringe")
  tb <- make_baseline_table(raw, "maringe", spec_m)
  expect_equal(as.integer(tb[label == "Total", n1]), sum(!is.na(raw$timetosurgery)))
})

test_that("bootstrap-t CI is asymmetric when T is skewed and recovers nominal coverage", {
  set.seed(2)
  # SE correlated with the estimate, as in the Maringe replicates
  e <- rnorm(4000, -0.5, 0.28)
  sim <- data.table(dataset = "m", theta = rep(c(-0.7, 0), each = 2000), rep = rep(1:2000, 2),
                    estimate = c(e[1:2000], e[2001:4000] + 0.5), status = "ok")
  sim[, se := 0.27 + 0.1 * (estimate - mean(estimate)), by = theta]
  tt <- data.table(dataset = "m", theta = c(-0.7, 0), scenario = c("observed", "null"),
                   truth = c(-0.5, 0))
  q <- boot_t_quantiles(sim[theta == -0.7, estimate], sim[theta == -0.7, se], -0.5)
  expect_lt(q[["q_lo"]], 0)
  expect_gt(q[["q_hi"]], 0)
  expect_false(isTRUE(all.equal(-q[["q_lo"]], q[["q_hi"]], tolerance = 0.02)))
  ci <- boot_t_ci(-0.5, 0.27, q)
  expect_lt(ci$lo, -0.5)
  expect_gt(ci$hi, -0.5)
  cm <- compare_ci_methods(sim, tt)
  expect_equal(nrow(cm), 10L)
  expect_setequal(unique(cm$method), c("Wald (robust SE)", "QBA widening", "Bootstrap-t"))
  bt <- cm[method == "Bootstrap-t" & check == "cross-fitted"]
  expect_true(all(abs(bt$coverage - 0.95) < 0.02))
  # Misses add up to non-coverage
  expect_equal(cm$miss_below + cm$miss_above, 1 - cm$coverage)
  expect_equal(nrow(make_ci_table(cm)), 10L)
  expect_s3_class(plot_ci_methods(cm), "ggplot")
})

test_that("balance at the grace period uses the clones still at risk and both weightings", {
  b <- ccw_balance(base_test, spec_m$pipeline, spec_m$dgm$covariates)
  expect_setequal(unique(b$weighting), c("Unweighted", "Weighted"))
  at_risk <- base_test[time > spec_m$pipeline$grace]
  expect_equal(b$n1[1], sum(at_risk$trt_time <= spec_m$pipeline$grace))
  expect_equal(b$n0[1], sum(at_risk$trt_time > spec_m$pipeline$grace))
  expect_true(all(is.finite(b$smd)))
})

test_that("results tables and figures are built from the adjustment outputs", {
  set.seed(5)
  sim <- data.table(dataset = "m", theta = rep(c(-0.7, 0), each = 400), rep = rep(1:400, 2),
                    estimate = c(rnorm(400, -0.49, 0.28), rnorm(400, -0.04, 0.28)), se = 0.27,
                    status = "ok", n_outcomes = 60)
  tt <- data.table(dataset = "m", theta = c(-0.7, 0), scenario = c("observed", "null"),
                   truth = c(-0.475, 0), truth_mcse = c(0.002, 0), psi_large_sample = c(-0.475, -0.03))
  obs <- data.table(estimate = -0.477, se = 0.265, ci_lo = -1, ci_hi = 0.04, p_value = 0.07)
  adj <- adjust_simple(obs, sim, tt, "m")
  val <- validate_adjustment(adj, copy(sim)[, estimate := estimate + rnorm(.N, 0, 0.01)], tt)
  nest_o <- data.table(estimate = rnorm(30, -0.48, 0.28), se = 0.27, bias_r = rnorm(30, -0.015, 0.01),
                       mean_inner_mcse = 0.018, truth_r_mcse = 0.01, n_inner_ok = 250, n_eval = 2,
                       at_bound = 0, theta_r = -0.7, seconds = 1, status = "ok",
                       mean_inner = rnorm(30, -0.48, 0.28))
  nest_o[, est_adj := estimate - bias_r]
  nest_o[, `:=`(ci_lo_np = est_adj - 0.55, ci_hi_np = est_adj + 0.55, ci_lo_np_refl = est_adj - 0.56,
                ci_hi_np_refl = est_adj + 0.54)]
  nest_n <- data.table(estimate = rnorm(30, 0, 0.28), se = 0.27, reject = rbinom(30, 1, 0.05),
                       n_inner_ok = 250, seconds = 1, status = "ok")
  nest_n[, reject_fav := reject * (estimate < 0)]
  ns <- summarise_nested(nest_o, nest_n, tt, adj)
  perf <- compare_performance(sim, tt, val, ns)
  expect_setequal(unique(main_analyses(perf)$analysis), names(paper_cols))
  perfm <- summarise_performance(sim, tt)
  tab <- make_results_table(obs, list(theta_hat = -0.4, gamma = c(age_c = 0.1), theta_hat = -0.4,
                                      se = c(age_c = 0.05, theta = 0.2), rho = 1.3, n_events = 48,
                                      cens = list(model = "none")),
                            list(theta_cal = -0.7), tt, perfm, adj, perf, ns,
                            list(arm1 = "surgery", outcome_label = "Death", dgm = list(age_unit = 5),
                                 pipeline = list(endpoint = "hr")))
  expect_named(tab, c("analysis", "performance", "ci_variants", "dgm"))
  expect_true(all(c("RMSE (MCSE)", "MSE = RMSE^2 (MCSE), log HR^2") %in% tab$performance$measure))
  expect_false(any(grepl("Power", tab$performance$measure)))
  expect_false("Power" %in% unique(perf$measure))
  expect_named(tab$performance, c("measure", "ccw", "simple", "nested"))
  expect_true(all(nzchar(tab$analysis$value)))
  expect_s3_class(plot_paper_performance(perf), "ggplot")
  expect_s3_class(plot_paper_distribution(sim, adj), "ggplot")
  expect_s3_class(plot_nested(nest_o, adj), "ggplot")
  # Survival-difference analysis: same tables, shown in percentage points
  spec_sd <- list(arm1 = "surgery", outcome_label = "Death", dgm = list(age_unit = 5),
                  pipeline = list(endpoint = "surv12"))
  tab2 <- make_results_table(obs, list(theta_hat = -0.4, gamma = c(age_c = 0.1), se = c(age_c = 0.05, theta = 0.2),
                                       rho = 1.3, n_events = 48, cens = list(model = "none")),
                             list(theta_cal = -0.7), tt, perfm, adj, perf, ns, spec_sd)
  expect_true(any(grepl("percentage points", tab2$analysis$item)))
  expect_true("MSE = RMSE^2 (MCSE), percentage points^2" %in% tab2$performance$measure)
})
