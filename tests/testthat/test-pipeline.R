contract <- c("estimate", "se", "ci_lo", "ci_hi", "p_value", "n_events", "status", "message")

test_that("both pipelines return the contract columns on the real data", {
  for (x in list(list(base_test, spec_m), list(base_syn, spec_s))) {
    res <- run_pipeline(x[[1]], x[[2]]$pipeline)
    expect_named(res, contract)
    expect_true(is.finite(res$estimate))
  }
})

test_that("Maringe cloning keeps two clones per patient and censors deviators", {
  cl <- build_clones(base_test, spec_m$pipeline$grace)
  expect_equal(nrow(cl), 2 * nrow(base_test))
  op <- cl[arm == 0L & is.finite(trt_time)]
  expect_true(all(op$art_cens == 1L & op$end == op$trt_time))
})

test_that("tutorial design: non-starters followed past the window are censored at it", {
  long <- ccw_window_arm(base_syn, 30, spec_s$pipeline)
  ns <- base_syn[time > 30 & trt_time > 30, id]
  expect_true(all(long[id %in% ns, max(stop), by = id]$V1 == 30))
  expect_true(all(long[start < 30, w] == 1))
})

test_that("results do not depend on batch size", {
  dgm <- fit_dgm(base_test, spec_m$dgm)
  seeds <- make_crn_seeds(1, 4)
  r2 <- rbind(
    run_batch(base_test, dgm, -0.5, seeds, 1, spec_m$pipeline, 2),
    run_batch(base_test, dgm, -0.5, seeds, 2, spec_m$pipeline, 2)
  )
  r4 <- run_batch(base_test, dgm, -0.5, seeds, 1, spec_m$pipeline, 4)
  expect_equal(r2$estimate, r4$estimate)
})

test_that("grace period is 180 days in both Maringe analyses", {
  expect_equal(spec_m$pipeline$grace, 180)
  expect_equal(spec_d$pipeline$grace, 180)
})

test_that("survival difference: contract columns, plausible value, agrees with weighted KM", {
  res <- run_pipeline(base_test, spec_d$pipeline)
  expect_named(res, contract)
  expect_gt(res$estimate, 0) # surgery arm survives longer in the Maringe data
  expect_lt(res$estimate, 0.25)
  expect_gt(res$se, 0)
  expect_equal(res$ci_hi - res$estimate, qnorm(0.975) * res$se)
  long <- ccw_grace_vs_none(base_test, spec_d$pipeline)
  km <- survival::survfit(survival::Surv(start, stop, ev) ~ arm, data = long, weights = w)
  s365 <- summary(km, times = spec_d$pipeline$surv_time, extend = TRUE)$surv
  expect_equal(res$estimate, s365[2] - s365[1], tolerance = 0.03)
})

test_that("survival difference: arm-wise Hajek proportions and patient-clustered SE", {
  n <- 100
  # Two clones per patient, unit weights, nobody censored: proportions alive and the SE are known
  alive1 <- rep(c(1L, 0L), c(80, 20))
  alive0 <- rep(c(1L, 0L), c(30, 70))
  mk <- function(arm, alive) data.table(
    id = as.character(seq_len(n)), arm = arm, start = 0,
    stop = ifelse(alive == 1L, 365.24, 100), ev = 1L - alive, w = 1)
  long <- rbind(mk(1L, alive1), mk(0L, alive0))
  sp <- list(design = "grace_vs_none", surv_time = 365, conf_level = 0.95)
  r <- fit_survival_difference(long, sp)
  expect_equal(r$estimate, 0.8 - 0.3)
  infl <- (alive1 - 0.8) / n - (alive0 - 0.3) / n
  expect_equal(r$se, sqrt(sum(infl^2)))
  # Weights enter as ratios: scaling all weights leaves the estimate unchanged
  long2 <- copy(long)[, w := 5]
  expect_equal(fit_survival_difference(long2, sp)$estimate, r$estimate)
})

test_that("survival difference: artificially censored clones are excluded, deaths are not", {
  long <- data.table(id = c("a", "b", "c", "a", "b", "c"), arm = c(1L, 1L, 1L, 0L, 0L, 0L),
                     start = 0, stop = c(365.24, 100, 90, 365.24, 365.24, 365.24),
                     ev = c(0L, 1L, 0L, 0L, 0L, 0L), w = c(1, 1, 7, 1, 1, 1))
  sp <- list(design = "grace_vs_none", surv_time = 365, conf_level = 0.95)
  r <- fit_survival_difference(long, sp)
  # Arm 1: c is censored at day 90 (unknown), a alive, b dead -> 1/2; arm 0 all alive -> 1
  expect_equal(r$estimate, 0.5 - 1)
})

test_that("endpoint_info gives the direction that favours arm 1", {
  expect_true(endpoint_info(spec_m$pipeline)$favours_low)
  expect_false(endpoint_info(spec_d$pipeline)$favours_low)
  expect_true(endpoint_info(list())$favours_low)
})
