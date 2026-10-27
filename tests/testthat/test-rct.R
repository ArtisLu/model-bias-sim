cfg_r <- cfg_test
cfg_r$rct <- list(control_censor_day = 181)

test_that("surgery arm keeps own dates within grace and draws donor dates otherwise", {
  s <- c(10, 50, 179, 181, 300, Inf)
  pool <- c(5, 20, 40)
  out <- rct_surgery_times(s, 180, pool, 1)
  expect_equal(out[1:3], s[1:3])
  expect_true(all(out[4:6] %in% pool))
  expect_identical(out, rct_surgery_times(s, 180, pool, 1))
})

test_that("hypothetical trial: null truth is zero and a protective effect gives a negative log HR", {
  dgm <- fit_dgm(base_test, spec_m$dgm, spec_m$censoring)
  pl <- spec_m$pipeline
  n0 <- compute_truth_rct(dgm, base_test, 0, pl, cfg_r, n_rep = 4, n_target = 20000, seed = 3)
  expect_lt(abs(n0$truth), 4 * n0$truth_mcse + 0.02)
  t1 <- compute_truth_rct(dgm, base_test, -0.7, pl, cfg_r, n_rep = 2, n_target = 20000, seed = 3)
  expect_lt(t1$truth, -0.05)
  # No weights and no confounding: stronger effect, stronger log HR (CRN across theta)
  t2 <- compute_truth_rct(dgm, base_test, -1.4, pl, cfg_r, n_rep = 2, n_target = 20000, seed = 3)
  expect_lt(t2$truth, t1$truth)
  # Follow-up of the arm without surgery to the end changes the value
  cfg_f <- cfg_r
  cfg_f$rct <- list(control_censor_day = NULL)
  t3 <- compute_truth_rct(dgm, base_test, -0.7, pl, cfg_f, n_rep = 2, n_target = 20000, seed = 3)
  expect_false(isTRUE(all.equal(t1$truth, t3$truth, tolerance = 1e-6)))
})

test_that("hypothetical trial is only defined for the hazard-ratio endpoint", {
  dgm <- fit_dgm(base_test, spec_d$dgm, spec_d$censoring)
  expect_error(compute_truth_rct(dgm, base_test, -0.5, spec_d$pipeline, cfg_r, 1, 1000, 1))
})

test_that("nested update replaces the truth and keeps the inner runs", {
  dgm <- fit_dgm(base_test, spec_m$dgm, spec_m$censoring)
  cfg_n <- cfg_r
  cfg_n$nested <- list(outer = 2, inner = 6, batch_size = 2, max_iter = 1, tol = 0.05,
                       truth_n_target = 600, truth_n_rep = 1)
  cal <- list(theta_cal = -0.7, slope = 0.55)
  seeds <- make_crn_seeds(3, 2)
  r <- run_nested_batch(base_test, dgm, cal, seeds, 1, spec_m, cfg_n)
  u <- rct_nested_update(r, base_test, dgm, cal, seeds, spec_m, cfg_n)
  ok <- u[!is.na(est_adj)]
  expect_gt(nrow(ok), 0)
  expect_equal(ok$mean_inner, r[!is.na(est_adj)]$mean_inner)
  expect_equal(ok$bias_r, ok$mean_inner - ok$truth_r)
  expect_equal(ok$est_adj, ok$estimate - ok$bias_r)
  expect_false(isTRUE(all.equal(ok$truth_r, r[!is.na(est_adj)]$truth_r)))
  expect_true(all(ok$ci_lo_np < ok$est_adj & ok$est_adj < ok$ci_hi_np))
})
