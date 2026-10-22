cfg_n <- cfg_test
cfg_n$nested <- list(outer = 2, inner = 6, batch_size = 2, max_iter = 1, tol = 0.05,
                     truth_n_target = 600, truth_n_rep = 1)

test_that("nested outer replicate: same outer data as validation, own calibration and CI", {
  dgm <- fit_dgm(base_test, spec_m$dgm, spec_m$censoring)
  cal <- list(theta_cal = -0.7, slope = 0.55)
  seeds <- make_crn_seeds(3, 2)
  r <- run_nested_batch(base_test, dgm, cal, seeds, 1, spec_m, cfg_n)
  expect_equal(nrow(r), 2L)
  v <- run_batch(base_test, dgm, -0.7, seeds, 1, spec_m$pipeline, 2)
  expect_equal(r$estimate, v$estimate)
  ok <- r[!is.na(est_adj)]
  expect_gt(nrow(ok), 0)
  expect_equal(ok$est_adj, ok$estimate - ok$bias_r)
  expect_equal(ok$bias_r, ok$mean_inner - ok$truth_r)
  expect_true(all(ok$ci_lo_np < ok$est_adj & ok$est_adj < ok$ci_hi_np))
  expect_true(all(ok$n_eval >= 2))
})

test_that("nested null replicate returns a rejection decision from its own limits", {
  dgm <- fit_dgm(base_test, spec_m$dgm, spec_m$censoring)
  r <- run_nested_batch(base_test, dgm, NULL, make_crn_seeds(4, 2), 1, spec_m, cfg_n)
  expect_true(all(r$reject %in% c(0, 1)))
  expect_true(all(r$reject_lower <= r$reject))
})

test_that("calibration finds the mean-estimate root under CRN", {
  dgm <- fit_dgm(base_test, spec_m$dgm, spec_m$censoring)
  seeds <- make_crn_seeds(9, 15)
  cal <- calibrate_mean(dgm, base_test, -0.477, seeds, spec_m$pipeline,
                        list(max_iter = 4, tol = 0.005, bound = 3))
  expect_gt(cal$slope, 0) # Maringe: m(theta) increases in theta
  expect_lt(abs(cal$mean_est - (-0.477)), 0.02)
  expect_equal(cal$mean_est, mean(replicate_estimates(base_test, dgm, cal$theta_cal, seeds,
                                                      spec_m$pipeline), na.rm = TRUE))
})
