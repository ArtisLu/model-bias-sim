test_that("stacked copies have unique ids", {
  big <- stack_copies(base_test, 3)
  expect_equal(nrow(big), 3 * nrow(base_test))
  expect_equal(uniqueN(big$id), nrow(big))
})

test_that("truth points reuse earlier evaluations and match the serial truth", {
  dgm <- fit_dgm(base_test, spec_m$dgm, spec_m$censoring)
  seeds <- make_crn_seeds(cfg_test$seed_truth, 2)
  rows <- rbindlist(lapply(1:2, function(j) {
    truth_point(dgm, base_test, -0.5, j, seeds, spec_m$pipeline, cfg_test)
  }))
  expect_equal(summarise_truth(rows)$truth,
               compute_truth(dgm, base_test, -0.5, spec_m$pipeline, cfg_test, n_rep = 2)$truth)
  # A known (theta, seed) is returned as is, without simulating
  other <- data.table(theta = 0, seed_id = 2L, estimate = 7, n_truth = 1L)
  fake <- rbind(copy(rows)[, estimate := c(98, 99)], other)
  expect_equal(truth_point(NULL, NULL, -0.5, 2L, seeds, NULL, NULL, known = fake)$estimate, 99)
  expect_equal(truth_point(NULL, NULL, 0, 2L, seeds, NULL, NULL, known = fake)$estimate, 7)
})

test_that("calibration steps converge on a smooth psi(theta)", {
  psi <- function(th) 0.1 + 0.7 * th + 0.05 * th^2
  ev <- function(th) data.table(theta = th, seed_id = 1L, estimate = psi(th), n_truth = 1L)
  cc <- list(tol = 0.005, bound = 3)
  target <- psi(-0.9)
  expect_equal(calib_start(list(theta_hat = -1.2)), c(-1.2, 0))
  expect_equal(calib_start(list(theta_hat = 0.1)), c(0.1, -0.4))
  rows <- rbind(ev(-1.2), ev(0))
  for (i in 1:2) rows <- rbind(rows, ev(calib_next(rows, target, cc)))
  best <- summarise_truth(rows)[which.min(abs(truth - target))]
  expect_lt(abs(best$truth - target), cc$tol)
  # Once within tolerance the best evaluated point is returned, so the next stage is free
  expect_equal(calib_next(rows, target, cc), best$theta)
})

test_that("psi_to_theta inverts the monotone psi curve and stops at the range ends", {
  pts <- data.table(theta = c(-4, -1, 0, 1, 4), psi = c(-1.4, -0.8, 0, 0.8, 1.6))
  f <- psi_curve(pts)
  expect_equal(f(c(-1, 0, 1)), c(-0.8, 0, 0.8))
  th <- psi_to_theta(c(-2, -0.4, 0.8, 3), f, c(-4, 4))
  expect_equal(th[c(1, 4)], c(-4, 4))
  expect_equal(f(th[2:3]), c(-0.4, 0.8), tolerance = 1e-6)
  expect_true(all(diff(f(seq(-4, 4, 0.1))) >= 0))
})
