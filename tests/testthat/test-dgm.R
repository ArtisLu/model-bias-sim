test_that("DGM recovers its own parameters on a large simulated dataset", {
  dgm <- fit_dgm(base_test, spec_m$dgm)
  big <- stack_copies(base_test, 50)
  u <- draw_replicate_rng(7, nrow(big))$u_event
  sim <- simulate_outcomes(big, dgm, theta = -0.8, u = u)
  refit <- fit_dgm(sim, spec_m$dgm)
  expect_equal(refit$theta_hat, -0.8, tolerance = 0.1)
  expect_equal(refit$rho, dgm$rho, tolerance = 0.05)
})

test_that("CRN: same u gives event times that are monotone in theta", {
  dgm <- fit_dgm(base_test, spec_m$dgm)
  u <- draw_replicate_rng(3, nrow(base_test))$u_event
  d <- copy(base_test)[, cens_time := Inf]
  t_lo <- simulate_outcomes(d, dgm, -1, u)$time
  t_hi <- simulate_outcomes(d, dgm, 0, u)$time
  expect_true(all(t_lo >= t_hi - 1e-9))
  untreated <- !is.finite(d$s)
  expect_equal(t_lo[untreated], t_hi[untreated])
})

test_that("simulation keeps each patient's censoring time", {
  dgm <- fit_dgm(base_syn, spec_s$dgm)
  sim <- simulate_outcomes(base_syn, dgm, 0, draw_replicate_rng(1, nrow(base_syn))$u_event)
  expect_true(all(sim$time <= base_syn$cens_time + 1e-9))
})
