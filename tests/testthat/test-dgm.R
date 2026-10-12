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

test_that("fixed censoring mode keeps each patient's observed censoring time", {
  dgm <- fit_dgm(base_syn, spec_s$dgm)
  sim <- simulate_outcomes(base_syn, dgm, 0, draw_replicate_rng(1, nrow(base_syn))$u_event)
  expect_true(all(sim$time <= base_syn$cens_time + 1e-9))
})

test_that("censoring model: simulated censoring respects the horizon and recovers its rate", {
  spec_c <- list(model = "weibull", covariates = spec_s$censoring$covariates)
  dgm <- fit_dgm(base_syn, spec_s$dgm, spec_c)
  expect_equal(dgm$cens$n_events, sum(base_syn$cens_event))
  big <- stack_copies(base_syn, 20)
  rng <- draw_replicate_rng(11, nrow(big))
  sim <- simulate_outcomes(big, dgm, dgm$theta_hat, rng$u_event, rng$u_cens)
  expect_true(all(sim$cens_time <= sim$admin_time + 1e-9))
  lost <- mean(sim$status == 0L & sim$time < sim$admin_time)
  expect_equal(lost, mean(base_syn$cens_event), tolerance = 0.3)
})

test_that("adding u_cens leaves the event-time draws unchanged (appended, not reordered)", {
  r <- draw_replicate_rng(5, 10)
  expect_identical(r$u_event, withr::with_seed(5, runif(10), .rng_kind = "L'Ecuyer-CMRG"))
  expect_length(r$u_cens, 10)
})
