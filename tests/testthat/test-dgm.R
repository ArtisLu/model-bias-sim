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

test_that("simulated data carry their own loss-to-follow-up indicator", {
  spec_c <- list(model = "weibull", covariates = spec_s$censoring$covariates)
  dgm <- fit_dgm(base_syn, spec_s$dgm, spec_c)
  rng <- draw_replicate_rng(12, nrow(base_syn))
  sim <- simulate_outcomes(base_syn, dgm, 0, rng$u_event, rng$u_cens)
  expect_equal(sim$cens_event, as.integer(sim$status == 0L & sim$time < sim$admin_time))
  expect_true(all(simulate_outcomes(base_test, fit_dgm(base_test, spec_m$dgm), 0,
                                    draw_replicate_rng(1, nrow(base_test))$u_event)$cens_event == 0L))
})

test_that("generalized gamma baseline: inverse, hazard and Weibull special case", {
  for (Q in c(-1, -0.4, 0.3, 1, 2)) for (sg in c(0.7, 1.5)) {
    t <- c(0.2, 0.9, 2, 6)
    H <- gg_cumhaz0(t, sg, Q)
    expect_equal(gg_time0(H, sg, Q), t, tolerance = 1e-8)
    h_num <- (gg_cumhaz0(t * (1 + 1e-6), sg, Q) - gg_cumhaz0(t * (1 - 1e-6), sg, Q)) / (2e-6 * t)
    expect_equal(exp(gg_loghaz0(t, sg, Q)), h_num, tolerance = 1e-5)
  }
  # Q = 1 is the Weibull with shape 1 / sigma; Q -> 0 is the log-normal
  expect_equal(gg_cumhaz0(c(0.3, 2), 0.5, 1), c(0.3, 2)^2, tolerance = 1e-10)
  expect_equal(gg_cumhaz0(2, 1, 1e-4), gg_cumhaz0(2, 1, 0), tolerance = 1e-3)
  expect_equal(gg_time0(0, 1, 0.5), 0)
  expect_equal(gg_time0(Inf, 1, 0.5), Inf)
})

test_that("generalized gamma PH model nests the Weibull and simulates by inversion", {
  sp <- spec_m$dgm
  sp$family <- "gengamma"
  wb <- fit_dgm(base_test, spec_m$dgm, spec_m$censoring)
  gg <- fit_dgm(base_test, sp, spec_m$censoring)
  expect_equal(gg$family, "gengamma")
  # Started at the Weibull fit, the likelihood cannot be lower
  expect_gte(gg$loglik, gg$loglik_weibull - 1e-6)
  expect_equal(gg$loglik_weibull, wb$loglik, tolerance = 1e-6)
  expect_equal(gg$theta_hat, wb$theta_hat, tolerance = 0.1)
  # Q = 1 reproduces the Weibull simulation exactly
  g1 <- gg
  g1$Q <- 1
  g1$sigma <- 1 / wb$rho
  g1$alpha <- wb$alpha + wb$rho * log(g1$tau)
  g1$gamma <- wb$gamma
  u <- draw_replicate_rng(3, nrow(base_test))$u_event
  expect_equal(simulate_outcomes(base_test, g1, -0.5, u)$time,
               simulate_outcomes(base_test, wb, -0.5, u)$time, tolerance = 1e-8)
  # CRN: same u gives times that are monotone in theta, untreated patients unchanged
  t_lo <- simulate_outcomes(base_test, gg, -1, u)$time
  t_hi <- simulate_outcomes(base_test, gg, 0, u)$time
  expect_true(all(t_lo >= t_hi - 1e-9))
  untreated <- !is.finite(base_test$s)
  expect_equal(t_lo[untreated], t_hi[untreated])
})

test_that("generalized gamma refit recovers the treatment effect of a large simulated dataset", {
  sp <- spec_m$dgm
  sp$family <- "gengamma"
  gg <- fit_dgm(base_test, sp)
  big <- stack_copies(base_test, 40)
  sim <- simulate_outcomes(big, gg, -0.6, draw_replicate_rng(7, nrow(big))$u_event)
  re <- fit_dgm(sim, sp)
  expect_equal(re$theta_hat, -0.6, tolerance = 0.1)
  # marginal survival of the fitted model follows the Kaplan-Meier curve of the simulated data
  tm <- c(100, 200, 300)
  km <- summary(survival::survfit(survival::Surv(time, status) ~ 1, data = sim), times = tm)$surv
  expect_equal(dgm_marginal_surv(re, sim, tm), km, tolerance = 0.03)
})
