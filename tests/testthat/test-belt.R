test_that("belt inversion recovers a known non-pivotal T and covers at an interpolated effect", {
  set.seed(3)
  # T ~ N(m(psi), 1) with a mean that drifts with psi: bootstrap-t at one psi would be off elsewhere
  m <- function(psi) -0.4 * psi
  grid <- c(-1.5, -1, -0.5, 0.5, 1, 1.5)
  sims <- data.table::rbindlist(lapply(grid, function(p) {
    data.table::data.table(theta = p, rep = 1:4000, se = 0.3, estimate = p + 0.3 * (m(p) + rnorm(4000)))
  }))
  belt <- fit_belt(sims, data.table::data.table(theta = grid, truth = grid), list(degree = 1))
  expect_equal(unname(belt$curve$lo), c(stats::qnorm(0.025), -0.4), tolerance = 0.08)
  # Check at psi = 0, which is not on the grid
  x <- 0.3 * (m(0) + rnorm(4000))
  ci <- belt_invert(x, rep(0.3, 4000), belt$curve)
  expect_true(all(ci$lo < ci$hi))
  expect_equal(mean(ci$lo <= 0 & 0 <= ci$hi), 0.95, tolerance = 0.012)
})

test_that("belt nodes span the reachable psi range, skip the main scenarios and map to theta", {
  expect_equal(cheb_lobatto(c(-1, 1), 3), c(-1, 0, 1))
  nodes <- cheb_lobatto(c(-2, 1), 5)
  expect_equal(range(nodes), c(-2, 1))
  expect_true(all(diff(nodes) > 0))
  # psi(theta) levels off for strong effects, as the CCW limit does; anchors at -4 and 4
  f <- function(th) 1.5 * tanh(th / 2) - 0.04
  th <- rep(c(-4, -3, -2, -1, -0.7, 0, 1, 2, 3, 4), each = 2)
  pts <- data.table::data.table(theta = th, seed_id = 1:2, n_truth = 10L, estimate = f(th))
  tt <- data.table::data.table(theta = c(-0.7, 0), truth = c(f(-0.7), 0),
                               psi_large_sample = f(c(-0.7, 0)))
  cb <- list(n_nodes = 5, min_gap = 0.05)
  out <- belt_theta_plan(c(-1.96, 1.5), tt, pts, cb)
  # Both requested ends lie beyond psi(-4) and psi(4): the grid runs from anchor to anchor
  expect_equal(range(out), c(-4, 4))
  nodes <- cheb_lobatto(f(c(-4, 4)), 5)
  # The middle node coincides with psi(0) = -0.04 and is dropped
  expect_equal(f(out), nodes[-3], tolerance = 0.01)
  # A node at a main-scenario psi is dropped
  out2 <- belt_theta_plan(c(f(-0.7), 1.5), tt, pts, list(n_nodes = 3, min_gap = 0.05))
  expect_length(out2, 2L)
})

test_that("belt range covers the observed CI and the replicate endpoints", {
  sims <- data.table::data.table(ci_lo = c(-2, -1, 0), ci_hi = c(0, 1, 2))
  r <- belt_psi_range(list(ci_lo = -3, ci_hi = 0.5), sims, tail = 0, expand = 0)
  expect_equal(r, c(-3, 2))
  expect_equal(belt_psi_range(list(ci_lo = -3, ci_hi = 0.5), sims, 0, 0.1), c(-3.5, 2.5))
})

test_that("quantile MCSE matches the normal-theory value", {
  set.seed(5)
  m <- replicate(200, quantile_mcse(rnorm(2000), 0.975))
  theory <- sqrt(0.975 * 0.025 / 2000) / stats::dnorm(stats::qnorm(0.975))
  expect_equal(mean(m), theory, tolerance = 0.15)
})

test_that("top-up plan asks for more replicates only where the quantiles are imprecise", {
  set.seed(6)
  sims <- data.table::data.table(theta = rep(c(-1, 1), c(1000, 1000)), rep = c(1:1000, 1:1000),
                                 se = 1, estimate = c(rnorm(1000), 3 * rnorm(1000)))
  truths <- data.table::data.table(theta = c(-1, 1), truth = 0)
  cb <- list(B_start = 1000, B_max = 5000, target_mcse = 0.15)
  pl <- belt_topup_plan(sims, truths, cb, batch_size = 250)
  # MCSE ~0.09 at SD 1 (no top-up) and ~0.26 at SD 3 (about 3,000 replicates)
  expect_equal(pl[theta == -1, n_batches], 4L)
  expect_gt(pl[theta == 1, n_batches], 8L)
  # An unneeded batch returns no rows without running the pipeline
  out <- run_belt_topup(NULL, NULL, -1, 1:5000, 5L, pl, NULL, 250)
  expect_equal(nrow(out), 0L)
})

test_that("belt coverage rows fit the CI comparison format", {
  set.seed(4)
  mk <- function(th, psi) data.table::data.table(dataset = "m", theta = th, rep = 1:400, se = 0.3,
                                                 estimate = psi + 0.3 * rnorm(400), status = "ok")
  sim <- rbind(mk(-0.7, -0.5), mk(0, 0))
  extra <- rbind(mk(-1.5, -1), mk(1, 0.6))
  tt <- data.table::data.table(dataset = "m", theta = c(-0.7, 0), scenario = c("observed", "null"),
                               truth = c(-0.5, 0), psi_large_sample = c(-0.5, 0))
  belt <- fit_belt(rbind(sim, extra), data.table::data.table(theta = c(-0.7, 0, -1.5, 1),
                                                             truth = c(-0.5, 0, -1, 0.6)),
                   list(degree = 1))
  cm <- compare_ci_methods(sim, tt, belt)
  expect_equal(nrow(cm[method == "Confidence belt"]), 2L)
  expect_true(all(cm[method == "Confidence belt", coverage] > 0.9))
})

test_that("belt truths use the CCW limit at every grid point, including no effect", {
  tt <- data.table::data.table(theta = c(-0.7, 0), truth = c(-0.5, 0),
                               scenario = c("observed", "null"), psi_large_sample = c(-0.5, -0.04))
  bt <- data.table::data.table(theta = c(-2, -0.2, 1.5), truth = c(-1.4, -0.15, 1), n_rep = 12L)
  out <- belt_truths(tt, bt)
  expect_named(out, c("theta", "truth"))
  expect_equal(out$truth, c(-0.5, -0.04, -1.4, -0.15, 1))
})

test_that("belt coverage of the null scenario is checked at psi(0)", {
  set.seed(7)
  mk <- function(th, psi) data.table::data.table(dataset = "m", theta = th, rep = 1:2000, se = 0.25,
                                                 estimate = psi + 0.25 * rnorm(2000), status = "ok")
  # CCW limit at no effect is -0.1, not 0: CIs built around psi(0) should cover psi(0)
  sim <- rbind(mk(-0.7, -0.5), mk(0, -0.1))
  extra <- rbind(mk(-1.5, -1), mk(1, 0.6))
  tt <- data.table::data.table(dataset = "m", theta = c(-0.7, 0), scenario = c("observed", "null"),
                               truth = c(-0.5, 0), psi_large_sample = c(-0.5, -0.1))
  bt <- data.table::data.table(theta = c(-1.5, 1), truth = c(-1, 0.6))
  belt <- fit_belt(rbind(sim, extra), belt_truths(tt, bt), list(degree = 1))
  cv <- belt_coverage(sim, tt, belt)
  expect_equal(cv[scenario == "null", coverage], 0.95, tolerance = 0.015)
  expect_equal(cv[scenario == "null", miss_below], 0.025, tolerance = 0.01)
})
