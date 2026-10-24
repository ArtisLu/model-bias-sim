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

test_that("belt nodes span the range, skip the main scenarios and map to theta", {
  expect_equal(cheb_lobatto(c(-1, 1), 3), c(-1, 0, 1))
  nodes <- cheb_lobatto(c(-2, 1), 5)
  expect_equal(range(nodes), c(-2, 1))
  expect_true(all(diff(nodes) > 0))
  # psi = 0.8 * theta - 0.1 on the evaluated points
  th <- rep(c(-1, 0, -0.5), each = 2)
  pts <- data.table::data.table(theta = th, seed_id = 1:2, n_truth = 10L, estimate = 0.8 * th - 0.1)
  tt <- data.table::data.table(theta = c(-1, 0), truth = c(-0.9, 0))
  cb <- list(n_nodes = 3, min_gap = 0.05, bound = 4)
  # Nodes at psi = -2, -0.5, 1; none within min_gap of the main truths
  expect_equal(belt_theta_plan(c(-2, 1), tt, pts, cb), (c(-2, -0.5, 1) + 0.1) / 0.8)
  # Nodes at the main-scenario truths (-0.9 and 0) are dropped
  expect_equal(belt_theta_plan(c(-1.8, 0), tt, pts, cb), (-1.8 + 0.1) / 0.8)
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
                               truth = c(-0.5, 0))
  belt <- fit_belt(rbind(sim, extra), data.table::data.table(theta = c(-0.7, 0, -1.5, 1),
                                                             truth = c(-0.5, 0, -1, 0.6)),
                   list(degree = 1))
  cm <- compare_ci_methods(sim, tt, belt)
  expect_equal(nrow(cm[method == "Confidence belt"]), 2L)
  expect_true(all(cm[method == "Confidence belt", coverage] > 0.9))
})

test_that("belt truths stack the main scenarios and the nodes", {
  tt <- data.table::data.table(theta = c(-0.7, 0), truth = c(-0.5, 0), scenario = c("observed", "null"))
  bt <- data.table::data.table(theta = c(-2, -0.2, 1.5), truth = c(-1.4, -0.15, 1), n_rep = 12L)
  out <- belt_truths(tt, bt)
  expect_named(out, c("theta", "truth"))
  expect_equal(out$truth, c(-0.5, 0, -1.4, -0.15, 1))
})
