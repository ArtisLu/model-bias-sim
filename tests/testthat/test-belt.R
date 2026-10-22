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

test_that("belt grid drops values that coincide with the main scenarios", {
  expect_equal(belt_theta_grid(list(-1, -0.7, 0, 0.02, 0.5), c(-0.7, 0)), c(-1, 0.5))
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
  ex <- extra_cis(data.table::data.table(estimate = -0.5, se = 0.3), sim, tt, belt)
  expect_true(ex$ci_lo_boot_t < -0.5 && -0.5 < ex$ci_hi_boot_t)
  expect_true(ex$ci_lo_belt < -0.5 && -0.5 < ex$ci_hi_belt)
})
