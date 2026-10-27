test_that("ANOVA components remove the inner noise from the spread of the outer estimates", {
  set.seed(11)
  B <- 4000; R <- 8; sigma_w <- 0.3; sd_true <- 0.2
  truth_b <- rnorm(B, 0, sd_true)
  # Outer-level estimate = true outer value + mean of R inner runs
  inner <- matrix(rnorm(B * R, rep(truth_b, R), sigma_w), B, R)
  x <- rowMeans(inner)
  w2 <- apply(inner, 1, var)
  cp <- nested_variance_components(x, w2, R)
  expect_equal(cp$var_obs, sd_true^2 + sigma_w^2 / R, tolerance = 0.05)
  expect_equal(cp$sigma_w2, sigma_w^2, tolerance = 0.05)
  expect_equal(cp$V, sd_true^2, tolerance = 0.1)
  expect_equal(cp$k, sigma_w^2 / sd_true^2, tolerance = 0.1)
  expect_equal(cp$shrink, sqrt(cp$V / cp$var_obs))
  xs <- shrink_deviations(x, cp$shrink)
  expect_equal(mean(xs), mean(x))
  expect_equal(var(xs), cp$V, tolerance = 1e-8)
  # Percentiles after shrinkage are closer to those of the noise-free estimator
  q <- function(v) quantile(v - mean(v), c(0.025, 0.975), names = FALSE)
  expect_lt(max(abs(q(xs) - q(truth_b))), max(abs(q(x) - q(truth_b))))
  expect_equal(shrink_deviations(x, NA_real_), x)
})

test_that("optimal inner number matches the numerical optimum of variance x cost", {
  # Var(V-hat) x cost per outer replicate for B outer and R inner runs (units V = 1, c_i = 1)
  f <- function(R, k, r) (2 * (1 + k / R)^2 + 2 * k^2 / (R^2 * (R - 1))) * (r + R)
  for (k in c(0.5, 1, 2)) for (r in c(20, 50, 100)) {
    num <- optimize(function(R) f(R, k, r), c(1.5, 200))$minimum
    expect_equal(optimal_inner(k, r, 1), num, tolerance = 0.05)
    expect_equal(optimal_inner(k, r, 1), sqrt(2 * k * r), tolerance = 0.25)
  }
  # The paper's example: k about 1 and c_o / c_i = 20 to 100 give R* of about 6 to 14
  expect_equal(round(c(optimal_inner(1, 20, 1), optimal_inner(1, 100, 1))), c(7, 15))
  # Costs enter as a ratio
  expect_equal(optimal_inner(1, 50, 1), optimal_inner(1, 500, 10))
})

test_that("outer replicates for a relative precision and percentile scaling of the allocation", {
  expect_equal(outer_for_precision(0.02), 1251L)
  ref <- list(R = 10, B = 1250, budget = 1250 * (50 + 10))
  a <- scale_percentile_allocation(ref$budget * 8, 50, 1, ref)
  expect_equal(a$R, 20) # 8^(1/3) = 2
  expect_lte(a$cost, ref$budget * 8)
  expect_gt(a$B, ref$B) # more outer replicates too, B grows like C^(2/3)
})

test_that("recommend_allocation reads the pilot and reports costs and the current design", {
  set.seed(12)
  B <- 40; R <- 20
  p <- data.table(
    rep = seq_len(B), estimate = rnorm(B, -0.48, 0.28), n_inner_ok = R, n_eval = 2,
    seconds = 30, seconds_cal = 20, mean_inner_mcse = 0.28 / sqrt(R)
  )
  p[, est_adj := estimate + rnorm(B, 0, 0.1)]
  pc <- list(precision_rel_se = 0.02, min_outer_percentile = 1000)
  a <- recommend_allocation(p, pc, list(outer = 500, inner = 250), c_study = 0.25)
  expect_equal(a$c_o_sec, 10)
  expect_equal(a$c_i_sec, 20 / R)
  expect_equal(a$B_precision, 1251L)
  expect_equal(a$B_rec, 1251L)
  expect_gte(a$R_rec, 2L)
  expect_equal(a$R_rec, a$R_variance)
  # Enough inner runs for the limits of the outer CIs override the variance optimum
  pc$min_inner_limits <- 50
  a2 <- recommend_allocation(p, pc, list(outer = 500, inner = 250), c_study = 0.25)
  expect_equal(a2$R_rec, max(a$R_variance, 50L))
  expect_equal(a2$R_variance, a$R_variance)
  expect_gt(a2$cost_rec_cpu_h, a$cost_rec_cpu_h)
  expect_equal(a$R_star, optimal_inner(a$k, a$c_o_sec, a$c_i_sec))
  expect_gt(a$cost_current_cpu_h, 0)
  expect_lt(a$relvar_rec_same_cost, a$relvar_current)
})

test_that("pilot batch keeps the inner estimates and the ANOVA is stable across sub-samples", {
  cfg_p <- cfg_test
  cfg_p$pilot <- list(outer = 3, inner = 6, batch_size = 3)
  cfg_p$nested$truth_n_target <- 600
  cfg_p$nested$truth_n_rep <- 1
  dgm <- fit_dgm(base_test, spec_m$dgm, spec_m$censoring)
  cal <- list(theta_cal = -0.7, slope = 0.55)
  r <- run_pilot_batch(base_test, dgm, cal, make_crn_seeds(5, 3), 1, spec_m, cfg_p)
  ok <- r[!is.na(est_adj)]
  expect_gt(nrow(ok), 0)
  expect_true(all(lengths(ok$inner_est) == 6))
  expect_equal(ok$mean_inner, vapply(ok$inner_est, mean, numeric(1)))
  an <- anova_by_inner_runs(r, c(2, 4, 6))
  expect_equal(an$r, c(2, 4, 6))
  expect_true(all(is.finite(an$var_obs)))
  # With all inner runs the check reproduces the summary of the nested function
  full <- nested_variance_components(ok$est_adj, ok$mean_inner_mcse^2 * ok$n_inner_ok, 6)
  expect_equal(an[r == 6, var_obs], full$var_obs, tolerance = 1e-8)
})

test_that("limits from few inner runs are too narrow", {
  set.seed(13)
  B <- 60; R <- 50
  inner <- lapply(seq_len(B), function(i) rnorm(R, 0, 0.3))
  p <- data.table(est_adj = rnorm(B), inner_est = inner)
  w <- limit_width_by_inner_runs(p, c(5, 11, 25, 50))
  expect_equal(w$r, c(5, 11, 25, 50))
  expect_true(all(diff(w$width[1:3]) > 0))
  expect_equal(w$width_rel[4], 1)
  expect_lt(w$width_rel[2], 0.9)
  # The range of 11 draws covers a new draw with probability 10/12
  expect_equal(w$width[4], 2 * 1.96 * 0.3, tolerance = 0.1)
})
