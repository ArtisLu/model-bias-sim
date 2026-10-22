mk_sims <- function(err_obs, err_null, truth = -0.475, theta = -0.7, se = 0.27) {
  n <- length(err_obs)
  data.table(dataset = "m", theta = rep(c(theta, 0), each = n), rep = rep(seq_len(n), 2),
             estimate = c(truth + err_obs, err_null), se = se, status = "ok")
}
tt_m <- data.table(dataset = "m", theta = c(-0.7, 0), scenario = c("observed", "null"),
                   truth = c(-0.475, 0), truth_mcse = c(0.002, 0), psi_large_sample = c(-0.475, -0.03))
obs_m <- data.table(estimate = -0.477, se = 0.265, ci_lo = -0.996, ci_hi = 0.042, p_value = 0.072)

test_that("simple adjustment follows paper 2.7.1-2.7.3", {
  set.seed(1)
  sim <- mk_sims(rnorm(2000, -0.015, 0.28), rnorm(2000, -0.04, 0.28))
  adj <- adjust_simple(obs_m, sim, tt_m, "m")
  o <- sim[theta == -0.7, estimate]
  expect_equal(adj$bias, mean(o) - (-0.475))
  # psi(theta*) = observed - b(theta*)
  expect_equal(adj$est_adj, obs_m$estimate - adj$bias)
  q <- stats::quantile(o - mean(o), c(0.025, 0.975), names = FALSE, type = 8)
  expect_equal(c(adj$ci_lo_np, adj$ci_hi_np), adj$est_adj + q)
  expect_equal(c(adj$ci_lo_np_refl, adj$ci_hi_np_refl), adj$est_adj - rev(q))
  expect_gt(adj$bias_mcse, adj$mean_est_mcse)
  expect_true(adj$p_np >= 0 && adj$p_np <= 1)
  expect_equal(adj$reject_np, adj$p_np < 0.05 || abs(adj$p_np - 0.05) < 1e-3)
})

test_that("non-parametric CI and test hold their level in independent replicates", {
  set.seed(2)
  # Skewed errors: the reflected form is exact, the paper's form is not
  draw <- function(n) 0.28 * (stats::rgamma(n, 4, 2) - 2) / sqrt(1)
  adj_sk <- adjust_simple(obs_m, mk_sims(draw(20000), draw(20000) - 0.04), tt_m, "m")
  val_sk <- validate_adjustment(adj_sk, mk_sims(draw(20000), draw(20000) - 0.04), tt_m)
  cov <- function(v, m) v[method == m & measure == "95% CI coverage", value]
  expect_equal(cov(val_sk, "Non-parametric (reflected)"), 0.95, tolerance = 0.006)
  expect_true(abs(cov(val_sk, "Non-parametric") - 0.95) > 0.01)
  # Symmetric errors: both forms are nominal; null shifted and skewed but test still at 5%
  adj <- adjust_simple(obs_m, mk_sims(rnorm(20000, -0.015, 0.28), draw(20000) - 0.04), tt_m, "m")
  val <- validate_adjustment(adj, mk_sims(rnorm(20000, -0.015, 0.28), draw(20000) - 0.04), tt_m)
  expect_equal(cov(val, "Non-parametric"), 0.95, tolerance = 0.006)
  t1 <- val[method == "Non-parametric" & measure == "Type I error", value]
  expect_equal(t1, 0.05, tolerance = 0.005)
  expect_equal(val[method == "Non-parametric" & measure == "Rejections favouring arm 1", value],
               0.025, tolerance = 0.004)
  expect_equal(val[method == "Non-parametric" & measure == "Bias (log HR)", value], 0, tolerance = 0.01)
})

test_that("performance rows: misses add up to non-coverage", {
  set.seed(3)
  est <- rnorm(500)
  r <- perf_rows("A", "x", est, est - 1.5, est + 1.5, 0, rej_null = est > 1.5, low_null = est < -1.5)
  cv <- r[measure == "95% CI coverage", value]
  expect_equal(r[measure == "CI entirely below truth", value] + r[measure == "CI entirely above truth", value],
               1 - cv)
  expect_setequal(r$measure, c("Mean estimate", "Bias (log HR)", "Empirical SE", "RMSE",
                               "95% CI coverage", "CI entirely below truth", "CI entirely above truth",
                               "Type I error", "Rejections favouring arm 1"))
})
