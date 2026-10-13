test_that("baseline tables have Table S2 columns and counts that add up", {
  for (x in list(list("maringe", spec_m), list("synpuf", spec_s))) {
    raw <- read_raw_data(test_path("..", "..", x[[2]]$file), x[[1]])
    tb <- make_baseline_table(raw, x[[1]], x[[2]])
    expect_named(tb, c("label", "indent", "n1", "p1", "n0", "p0", "smd"))
    tot <- tb[label == "Total"]
    expect_equal(as.integer(tot$n1) + as.integer(tot$n0), nrow(raw))
    men <- tb[label == "Men"]; women <- tb[label == "Women"]
    expect_equal(as.integer(men$n1) + as.integer(women$n1), as.integer(tot$n1))
  }
})

test_that("Maringe table reproduces the surgery split of the data", {
  raw <- read_raw_data(test_path("..", "..", spec_m$file), "maringe")
  tb <- make_baseline_table(raw, "maringe", spec_m)
  expect_equal(as.integer(tb[label == "Total", n1]), sum(!is.na(raw$timetosurgery)))
})

test_that("cross-fitted adjustment returns both analyses and sensible proportions", {
  set.seed(1)
  sim <- data.table(
    dataset = "synpuf", theta = rep(c(-0.2, 0), each = 200), rep = rep(1:200, 2),
    estimate = c(rnorm(200, -0.15, 0.1), rnorm(200, 0, 0.1)), se = 0.09, status = "ok"
  )
  tt <- data.table(dataset = "synpuf", theta = c(-0.2, 0), scenario = c("observed", "null"),
                   truth = c(-0.2, 0))
  pa <- evaluate_adjustment(sim, tt)
  expect_setequal(unique(pa$analysis), c("CCW", "CCW + QBA adjustment"))
  adj_bias <- pa[analysis != "CCW" & measure == "Bias (log HR)", value]
  expect_lt(abs(adj_bias), abs(pa[analysis == "CCW" & measure == "Bias (log HR)", value]))
  expect_true(all(pa[measure != "Bias (log HR)", value] >= 0 & pa[measure != "Bias (log HR)", value] <= 1))
})

test_that("bootstrap-t CI is asymmetric when T is skewed and recovers nominal coverage", {
  set.seed(2)
  # SE correlated with the estimate, as in the Maringe replicates
  e <- rnorm(4000, -0.5, 0.28)
  sim <- data.table(dataset = "m", theta = rep(c(-0.7, 0), each = 2000), rep = rep(1:2000, 2),
                    estimate = c(e[1:2000], e[2001:4000] + 0.5), status = "ok")
  sim[, se := 0.27 + 0.1 * (estimate - mean(estimate)), by = theta]
  tt <- data.table(dataset = "m", theta = c(-0.7, 0), scenario = c("observed", "null"),
                   truth = c(-0.5, 0))
  q <- boot_t_quantiles(sim[theta == -0.7, estimate], sim[theta == -0.7, se], -0.5)
  expect_lt(q[["q_lo"]], 0)
  expect_gt(q[["q_hi"]], 0)
  expect_false(isTRUE(all.equal(-q[["q_lo"]], q[["q_hi"]], tolerance = 0.02)))
  ci <- boot_t_ci(-0.5, 0.27, q)
  expect_lt(ci$lo, -0.5)
  expect_gt(ci$hi, -0.5)
  pb <- evaluate_boot_t(sim, tt)
  expect_equal(nrow(pb), 4L)
  expect_true(all(abs(pb[check == "cross-fitted", coverage] - 0.95) < 0.02))
})

test_that("balance at the grace period uses the clones still at risk and both weightings", {
  b <- ccw_balance(base_test, spec_m$pipeline, spec_m$dgm$covariates)
  expect_setequal(unique(b$weighting), c("Unweighted", "Weighted"))
  at_risk <- base_test[time > spec_m$pipeline$grace]
  expect_equal(b$n1[1], sum(at_risk$trt_time <= spec_m$pipeline$grace))
  expect_equal(b$n0[1], sum(at_risk$trt_time > spec_m$pipeline$grace))
  expect_true(all(is.finite(b$smd)))
})
