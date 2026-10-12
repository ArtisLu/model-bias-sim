test_that("baseline table has Table S2 columns and counts that add up", {
  raw <- haven::read_sas(test_path("..", "..", spec_s$file))
  tb <- make_baseline_table(raw, spec_s)
  expect_named(tb, c("label", "indent", "n1", "p1", "n0", "p0", "smd"))
  tot <- tb[label == "Total"]
  expect_equal(as.integer(tot$n1) + as.integer(tot$n0), nrow(raw))
  men <- tb[label == "Men"]; women <- tb[label == "Women"]
  expect_equal(as.integer(men$n1) + as.integer(women$n1), as.integer(tot$n1))
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
