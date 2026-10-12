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
