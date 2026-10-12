contract <- c("estimate", "se", "ci_lo", "ci_hi", "p_value", "n_events", "status", "message")

test_that("both pipelines return the contract columns on the real data", {
  for (x in list(list(base_test, spec_m), list(base_syn, spec_s))) {
    res <- run_pipeline(x[[1]], x[[2]]$pipeline)
    expect_named(res, contract)
    expect_true(is.finite(res$estimate))
  }
})

test_that("Maringe cloning keeps two clones per patient and censors deviators", {
  cl <- build_clones(base_test, spec_m$pipeline$grace)
  expect_equal(nrow(cl), 2 * nrow(base_test))
  op <- cl[arm == 0L & is.finite(trt_time)]
  expect_true(all(op$art_cens == 1L & op$end == op$trt_time))
})

test_that("tutorial design: non-starters followed past the window are censored at it", {
  long <- ccw_window_arm(base_syn, 30, spec_s$pipeline)
  ns <- base_syn[time > 30 & trt_time > 30, id]
  expect_true(all(long[id %in% ns, max(stop), by = id]$V1 == 30))
  expect_true(all(long[start < 30, w] == 1))
})

test_that("results do not depend on batch size", {
  dgm <- fit_dgm(base_test, spec_m$dgm)
  seeds <- make_crn_seeds(1, 4)
  r2 <- rbind(
    run_batch(base_test, dgm, -0.5, seeds, 1, spec_m$pipeline, 2),
    run_batch(base_test, dgm, -0.5, seeds, 2, spec_m$pipeline, 2)
  )
  r4 <- run_batch(base_test, dgm, -0.5, seeds, 1, spec_m$pipeline, 4)
  expect_equal(r2$estimate, r4$estimate)
})
