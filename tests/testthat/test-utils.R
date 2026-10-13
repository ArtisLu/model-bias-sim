test_that("safely_run returns a valid row for an error and a warning", {
  err <- safely_run(function() stop("boom"))
  expect_equal(err$status, "error")
  expect_true(is.na(err$estimate))
  wrn <- safely_run(function() {
    warning("careful")
    empty_result("ok", NA_character_)
  })
  expect_equal(wrn$status, "warning")
  expect_match(wrn$message, "careful")
})

test_that("same seed gives identical draws; seeds are reproducible", {
  expect_identical(make_crn_seeds(1, 5), make_crn_seeds(1, 5))
  expect_identical(draw_replicate_rng(42, 10), draw_replicate_rng(42, 10))
})

test_that("observed_workers counts sustained parallelism and ignores recording lag", {
  mk <- function(workers, n = 60, dur = 14, lag = 1) {
    start <- rep(seq(0, by = dur, length.out = n / workers), each = workers)
    data.frame(name = sprintf("sim_results_maringe_%016x", seq_len(n)), type = "branch",
               time = as.POSIXct(start + dur + lag, origin = "1970-01-01"), seconds = dur)
  }
  expect_equal(observed_workers(mk(1)), 1L)
  expect_equal(observed_workers(mk(3)), 3L)
  expect_equal(observed_workers(mk(8, n = 160)), 8L)
})
