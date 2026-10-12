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
