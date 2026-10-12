test_that("stacked copies have unique ids", {
  big <- stack_copies(base_test, 3)
  expect_equal(nrow(big), 3 * nrow(base_test))
  expect_equal(uniqueN(big$id), nrow(big))
})
