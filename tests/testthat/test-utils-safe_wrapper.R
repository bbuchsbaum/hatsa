library(testthat)

 test_that("safe_wrapper returns default and warns", {
   expect_warning(res <- safe_wrapper(stop("boom"), "Failed with: %s", default = 42), "boom")
   expect_equal(res, 42)
 })
