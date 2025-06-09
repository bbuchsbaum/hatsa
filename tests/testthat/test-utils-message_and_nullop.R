context("message_stage and null-coalescing operator")

test_that("message_stage and %||% basic functionality", {
  expect_message(message_stage("hello", verbose = TRUE, interactive_only = FALSE), "hello")
  expect_silent(message_stage("quiet", verbose = FALSE, interactive_only = FALSE))
  expect_silent(message_stage("skip", verbose = TRUE, interactive_only = TRUE))

  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
  expect_null(NULL %||% NULL)
})
