library(testthat)

# Test hatsa_suggest with typical dataset

test_that("hatsa_suggest returns sensible parameters", {
  set.seed(123)
  data <- lapply(1:3, function(i) matrix(rnorm(60 * 200), 60, 200))
  params <- hatsa_suggest(data, verbose = FALSE)
  expect_true(is.list(params))
  expect_true(params$components <= 30)
  expect_true(params$n_anchors >= params$components)
  expect_true(params$n_refine >= 2)
})

# Test hatsa_validate_params error and warning conditions

test_that("hatsa_validate_params checks anchors and components", {
  data <- lapply(1:2, function(i) matrix(rnorm(10 * 20), 10, 20))
  expect_error(hatsa_validate_params(data, 1:5, spectral_rank_k = 11))
  expect_warning(hatsa_validate_params(data, 1:15, spectral_rank_k = 5), "unusually high")
})

# Test select_task_method heuristic

test_that("select_task_method chooses method based on data", {
  data <- list(matrix(rnorm(300), 100, 3))
  task_many <- list(matrix(rnorm(200), 20, 10 + 1))
  task_few <- list(matrix(rnorm(50), 25, 2))
  expect_equal(select_task_method(data, task_many), "gev")
  expect_equal(select_task_method(data, task_few), "augmented")
  expect_equal(select_task_method(data, list(matrix(rnorm(2000), 200, 5))), "blend")
})
