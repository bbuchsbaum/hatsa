library(testthat)

describe("build_augmented_anchor_matrix", {
  test_that("combines rows and names correctly", {
    A <- matrix(1:6, nrow=3, ncol=2,
                 dimnames = list(paste0("p",1:3), paste0("k",1:2)))
    Z <- matrix(7:10, nrow=2, ncol=2,
                 dimnames = list(paste0("t",1:2), paste0("k",1:2)))
    res <- build_augmented_anchor_matrix(A, Z)
    expect_equal(dim(res), c(5,2))
    expect_equal(rownames(res), c("p1","p2","p3","t1","t2"))
    expect_equal(colnames(res), colnames(A))
    expect_equal(res[1:3, ], A)
    expect_equal(res[4:5, ], Z)
  })

  test_that("errors on mismatched column names or counts", {
    A <- matrix(1:4, nrow=2, ncol=2, dimnames=list(NULL,c("a","b")))
    Z_bad <- matrix(5:8, nrow=2, ncol=2, dimnames=list(NULL,c("b","a")))
    expect_error(build_augmented_anchor_matrix(A, Z_bad), "Column names")
    Z_bad_dim <- matrix(5:9, nrow=2, ncol=3)
    expect_error(build_augmented_anchor_matrix(A, Z_bad_dim), "Number of columns")
  })

  test_that("handles NULL or empty inputs", {
    A <- matrix(1:4, nrow=2, ncol=2)
    res <- build_augmented_anchor_matrix(A, NULL)
    expect_identical(res, A)

    A0 <- matrix(numeric(0), nrow=0, ncol=2)
    Z <- matrix(5:8, nrow=2, ncol=2)
    res2 <- build_augmented_anchor_matrix(A0, Z)
    expect_identical(res2, Z)
  })
})
