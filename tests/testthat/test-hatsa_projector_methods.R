describe("hatsa_projector basic methods", {

  # helper to generate small subject data
  .gen_subject_data <- function(N, V, T_len = 20) {
    lapply(1:N, function(i) matrix(rnorm(T_len * V), ncol = V))
  }

  .default_params <- function(V, k) {
    list(
      anchor_indices = 1:min(V,2),
      spectral_rank_k = k,
      k_conn_pos = 2,
      k_conn_neg = 2,
      n_refine = 1
    )
  }

  test_that("constructor builds valid object and checks dimensions", {
    V <- 4; N <- 2; k <- 2
    data_list <- .gen_subject_data(N, V)
    params <- .default_params(V, k)
    proj <- suppressMessages(run_hatsa_core(
      subject_data_list = data_list,
      anchor_indices   = params$anchor_indices,
      spectral_rank_k  = params$spectral_rank_k,
      k_conn_pos       = params$k_conn_pos,
      k_conn_neg       = params$k_conn_neg,
      n_refine         = params$n_refine
    ))
    expect_s3_class(proj, "hatsa_projector")
    expect_equal(dim(proj$v), c(V, k))
    expect_equal(dim(proj$s), c(N * V, k))
    expect_length(proj$block_indices, N)

    bad_results <- list(
      U_aligned_list = list(matrix(1, V, k), matrix(1, V+1, k)),
      R_final_list = replicate(N, diag(k), simplify = FALSE),
      U_original_list = replicate(N, matrix(1, V, k), simplify = FALSE),
      Lambda_original_list = replicate(N, rep(1, k), simplify = FALSE),
      Lambda_original_gaps_list = replicate(N, rep(1, k-1), simplify = FALSE),
      T_anchor_final = matrix(1, length(params$anchor_indices), k)
    )
    expect_error(
      hatsa_projector(bad_results, list(k = k, N_subjects = N, V_p = V, method = "hatsa_core")),
      "inconsistent dimensions"
    )
  })

  test_that("predict.hatsa_projector handles new data and dimension checks", {
    V <- 4; N <- 2; k <- 2
    data_list <- .gen_subject_data(N, V)
    params <- .default_params(V, k)
    proj <- suppressMessages(run_hatsa_core(
      subject_data_list = data_list,
      anchor_indices   = params$anchor_indices,
      spectral_rank_k  = params$spectral_rank_k,
      k_conn_pos       = params$k_conn_pos,
      k_conn_neg       = params$k_conn_neg,
      n_refine         = params$n_refine
    ))

    new_list <- .gen_subject_data(1, V)
    pred <- suppressWarnings(predict(proj, newdata_list = new_list))
    expect_true(is.list(pred))
    expect_equal(length(pred), 1)
    expect_true(is.matrix(pred[[1]]))
    expect_equal(dim(pred[[1]]), c(V, k))

    wrong_list <- .gen_subject_data(1, V + 1)
    expect_warning(predict(proj, newdata_list = wrong_list), "model expects")
  })

  test_that("project_block.hatsa_projector retrieves stored block and errors", {
    V <- 4; N <- 2; k <- 2
    data_list <- .gen_subject_data(N, V)
    params <- .default_params(V, k)
    proj <- suppressMessages(run_hatsa_core(
      subject_data_list = data_list,
      anchor_indices   = params$anchor_indices,
      spectral_rank_k  = params$spectral_rank_k,
      k_conn_pos       = params$k_conn_pos,
      k_conn_neg       = params$k_conn_neg,
      n_refine         = params$n_refine
    ))

    stored <- project_block(proj, block = 1)
    expect_equal(dim(stored), c(V, k))

    expect_error(project_block(proj, block = N + 1), "Invalid block index")
  })

  test_that("summary.hatsa_projector returns structured list", {
    V <- 4; N <- 2; k <- 2
    data_list <- .gen_subject_data(N, V)
    params <- .default_params(V, k)
    proj <- suppressMessages(run_hatsa_core(
      subject_data_list = data_list,
      anchor_indices   = params$anchor_indices,
      spectral_rank_k  = params$spectral_rank_k,
      k_conn_pos       = params$k_conn_pos,
      k_conn_neg       = params$k_conn_neg,
      n_refine         = params$n_refine
    ))

    summ <- summary(proj)
    expect_s3_class(summ, "summary.hatsa_projector")
    expect_true(is.list(summ))
    expect_true("mean_anchor_alignment_error" %in% names(summ))
  })

  test_that("reconstruction_error.hatsa_projector computes anchor errors", {
    V <- 4; N <- 2; k <- 2
    data_list <- .gen_subject_data(N, V)
    params <- .default_params(V, k)
    proj <- suppressMessages(run_hatsa_core(
      subject_data_list = data_list,
      anchor_indices   = params$anchor_indices,
      spectral_rank_k  = params$spectral_rank_k,
      k_conn_pos       = params$k_conn_pos,
      k_conn_neg       = params$k_conn_neg,
      n_refine         = params$n_refine
    ))

    err <- reconstruction_error(proj, type = "anchors")
    expect_true(is.list(err))
    expect_equal(err$type, "anchors")
    expect_length(err$per_subject_error, N)

    expect_error(reconstruction_error(proj, type = "unknown"), "Unknown reconstruction error type")
  })

})
