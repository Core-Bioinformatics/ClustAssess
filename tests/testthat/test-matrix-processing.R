test_that("make_pca_matrix_processing exact method matches prcomp embedding", {
    set.seed(1)
    dt_mtx <- matrix(stats::rnorm(10 * 5), nrow = 10)
    rownames(dt_mtx) <- paste0("cell_", seq_len(nrow(dt_mtx)))
    colnames(dt_mtx) <- paste0("feature_", seq_len(ncol(dt_mtx)))

    processor <- make_pca_matrix_processing(actual_npcs = 3, method = "exact")
    observed <- processor(dt_mtx)
    expected <- stats::prcomp(dt_mtx, rank. = 3)$x[, seq_len(3), drop = FALSE]
    colnames(expected) <- paste0("PC_", seq_len(ncol(expected)))

    expect_equal(formals(processor)$actual_npcs, 3L)
    expect_equal(observed, expected)
})

test_that("make_pca_matrix_processing irlba method returns named embedding", {
    testthat::skip_if_not_installed("irlba")

    set.seed(1)
    dt_mtx <- matrix(stats::rnorm(50 * 30), nrow = 50)
    rownames(dt_mtx) <- paste0("cell_", seq_len(nrow(dt_mtx)))
    colnames(dt_mtx) <- paste0("feature_", seq_len(ncol(dt_mtx)))

    processor <- make_pca_matrix_processing(actual_npcs = 5, method = "irlba")
    observed <- processor(dt_mtx)

    expect_equal(dim(observed), c(50L, 5L))
    expect_equal(rownames(observed), rownames(dt_mtx))
    expect_equal(colnames(observed), paste0("PC_", seq_len(5)))
    expect_equal(formals(processor)$actual_npcs, 5L)
})

test_that("make_pca_matrix_processing irlba method is deterministic and preserves RNG state", {
    testthat::skip_if_not_installed("irlba")

    set.seed(123)
    dt_mtx <- matrix(stats::rnorm(60 * 40), nrow = 60)
    rownames(dt_mtx) <- paste0("cell_", seq_len(nrow(dt_mtx)))
    colnames(dt_mtx) <- paste0("feature_", seq_len(ncol(dt_mtx)))

    processor <- make_pca_matrix_processing(
        actual_npcs = 6,
        method = "irlba",
        seed = 99
    )

    set.seed(1)
    before_seed <- .Random.seed
    first <- processor(dt_mtx)
    after_seed <- .Random.seed

    set.seed(999)
    second <- processor(dt_mtx)

    expect_identical(after_seed, before_seed)
    expect_equal(second, first)

    different_seed_processor <- make_pca_matrix_processing(
        actual_npcs = 6,
        method = "irlba",
        seed = 100
    )
    set.seed(1)
    third <- different_seed_processor(dt_mtx)

    expect_false(isTRUE(all.equal(third, first)))
})

test_that("make_pca_matrix_processing validates inputs", {
    expect_error(
        make_pca_matrix_processing(actual_npcs = 0),
        "actual_npcs parameter should be at least 1",
        fixed = TRUE
    )
    expect_error(
        make_pca_matrix_processing(actual_npcs = NA),
        "actual_npcs parameter should be a single numeric value",
        fixed = TRUE
    )
    expect_error(
        make_pca_matrix_processing(method = "irlba", seed = NA),
        "seed parameter should be a single numeric value",
        fixed = TRUE
    )

    processor <- make_pca_matrix_processing(actual_npcs = 2)
    expect_error(
        processor(data.frame(a = 1:3)),
        "dt_mtx parameter should be a matrix",
        fixed = TRUE
    )
})
