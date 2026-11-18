test_that("a last chunk of one gene keeps its row", {
    set.seed(1)
    dense <- matrix(rpois(101 * 6, 3), nrow = 101, dimnames = list(paste0("g", 1:101), paste0("c", 1:6)))
    dense[7, ] <- 2 # no variance
    dense[101, ] <- c(0, 5, 1, 0, 9, 2)
    expected <- dense[matrixStats::rowVars(dense) > 0, , drop = FALSE]
    expect_true("g101" %in% rownames(expected))
    expect_false("g7" %in% rownames(expected))

    for (input in list(dense, Matrix::Matrix(dense, sparse = TRUE))) {
        result <- NULL
        expect_output(result <- ClustAssess:::dense_rows_above_variance(input, 100, 0))
        expect_true(is.matrix(result))
        expect_identical(dimnames(result), dimnames(expected))
        expect_equal(unname(result), unname(expected))
    }
})

test_that("a last chunk of one gene without variance is left out, as any other gene", {
    dense <- matrix(as.numeric(seq_len(101 * 4)), nrow = 101, dimnames = list(paste0("g", 1:101), paste0("c", 1:4)))
    dense[101, ] <- 1
    result <- NULL
    expect_output(result <- ClustAssess:::dense_rows_above_variance(Matrix::Matrix(dense, sparse = TRUE), 100, 0))
    expect_identical(rownames(result), paste0("g", 1:100))
    expect_equal(unname(result), unname(dense[1:100, ]))
})

test_that("chunks of several genes give the rows of one pass", {
    set.seed(2)
    dense <- matrix(rpois(250 * 5, 1), nrow = 250, dimnames = list(paste0("g", 1:250), paste0("c", 1:5)))
    expected <- dense[matrixStats::rowVars(dense) > 0.5, , drop = FALSE]
    for (size in c(1, 7, 100, 250, 1000)) {
        result <- NULL
        expect_output(result <- ClustAssess:::dense_rows_above_variance(Matrix::Matrix(dense, sparse = TRUE), size, 0.5))
        expect_identical(dimnames(result), dimnames(expected))
        expect_equal(unname(result), unname(expected))
    }
})
