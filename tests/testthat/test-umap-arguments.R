test_that("UMAP argument processing preserves an explicit neighbor count", {
    for (n_neighbors in c(3, 15, 30)) {
        args <- list(n_neighbors = n_neighbors, min_dist = 0.2, seed = 42L)
        original <- args
        processed <- ClustAssess:::process_umap_arguments(args, ncells = 50)

        expect_equal(processed$n_neighbors, n_neighbors)
        expect_equal(processed$min_dist, 0.2)
        expect_identical(processed$seed, 42L)
        expect_identical(args, original)
    }
})

test_that("UMAP neighbor defaults and cell-count bounds are preserved", {
    process <- ClustAssess:::process_umap_arguments

    expect_equal(process(list(), ncells = 50)$n_neighbors, 15)
    expect_equal(process(list(), ncells = 8)$n_neighbors, 7)
    expect_equal(process(list(n_neighbors = 30), ncells = 20)$n_neighbors, 19)
    expect_equal(process(list(n_neighbors = 7), ncells = 8)$n_neighbors, 7)
    expect_equal(process(list(n_neighbors = 3), ncells = 8)$n_neighbors, 3)

    processed <- process(list(n_threads = 4, n_sgd_threads = 2), ncells = 50)
    expect_equal(processed$n_threads, 1)
    expect_equal(processed$n_sgd_threads, 1)
})

test_that("UMAP connected components forward the requested neighbor count", {
    testthat::skip_if_not_installed("uwot")
    foreach::registerDoSEQ()
    observed <- list()
    testthat::local_mocked_bindings(
        umap = function(X, n_neighbors, n_threads, n_sgd_threads, ...) {
            observed[[length(observed) + 1L]] <<- list(
                n_neighbors = n_neighbors,
                n_threads = n_threads,
                n_sgd_threads = n_sgd_threads
            )
            X[, 1:2, drop = FALSE]
        },
        .package = "uwot"
    )
    embedding <- cbind(seq_len(20), (seq_len(20) %% 5)^2)

    result <- get_nn_conn_comps(
        embedding = embedding,
        n_neigh_sequence = c(2, 3),
        include_umap = TRUE,
        seed_sequence = c(1, 101),
        umap_arguments = list(n_neighbors = 3)
    )

    expect_length(observed, 2)
    expect_equal(observed, rep(list(list(
        n_neighbors = 3, n_threads = 1, n_sgd_threads = 1
    )), 2))
    expect_named(result, c("PCA", "UMAP"))
    expect_named(result$UMAP, c("2", "3"))
})
