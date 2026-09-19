test_that("getNNmatrix builds the expected directed NN adjacency", {
    nn_idx <- matrix(
        c(
            1, 2, 3, 4,
            1, 2, 4, 3,
            2, 3, 4, 1
        ),
        nrow = 4,
        byrow = TRUE
    )

    nn <- getNNmatrix(nn_idx, k = 2, start = 0, prune = -1)$nn

    expected <- matrix(
        c(
            1, 1, 0, 0,
            1, 0, 0, 0,
            0, 0, 1, 1,
            0, 1, 1, 1
        ),
        nrow = 4
    )

    expect_equal(as.matrix(nn), expected)
})

test_that("computeSNN normalizes shared-neighbor counts and prunes thresholded entries", {
    nn_idx <- matrix(
        c(
            1, 2, 3, 4,
            1, 2, 4, 3,
            2, 3, 4, 1
        ),
        nrow = 4,
        byrow = TRUE
    )
    nn <- getNNmatrix(nn_idx, k = 2, start = 0, prune = -1)$nn

    snn <- computeSNN(nn, k = 2, prune = 0)
    pruned <- computeSNN(nn, k = 2, prune = 0.5)

    expected_snn <- matrix(
        c(
            1, 1 / 3, 0, 0,
            1 / 3, 1, 1 / 3, 1 / 3,
            0, 1 / 3, 1, 1,
            0, 1 / 3, 1, 1
        ),
        nrow = 4
    )
    expected_pruned <- matrix(
        c(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 1,
            0, 0, 1, 1
        ),
        nrow = 4
    )

    expect_equal(as.matrix(snn), expected_snn)
    expect_equal(as.matrix(pruned), expected_pruned)
})

test_that("parallel_nn2_idx matches RANN::nn2 neighbor indices", {
    testthat::skip_if_not_installed("RANN")

    embedding <- matrix(
        c(
            0, 0,
            0, 0.1,
            5, 5,
            5.1, 5,
            10, 10,
            10.1, 10
        ),
        ncol = 2,
        byrow = TRUE
    )

    expected <- RANN::nn2(embedding, k = 3)$nn.idx
    observed <- ClustAssess:::parallel_nn2_idx(embedding, k = 3)

    expect_equal(observed, expected)
})

test_that("build_neigh_matrix matches existing NN and fixed-prune SNN primitives", {
    nn_idx <- matrix(
        c(
            1, 2, 3, 4,
            1, 2, 4, 3,
            2, 3, 4, 1
        ),
        nrow = 4,
        byrow = TRUE
    )
    cell_names <- paste0("cell_", seq_len(nrow(nn_idx)))

    observed <- ClustAssess:::build_neigh_matrix(
        nn_idx = nn_idx,
        n_neigh = 2,
        cell_names = cell_names,
        graph_type = 2,
        prune_value = 0.5
    )
    expected_nn <- getNNmatrix(nn_idx, k = 2, start = 0, prune = -1)$nn
    expected_snn <- computeSNN(expected_nn, k = 2, prune = 0.5)
    rownames(expected_nn) <- cell_names
    colnames(expected_nn) <- cell_names
    rownames(expected_snn) <- cell_names
    colnames(expected_snn) <- cell_names

    expect_equal(observed$nn, expected_nn)
    expect_equal(observed$snn, expected_snn)
})

test_that("build_neigh_matrix respects graph type and automatic pruning behavior", {
    nn_idx <- matrix(
        c(
            1, 2, 3, 4,
            1, 2, 4, 3,
            2, 3, 4, 1
        ),
        nrow = 4,
        byrow = TRUE
    )
    cell_names <- paste0("cell_", seq_len(nrow(nn_idx)))

    nn_only <- ClustAssess:::build_neigh_matrix(
        nn_idx = nn_idx,
        n_neigh = 2,
        cell_names = cell_names,
        graph_type = 0,
        prune_value = -1
    )
    auto_pruned <- ClustAssess:::build_neigh_matrix(
        nn_idx = nn_idx,
        n_neigh = 2,
        cell_names = cell_names,
        graph_type = 1,
        prune_value = -1
    )

    expected_nn <- getNNmatrix(nn_idx, k = 2, start = 0, prune = -1)$nn
    expected_snn <- get_highest_prune_param(expected_nn, n_neigh = 2)$adj_matrix
    rownames(expected_nn) <- cell_names
    colnames(expected_nn) <- cell_names
    rownames(expected_snn) <- cell_names
    colnames(expected_snn) <- cell_names

    expect_false("snn" %in% names(nn_only))
    expect_equal(nn_only$nn, expected_nn)
    expect_equal(auto_pruned$nn, expected_nn)
    expect_equal(auto_pruned$snn, expected_snn)
})

test_that("build_neigh_matrix can omit unused NN matrix for SNN-only callers", {
    nn_idx <- matrix(
        c(
            1, 2, 3, 4, 5,
            1, 3, 2, 4, 5,
            2, 3, 4, 5, 1,
            3, 4, 5, 1, 2,
            4, 5, 3, 2, 1
        ),
        nrow = 5,
        byrow = TRUE
    )
    cell_names <- paste0("cell_", seq_len(nrow(nn_idx)))

    retained <- ClustAssess:::build_neigh_matrix(
        nn_idx = nn_idx,
        n_neigh = 3,
        cell_names = cell_names,
        graph_type = 1,
        prune_value = 0.25
    )
    omitted <- ClustAssess:::build_neigh_matrix(
        nn_idx = nn_idx,
        n_neigh = 3,
        cell_names = cell_names,
        graph_type = 1,
        prune_value = 0.25,
        keep_nn = FALSE
    )

    expect_named(omitted, "snn")
    expect_equal(omitted$snn, retained$snn)
})

test_that("automatic SNN pruning is unchanged by optional GC mode", {
    nn_idx <- matrix(
        c(
            1, 2, 3, 4,
            1, 2, 4, 3,
            2, 3, 4, 1
        ),
        nrow = 4,
        byrow = TRUE
    )
    nn_matrix <- getNNmatrix(nn_idx, k = 2, start = 0, prune = -1)$nn

    old_options <- options(ClustAssess.force_gc = FALSE)
    on.exit(options(old_options), add = TRUE)
    without_gc <- get_highest_prune_param(nn_matrix, n_neigh = 2)

    options(ClustAssess.force_gc = TRUE)
    with_gc <- get_highest_prune_param(nn_matrix, n_neigh = 2)

    expect_equal(with_gc$prune_value, without_gc$prune_value)
    expect_equal(with_gc$adj_matrix, without_gc$adj_matrix)
})

test_that("assess_nn_stability returns stable object structure for PCA graphs", {
    testthat::skip_if_not_installed("Seurat")

    embedding <- matrix(
        c(
            0, 0,
            0, 0.1,
            5, 5,
            5.1, 5,
            10, 10,
            10.1, 10,
            15, 15,
            15.1, 15
        ),
        ncol = 2,
        byrow = TRUE
    )
    rownames(embedding) <- paste0("cell_", seq_len(nrow(embedding)))
    colnames(embedding) <- c("PC_1", "PC_2")

    result <- suppressWarnings(
        assess_nn_stability(
            embedding = embedding,
            n_neigh_sequence = c(2),
            seed_sequence = c(1, 101),
            graph_reduction_type = "PCA",
            graph_type = 1,
            clustering_algorithm = 1
        )
    )

    expect_named(
        result,
        c("n_neigh_k_corresp", "n_neigh_ec_consistency", "n_different_partitions", "nn_idx")
    )
    expect_equal(dim(result$nn_idx), c(nrow(embedding), 2L))
    expect_named(result$n_neigh_ec_consistency, "PCA_snn")
    expect_named(result$n_neigh_ec_consistency$PCA_snn, "2")
    expect_true(is.numeric(result$n_neigh_ec_consistency$PCA_snn[["2"]]))
})

test_that("assess_nn_stability return_nn_idx flag preserves positional argument order", {
    stability_args <- names(formals(assess_nn_stability))
    legacy_tail <- c(
        "prune_value", "clustering_algorithm", "clustering_arguments",
        "umap_arguments"
    )
    expect_equal(
        stability_args[
            seq(
                match("graph_type", stability_args) + 1,
                match("umap_arguments", stability_args)
            )
        ],
        legacy_tail
    )
    expect_gt(
        match("return_nn_idx", stability_args),
        match("umap_arguments", stability_args)
    )
})

test_that("assess_nn_stability PCA graph_type controls returned graph families", {
    testthat::skip_if_not_installed("Seurat")

    embedding <- matrix(
        c(
            0, 0,
            0, 0.1,
            5, 5,
            5.1, 5,
            10, 10,
            10.1, 10,
            15, 15,
            15.1, 15
        ),
        ncol = 2,
        byrow = TRUE
    )
    rownames(embedding) <- paste0("cell_", seq_len(nrow(embedding)))
    colnames(embedding) <- c("PC_1", "PC_2")

    run_pca <- function(graph_type) {
        suppressWarnings(
            assess_nn_stability(
                embedding = embedding,
                n_neigh_sequence = c(2),
                seed_sequence = c(1, 101),
                graph_reduction_type = "PCA",
                graph_type = graph_type,
                clustering_algorithm = 1
            )
        )
    }

    nn_only <- run_pca(graph_type = 0)
    snn_only <- run_pca(graph_type = 1)
    both <- run_pca(graph_type = 2)

    expect_named(nn_only$n_neigh_ec_consistency, "PCA_nn")
    expect_named(nn_only$n_neigh_k_corresp, "PCA_nn")
    expect_named(nn_only$n_different_partitions, "PCA_nn")
    expect_named(snn_only$n_neigh_ec_consistency, "PCA_snn")
    expect_named(snn_only$n_neigh_k_corresp, "PCA_snn")
    expect_named(snn_only$n_different_partitions, "PCA_snn")
    expect_named(both$n_neigh_ec_consistency, c("PCA_snn", "PCA_nn"))
    expect_named(both$n_neigh_k_corresp, c("PCA_snn", "PCA_nn"))
    expect_named(both$n_different_partitions, c("PCA_snn", "PCA_nn"))
})

test_that("assess_nn_stability can omit returned PCA nn_idx without changing summaries", {
    testthat::skip_if_not_installed("Seurat")

    embedding <- matrix(
        c(
            0, 0,
            0, 0.1,
            5, 5,
            5.1, 5,
            10, 10,
            10.1, 10,
            15, 15,
            15.1, 15
        ),
        ncol = 2,
        byrow = TRUE
    )
    rownames(embedding) <- paste0("cell_", seq_len(nrow(embedding)))
    colnames(embedding) <- c("PC_1", "PC_2")

    with_nn_idx <- suppressWarnings(assess_nn_stability(
        embedding = embedding,
        n_neigh_sequence = c(2, 3),
        seed_sequence = c(1, 101),
        graph_reduction_type = "PCA",
        graph_type = 1,
        clustering_algorithm = 1,
        return_nn_idx = TRUE
    ))
    without_nn_idx <- suppressWarnings(assess_nn_stability(
        embedding = embedding,
        n_neigh_sequence = c(2, 3),
        seed_sequence = c(1, 101),
        graph_reduction_type = "PCA",
        graph_type = 1,
        clustering_algorithm = 1,
        return_nn_idx = FALSE
    ))

    expect_true("nn_idx" %in% names(with_nn_idx))
    expect_false("nn_idx" %in% names(without_nn_idx))
    with_nn_idx$nn_idx <- NULL
    expect_equal(without_nn_idx, with_nn_idx)
})

test_that("assess_nn_stability PCA works with PSOCK workers when package-installed", {
    testthat::skip_if_not_installed("ClustAssess")
    testthat::skip_if_not_installed("doParallel")
    testthat::skip_if_not_installed("Seurat")

    package_path <- system.file(package = "ClustAssess")
    testthat::skip_if(!nzchar(package_path), "ClustAssess is not installed in a library path")
    testthat::skip_if(
        !file.exists(file.path(package_path, "Meta", "package.rds")),
        "ClustAssess is loaded from source, not an installed package"
    )

    embedding <- matrix(
        c(
            0, 0,
            0, 0.1,
            5, 5,
            5.1, 5,
            10, 10,
            10.1, 10,
            15, 15,
            15.1, 15
        ),
        ncol = 2,
        byrow = TRUE
    )
    rownames(embedding) <- paste0("cell_", seq_len(nrow(embedding)))
    colnames(embedding) <- c("PC_1", "PC_2")

    cl <- parallel::makePSOCKcluster(2)
    on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
    }, add = TRUE)
    parallel::clusterCall(cl, function(paths) {
        .libPaths(paths)
        NULL
    }, .libPaths())
    doParallel::registerDoParallel(cl)

    result <- suppressWarnings(
        ClustAssess::assess_nn_stability(
            embedding = embedding,
            n_neigh_sequence = c(2),
            seed_sequence = c(1, 101),
            graph_reduction_type = "PCA",
            graph_type = 1,
            clustering_algorithm = 1
        )
    )

    expect_named(result$n_neigh_ec_consistency, "PCA_snn")
    expect_named(result$n_neigh_ec_consistency$PCA_snn, "2")
})

test_that("assess_nn_stability PCA multi-neighbor results match separate runs", {
    testthat::skip_if_not_installed("Seurat")

    embedding <- matrix(
        c(
            0, 0,
            0, 0.1,
            5, 5,
            5.1, 5,
            10, 10,
            10.1, 10,
            15, 15,
            15.1, 15
        ),
        ncol = 2,
        byrow = TRUE
    )
    rownames(embedding) <- paste0("cell_", seq_len(nrow(embedding)))
    colnames(embedding) <- c("PC_1", "PC_2")

    run_nn <- function(n_neigh_sequence) {
        suppressWarnings(
            assess_nn_stability(
                embedding = embedding,
                n_neigh_sequence = n_neigh_sequence,
                seed_sequence = c(1, 101),
                graph_reduction_type = "PCA",
                graph_type = 1,
                clustering_algorithm = 1
            )
        )
    }

    combined <- run_nn(c(2, 3))
    separate_2 <- run_nn(c(2))
    separate_3 <- run_nn(c(3))

    expect_equal(
        combined$n_neigh_ec_consistency$PCA_snn[["2"]],
        separate_2$n_neigh_ec_consistency$PCA_snn[["2"]]
    )
    expect_equal(
        combined$n_neigh_ec_consistency$PCA_snn[["3"]],
        separate_3$n_neigh_ec_consistency$PCA_snn[["3"]]
    )
    expect_equal(
        combined$n_neigh_k_corresp$PCA_snn[["2"]],
        separate_2$n_neigh_k_corresp$PCA_snn[["2"]]
    )
    expect_equal(
        combined$n_neigh_k_corresp$PCA_snn[["3"]],
        separate_3$n_neigh_k_corresp$PCA_snn[["3"]]
    )
    expect_equal(
        combined$n_different_partitions$PCA_snn[["2"]],
        separate_2$n_different_partitions$PCA_snn[["2"]]
    )
    expect_equal(
        combined$n_different_partitions$PCA_snn[["3"]],
        separate_3$n_different_partitions$PCA_snn[["3"]]
    )
})

test_that("assess_nn_stability UMAP custom seed_sequence ignores n_repetitions", {
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if_not_installed("uwot")

    set.seed(42)
    embedding <- matrix(stats::rnorm(24 * 4), ncol = 4)
    rownames(embedding) <- paste0("cell_", seq_len(nrow(embedding)))
    colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))

    run_umap <- function(graph_type = 1, ...) {
        suppressWarnings(
            assess_nn_stability(
                embedding = embedding,
                n_neigh_sequence = c(2, 3),
                seed_sequence = c(1, 101),
                graph_reduction_type = "UMAP",
                graph_type = graph_type,
                clustering_algorithm = 1,
                umap_arguments = list(
                    n_neighbors = 5,
                    n_epochs = 20,
                    n_threads = 1,
                    n_sgd_threads = 1
                ),
                ...
            )
        )
    }

    implicit_repetitions <- run_umap()
    explicit_repetitions <- run_umap(n_repetitions = 2)

    expect_named(implicit_repetitions$n_neigh_ec_consistency, "UMAP_snn")
    expect_equal(
        implicit_repetitions$n_neigh_ec_consistency,
        explicit_repetitions$n_neigh_ec_consistency
    )
    expect_equal(
        implicit_repetitions$n_neigh_k_corresp,
        explicit_repetitions$n_neigh_k_corresp
    )
    expect_equal(
        implicit_repetitions$n_different_partitions,
        explicit_repetitions$n_different_partitions
    )

    nn_result <- run_umap(graph_type = 0)
    expect_named(nn_result$n_neigh_ec_consistency, "UMAP_nn")
})

test_that("assess_nn_stability UMAP works with PSOCK workers when package-installed", {
    testthat::skip_if_not_installed("ClustAssess")
    testthat::skip_if_not_installed("doParallel")
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if_not_installed("uwot")

    package_path <- system.file(package = "ClustAssess")
    testthat::skip_if(!nzchar(package_path), "ClustAssess is not installed in a library path")
    testthat::skip_if(
        !file.exists(file.path(package_path, "Meta", "package.rds")),
        "ClustAssess is loaded from source, not an installed package"
    )

    set.seed(42)
    embedding <- matrix(stats::rnorm(24 * 4), ncol = 4)
    rownames(embedding) <- paste0("cell_", seq_len(nrow(embedding)))
    colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))

    cl <- parallel::makePSOCKcluster(2)
    on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
    }, add = TRUE)
    parallel::clusterCall(cl, function(paths) {
        .libPaths(paths)
        NULL
    }, .libPaths())
    doParallel::registerDoParallel(cl)

    result <- suppressWarnings(
        ClustAssess::assess_nn_stability(
            embedding = embedding,
            n_neigh_sequence = c(2, 3),
            seed_sequence = c(1, 101),
            graph_reduction_type = "UMAP",
            graph_type = 1,
            clustering_algorithm = 1,
            umap_arguments = list(
                n_neighbors = 5,
                n_epochs = 20,
                n_threads = 1,
                n_sgd_threads = 1
            )
        )
    )

    expect_named(result$n_neigh_ec_consistency, "UMAP_snn")
    expect_named(result$n_neigh_ec_consistency$UMAP_snn, c("2", "3"))
})
