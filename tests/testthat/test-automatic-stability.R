tiny_automatic_expression_matrix <- function() {
    expression_matrix <- matrix(
        c(
            5, 5, 0, 0, 0, 0, 0, 0,
            4, 4, 0, 0, 0, 0, 0, 0,
            0, 0, 5, 5, 0, 0, 0, 0,
            0, 0, 4, 4, 0, 0, 0, 0,
            0, 0, 0, 0, 5, 5, 0, 0,
            0, 0, 0, 0, 4, 4, 0, 0,
            0, 0, 0, 0, 0, 0, 5, 5,
            0, 0, 0, 0, 0, 0, 4, 4
        ),
        nrow = 8,
        byrow = TRUE
    )
    rownames(expression_matrix) <- paste0("gene_", seq_len(nrow(expression_matrix)))
    colnames(expression_matrix) <- paste0("cell_", seq_len(ncol(expression_matrix)))
    expression_matrix
}

tiny_automatic_matrix_processing <- function(dt_mtx, actual_npcs = 2, ...) {
    npc <- min(actual_npcs, nrow(dt_mtx) - 1, ncol(dt_mtx))
    embedding <- stats::prcomp(dt_mtx, rank. = npc)$x[, seq_len(npc), drop = FALSE]
    rownames(embedding) <- rownames(dt_mtx)
    colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))
    embedding
}

run_tiny_automatic_stability <- function(force_gc = FALSE,
                                         compute_umap_embeddings = TRUE,
                                         return_nn_idx = TRUE,
                                         matrix_processing = tiny_automatic_matrix_processing) {
    expression_matrix <- tiny_automatic_expression_matrix()
    features_sets <- list(all = rownames(expression_matrix))
    steps <- list(all = c(4, 8))

    old_options <- options(ClustAssess.force_gc = force_gc)
    on.exit(options(old_options), add = TRUE)

    suppressWarnings(
        automatic_stability_assessment(
            expression_matrix = expression_matrix,
            n_repetitions = 2,
            n_neigh_sequence = c(2),
            resolution_sequence = c(0.5),
            features_sets = features_sets,
            steps = steps,
            seed_sequence = c(1, 101),
            n_top_configs = 1,
            matrix_processing = matrix_processing,
            compute_umap_embeddings = compute_umap_embeddings,
            return_nn_idx = return_nn_idx,
            verbose = FALSE,
            save_temp = FALSE
        )
    )
}

test_that("automatic_stability_assessment returns complete tiny pipeline object", {
    testthat::skip_if_not_installed("Seurat")

    result <- run_tiny_automatic_stability(force_gc = FALSE)

    expect_named(result, c("feature_stability", "all"))
    expect_true(length(result$all) >= 1)

    config <- result$all[[1]]
    expect_named(
        config,
        c(
            "pca", "umap", "stable_config", "nn_conn_comps",
            "nn_stability", "adj_matrix", "clustering_stability"
        )
    )
    expect_named(
        config$stable_config,
        c(
            "feature_set", "n_features", "n_pcs", "base_embedding",
            "graph_type", "n_neighbours", "prune_param"
        )
    )
    expect_equal(config$stable_config$feature_set, "all")
    expect_true(is.matrix(config$pca))
    expect_true(is.matrix(config$umap))
    expect_true(inherits(config$adj_matrix, "sparseMatrix"))
    expect_named(config$clustering_stability, c("split_by_resolution", "split_by_k"))
})

test_that("automatic_stability_assessment accepts PCA matrix-processing factory", {
    testthat::skip_if_not_installed("Seurat")

    expression_matrix <- tiny_automatic_expression_matrix()
    result <- suppressWarnings(
        automatic_stability_assessment(
            expression_matrix = expression_matrix,
            n_repetitions = 2,
            n_neigh_sequence = c(2),
            resolution_sequence = c(0.5),
            features_sets = list(all = rownames(expression_matrix)),
            steps = list(all = c(4, 8)),
            seed_sequence = c(1, 101),
            n_top_configs = 1,
            matrix_processing = make_pca_matrix_processing(actual_npcs = 2),
            compute_umap_embeddings = FALSE,
            return_nn_idx = FALSE,
            verbose = FALSE,
            save_temp = FALSE
        )
    )

    config <- result$all[[1]]
    expect_equal(config$stable_config$n_pcs, 2)
    expect_false("nn_idx" %in% names(config$nn_stability))
    expect_null(config$umap)
})

test_that("automatic stability records actual processor dimensions for small panels", {
    testthat::skip_if_not_installed("Seurat")

    processors <- list(
        factory = make_pca_matrix_processing(actual_npcs = 30),
        custom = function(dt_mtx, actual_npcs = 30, ...) {
            tiny_automatic_matrix_processing(dt_mtx, actual_npcs = 3)
        },
        without_npcs_formal = function(dt_mtx) {
            tiny_automatic_matrix_processing(dt_mtx, actual_npcs = 1)
        }
    )

    for (processor_name in names(processors)) {
        result <- run_tiny_automatic_stability(
            compute_umap_embeddings = FALSE,
            return_nn_idx = FALSE,
            matrix_processing = processors[[processor_name]]
        )
        config <- result$all[[1]]
        expected_npcs <- switch(
            processor_name,
            factory = min(7L, config$stable_config$n_features),
            custom = 3L,
            without_npcs_formal = 1L
        )

        expect_equal(nrow(config$pca), 8L, info = processor_name)
        expect_equal(ncol(config$pca), expected_npcs, info = processor_name)
        expect_equal(
            config$stable_config$n_pcs,
            ncol(config$pca),
            info = processor_name
        )
    }
})

test_that("automatic_stability_assessment behavior is unchanged by optional GC mode", {
    testthat::skip_if_not_installed("Seurat")

    without_gc <- run_tiny_automatic_stability(force_gc = FALSE)
    with_gc <- run_tiny_automatic_stability(force_gc = TRUE)

    without_config <- without_gc$all[[1]]
    with_config <- with_gc$all[[1]]

    expect_equal(with_config$stable_config, without_config$stable_config)
    expect_equal(with_config$adj_matrix, without_config$adj_matrix)
    expect_equal(with_config$nn_stability, without_config$nn_stability)
    expect_equal(with_config$clustering_stability, without_config$clustering_stability)
})

test_that("optional UMAP embedding storage leaves automatic stability results unchanged", {
    testthat::skip_if_not_installed("Seurat")

    with_embeddings <- run_tiny_automatic_stability(compute_umap_embeddings = TRUE)
    without_embeddings <- run_tiny_automatic_stability(compute_umap_embeddings = FALSE)

    with_config <- with_embeddings$all[[1]]
    without_config <- without_embeddings$all[[1]]
    with_config$umap <- NULL
    without_config$umap <- NULL

    expect_equal(without_config, with_config)
    expect_equal(
        without_embeddings$feature_stability$by_steps,
        with_embeddings$feature_stability$by_steps
    )
    expect_equal(
        without_embeddings$feature_stability$incremental,
        with_embeddings$feature_stability$incremental
    )
    expect_null(without_embeddings$feature_stability$embedding_list$all[["4"]])
    expect_null(without_embeddings$feature_stability$embedding_list$all[["8"]])
    expect_true(is.matrix(with_embeddings$feature_stability$embedding_list$all[["4"]]))
    expect_true(is.matrix(with_embeddings$feature_stability$embedding_list$all[["8"]]))
})

test_that("automatic stability can omit returned nn_idx without changing results", {
    testthat::skip_if_not_installed("Seurat")

    with_nn_idx <- run_tiny_automatic_stability(return_nn_idx = TRUE)
    without_nn_idx <- run_tiny_automatic_stability(return_nn_idx = FALSE)

    with_config <- with_nn_idx$all[[1]]
    without_config <- without_nn_idx$all[[1]]
    expect_false("nn_idx" %in% names(without_config$nn_stability))
    expect_true("nn_idx" %in% names(with_config$nn_stability))

    with_config$nn_stability$nn_idx <- NULL
    expect_equal(without_config, with_config)
    expect_equal(without_nn_idx$feature_stability, with_nn_idx$feature_stability)
})

test_that("optional UMAP embedding flag preserves exported positional argument order", {
    feature_args <- names(formals(assess_feature_stability))
    feature_legacy_tail <- c(
        "prune_value", "clustering_algorithm", "clustering_arguments", "verbose"
    )
    expect_equal(
        feature_args[
            seq(
                match("umap_arguments", feature_args) + 1,
                match("verbose", feature_args)
            )
        ],
        feature_legacy_tail
    )
    expect_gt(
        match("compute_umap_embeddings", feature_args),
        match("verbose", feature_args)
    )

    automatic_args <- names(formals(automatic_stability_assessment))
    automatic_legacy_tail <- c(
        "prune_value", "algorithm_dim_reduction", "algorithm_graph_construct",
        "algorithms_clustering_assessment", "clustering_arguments", "verbose",
        "temp_file", "save_temp"
    )
    expect_equal(
        automatic_args[
            seq(
                match("umap_arguments", automatic_args) + 1,
                match("save_temp", automatic_args)
            )
        ],
        automatic_legacy_tail
    )
    expect_gt(
        match("compute_umap_embeddings", automatic_args),
        match("save_temp", automatic_args)
    )
    expect_gt(
        match("return_nn_idx", automatic_args),
        match("compute_umap_embeddings", automatic_args)
    )
})

test_that("assess_feature_stability can skip optional UMAP embeddings without changing scores", {
    testthat::skip_if_not_installed("Seurat")

    expression_matrix <- tiny_automatic_expression_matrix()
    with_embeddings <- suppressWarnings(assess_feature_stability(
        data_matrix = expression_matrix,
        feature_set = rownames(expression_matrix),
        steps = c(4),
        feature_type = "all",
        resolution = c(0.5),
        n_repetitions = 2,
        seed_sequence = c(1, 101),
        matrix_processing = tiny_automatic_matrix_processing,
        compute_umap_embeddings = TRUE,
        verbose = FALSE
    ))
    without_embeddings <- suppressWarnings(assess_feature_stability(
        data_matrix = expression_matrix,
        feature_set = rownames(expression_matrix),
        steps = c(4),
        feature_type = "all",
        resolution = c(0.5),
        n_repetitions = 2,
        seed_sequence = c(1, 101),
        matrix_processing = tiny_automatic_matrix_processing,
        compute_umap_embeddings = FALSE,
        verbose = FALSE
    ))

    expect_equal(without_embeddings$by_steps, with_embeddings$by_steps)
    expect_equal(without_embeddings$incremental, with_embeddings$incremental)
    expect_equal(without_embeddings$pca_list, with_embeddings$pca_list)
    expect_null(without_embeddings$embedding_list$all[["4"]])
    expect_true(is.matrix(with_embeddings$embedding_list$all[["4"]]))
})
