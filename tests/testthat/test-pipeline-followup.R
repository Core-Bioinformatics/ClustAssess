mock_feature_cache <- function(data_matrix, feature_set, steps, feature_type,
                               matrix_processing, compute_umap_embeddings,
                               omit_pca = FALSE, ...) {
    steps <- as.integer(steps)
    steps[steps <= 0 | steps > length(feature_set)] <- length(feature_set)
    steps <- sort(unique(steps))
    keys <- as.character(steps)
    embedding <- matrix(seq_len(ncol(data_matrix) * 2), ncol = 2,
                        dimnames = list(colnames(data_matrix), c("PC_1", "PC_2")))
    by_steps <- lapply(seq_along(steps), function(i) {
        list(`0.5` = list(ecc = rep(if (i == 1) 1 else 0.5, ncol(data_matrix))))
    })
    names(by_steps) <- keys
    incremental <- setNames(rep(list(list(`0.5` = rep(1, ncol(data_matrix)))),
                                length(steps) - 1L), keys[-1])
    cached <- setNames(rep(list(embedding), length(keys)), keys)
    list(
        by_steps = setNames(list(by_steps), feature_type),
        incremental = setNames(list(incremental), feature_type),
        pca_list = if (omit_pca) NULL else setNames(list(cached), feature_type),
        embedding_list = setNames(list(if (compute_umap_embeddings) cached else list()),
                                  feature_type)
    )
}

local_mock_pipeline_graphs <- function(.local_envir = parent.frame()) {
    testthat::local_mocked_bindings(
        get_nn_conn_comps = function(...) list(PCA = list(`2` = 1)),
        assess_nn_stability = function(embedding, ...) {
            list(
                n_neigh_k_corresp = list(PCA_nn = list(`2` = 2)),
                n_neigh_ec_consistency = list(PCA_nn = list(`2` = rep(1, nrow(embedding)))),
                n_different_partitions = list(PCA_nn = c(`2` = 1)),
                nn_idx = RANN::nn2(embedding, k = 2)$nn.idx
            )
        },
        assess_clustering_stability = function(...) list(),
        .package = "ClustAssess", .env = .local_envir
    )
}

test_that("automatic feature counts share integer cache identities", {
    local_mock_pipeline_graphs()
    forwarded_steps <- NULL
    testthat::local_mocked_bindings(
        assess_feature_stability = function(..., steps) {
            forwarded_steps <<- steps
            mock_feature_cache(..., steps = steps)
        },
        .package = "ClustAssess"
    )
    expression <- matrix(0, nrow = 100002, ncol = 4,
                         dimnames = list(paste0("feature_", seq_len(100002)),
                                         paste0("cell_", seq_len(4))))
    for (steps in list(c(100000L, 100002L), c(100000, 100002),
                       c(1e5, 100002), c(100000.9, 100002))) {
        result <- automatic_stability_assessment(
            expression, 2, 2, 0.5, list(all = rownames(expression)), list(all = steps),
            n_top_configs = 1, save_temp = FALSE, verbose = FALSE,
            matrix_processing = function(x) stop("Unexpected cache miss")
        )
        expect_identical(forwarded_steps, c(100000L, 100002L))
        expect_named(result$all, c("100000", "feature_list"))
        config <- result$all[["100000"]]
        expect_identical(config$stable_config$n_features, 100000L)
        expect_equal(dim(config$pca), c(4L, 2L))
        expect_equal(config$umap, config$pca)
    }
})

test_that("automatic cache fallback passes a cell-by-feature matrix", {
    local_mock_pipeline_graphs()
    testthat::local_mocked_bindings(
        assess_feature_stability = function(...) mock_feature_cache(..., omit_pca = TRUE),
        .package = "ClustAssess"
    )
    dense <- matrix(seq_len(20), nrow = 5,
                    dimnames = list(paste0("feature_", seq_len(5)), paste0("cell_", seq_len(4))))
    for (expression in list(dense, methods::as(dense, "dgCMatrix"))) {
        for (count in c(1L, 3L)) {
            observed <- NULL
            result <- automatic_stability_assessment(
                expression, 2, 2, 0.5, list(all = rownames(expression)),
                list(all = c(count, 5L)), n_top_configs = 1,
                compute_umap_embeddings = FALSE, save_temp = FALSE, verbose = FALSE,
                matrix_processing = function(x) {
                    observed <<- x
                    matrix(seq_len(ncol(expression) * 2), ncol = 2,
                           dimnames = list(colnames(expression), c("PC_1", "PC_2")))
                }
            )
            expect_equal(dim(observed), c(4L, count))
            expect_equal(rownames(observed), colnames(expression))
            expect_equal(colnames(observed), rownames(expression)[seq_len(count)])
            expect_equal(as.matrix(observed), t(dense[seq_len(count), , drop = FALSE]))
            expect_equal(result$all[[as.character(count)]]$stable_config$n_features, count)
        }
    }
})

test_that("feature entrypoints preserve fractional steps and matrix dimensions", {
    expression <- matrix(seq_len(20), nrow = 5,
                         dimnames = list(paste0("feature_", seq_len(5)), paste0("cell_", seq_len(4))))
    for (step in c(0.5, 1, 2.9, -1, 6)) {
        expected_count <- if (step < 1 || step > 5) 5L else as.integer(step)
        observed <- NULL
        expect_error(assess_feature_stability(
            expression, rownames(expression), step, "all", 0.5,
            compute_umap_embeddings = FALSE,
            matrix_processing = function(x) {
                observed <<- x
                stop("processor reached")
            }
        ), "processor reached", fixed = TRUE)
        expect_equal(dim(observed), c(4L, expected_count))
        expect_equal(rownames(observed), colnames(expression))
        expect_equal(colnames(observed), rownames(expression)[seq_len(expected_count)])
    }
    expect_identical(normalize_feature_steps(c(0.5, 2.1, 2.9, 5, 6), 5L), c(2L, 5L))
})

for (pipeline_name in c("assess_feature_stability", "automatic_stability_assessment")) {
    test_that(paste(pipeline_name, "default PCA keeps a one-feature component"), {
        foreach::registerDoSEQ()
        pipeline <- get(pipeline_name)
        processor <- eval(formals(pipeline)$matrix_processing, envir = environment(pipeline))
        input <- matrix(sin(seq_len(48)), nrow = 8,
                        dimnames = list(paste0("cell_", seq_len(8)), paste0("gene_", seq_len(6))))
        expect_equal(ncol(processor(input, actual_npcs = 1)), 1L)
        for (invalid in c(0, -1)) {
            expect_error(processor(input[, 1, drop = FALSE], actual_npcs = invalid))
        }
        for (count in c(2L, 5L, 6L, 1L)) {
            subset <- input[, seq_len(count), drop = FALSE]
            observed <- processor(subset)
            expected <- stats::prcomp(subset, rank. = max(1L, count %/% 2))$x
            colnames(expected) <- paste0("PC_", seq_len(ncol(expected)))
            expect_equal(observed, expected)
            expect_equal(ncol(observed), max(1L, count %/% 2))
        }
    })
}

test_that("standalone feature assessment uses default PCA for one-column subsets", {
    foreach::registerDoSEQ()
    testthat::local_mocked_bindings(
        prepare_clustering_functions_for_workers = function() {
            function(object, resolution, seed, ...) rep(c(1L, 2L), length.out = nrow(object))
        },
        .package = "ClustAssess"
    )
    dense <- matrix(seq_len(8), nrow = 1,
                    dimnames = list("gene_1", paste0("cell_", seq_len(8))))
    expected <- stats::prcomp(t(dense), rank. = 1)$x
    colnames(expected) <- "PC_1"
    for (expression in list(dense, methods::as(dense, "dgCMatrix"))) {
        result <- assess_feature_stability(
            expression, rownames(expression), 1, "all", 0.5,
            seed_sequence = c(1, 101), compute_umap_embeddings = FALSE, prune_value = 0
        )
        expect_equal(result$pca_list$all[["1"]], expected)
        expect_equal(dim(result$pca_list$all[["1"]]), c(8L, 1L))
    }
})

test_that("automatic one-column cache fallbacks use default PCA", {
    foreach::registerDoSEQ()
    local_mock_pipeline_graphs()
    testthat::local_mocked_bindings(
        assess_feature_stability = function(...) mock_feature_cache(..., omit_pca = TRUE),
        .package = "ClustAssess"
    )
    dense <- matrix(seq_len(40), nrow = 5,
                    dimnames = list(paste0("gene_", seq_len(5)), paste0("cell_", seq_len(8))))
    expected <- stats::prcomp(t(dense[1, , drop = FALSE]), rank. = 1)$x
    colnames(expected) <- "PC_1"
    for (expression in list(dense, methods::as(dense, "dgCMatrix"))) {
        result <- automatic_stability_assessment(
            expression, 2, 2, 0.5, list(all = rownames(expression)), list(all = c(1L, 5L)),
            n_top_configs = 1, compute_umap_embeddings = FALSE, save_temp = FALSE, verbose = FALSE
        )
        config <- result$all[["1"]]
        expect_equal(config$pca, expected)
        expect_identical(config$stable_config$n_features, 1L)
        expect_equal(config$stable_config$n_pcs, 1L)
    }
})

test_that("feature counts handle unbounded values before integer conversion", {
    oversized <- c(Inf, -Inf, .Machine$integer.max + 1, -.Machine$integer.max - 1)
    expect_no_warning(observed <- normalize_feature_steps(c(2, oversized), 5L))
    expect_identical(observed, c(2L, 5L))
    expect_identical(normalize_feature_steps(c(NA_real_, NaN), 5L), integer())
    expect_identical(normalize_feature_steps(c(NA_real_, NaN, 1.9, 2.9), 5L), c(1L, 2L))
    expect_identical(normalize_feature_steps(c(-0.5, 0, 0.5, 0.999, 1, 1.001), 5L), c(1L, 5L))

    local_mock_pipeline_graphs()
    forwarded_steps <- NULL
    testthat::local_mocked_bindings(
        assess_feature_stability = function(..., steps) {
            forwarded_steps <<- steps
            mock_feature_cache(..., steps = steps)
        },
        .package = "ClustAssess"
    )
    expression <- matrix(seq_len(20), nrow = 5,
                         dimnames = list(paste0("gene_", seq_len(5)), paste0("cell_", seq_len(4))))
    expect_no_warning(automatic_stability_assessment(
        expression, 2, 2, 0.5, list(all = rownames(expression)),
        list(all = c(1, 2, oversized, NA_real_, NaN)), n_top_configs = 1,
        compute_umap_embeddings = FALSE, save_temp = FALSE, verbose = FALSE,
        matrix_processing = function(x) stop("Unexpected cache miss")
    ))
    expect_identical(forwarded_steps, c(1L, 2L, 5L))
})

test_that("combined assessment keeps named summaries and PCA graph ownership", {
    local_mock_pipeline_graphs()
    testthat::local_mocked_bindings(
        assess_feature_stability = mock_feature_cache,
        assess_nn_stability = function(embedding, graph_reduction_type, ...) {
            config <- paste0(graph_reduction_type, "_nn")
            result <- list(
                n_neigh_k_corresp = setNames(list(list(`2` = 2)), config),
                n_neigh_ec_consistency = setNames(list(list(`2` = rep(
                    if (graph_reduction_type == "PCA") 1 else 0.5, nrow(embedding)))), config),
                n_different_partitions = setNames(list(c(`2` = 1)), config)
            )
            if (graph_reduction_type == "PCA") result$nn_idx <- RANN::nn2(embedding, k = 2)$nn.idx
            result
        },
        .package = "ClustAssess"
    )
    expression <- matrix(0, nrow = 5, ncol = 4,
                         dimnames = list(paste0("gene_", seq_len(5)), paste0("cell_", seq_len(4))))
    result <- automatic_stability_assessment(
        expression, 2, 2, 0.5, list(all = rownames(expression)), list(all = c(3, 5)),
        include_umap_nn_assessment = TRUE, n_top_configs = 1,
        save_temp = FALSE, verbose = FALSE
    )
    config <- result$all[["3"]]
    expect_identical(config$stable_config$base_embedding, "pca")
    expect_named(config$nn_stability, c("n_neigh_k_corresp", "n_neigh_ec_consistency",
                                       "n_different_partitions", "nn_idx"))
    for (field in names(config$nn_stability)[seq_len(3)]) {
        expect_named(config$nn_stability[[field]], c("PCA_nn", "UMAP_nn"))
    }
    expect_identical(config$nn_stability$nn_idx, RANN::nn2(config$pca, k = 2)$nn.idx)
    expected_graph <- getNNmatrix(config$nn_stability$nn_idx, 2, 0, -1)$nn
    dimnames(expected_graph) <- list(colnames(expression), colnames(expression))
    expect_equal(config$adj_matrix, expected_graph)
})

test_that("combined winners own their neighbors, UMAP seed, and output options", {
    local_mock_pipeline_graphs()
    winner <- "UMAP_nn"
    graph_arguments <- list()
    replay_arguments <- NULL
    testthat::local_mocked_bindings(
        assess_feature_stability = mock_feature_cache,
        assess_nn_stability = function(embedding, graph_reduction_type,
                                        clustering_arguments, return_partition_seeds = FALSE, ...) {
            graph_arguments[[graph_reduction_type]] <<- clustering_arguments
            configs <- paste0(graph_reduction_type, c("_snn", "_nn"))
            scores <- setNames(lapply(configs, function(config) {
                setNames(lapply(c("2", "3"), function(k) {
                    rep(if (config == winner && k == "3") 1 else 0.25, nrow(embedding))
                }), c("2", "3"))
            }), configs)
            result <- list(
                n_neigh_k_corresp = lapply(scores, function(x) lapply(x, function(y) 2)),
                n_neigh_ec_consistency = scores,
                n_different_partitions = setNames(rep(list(c(`2` = 1, `3` = 1)), 2), configs)
            )
            if (graph_reduction_type == "PCA") {
                result$nn_idx <- RANN::nn2(embedding, k = 3)$nn.idx
            } else {
                expect_true(return_partition_seeds)
                result$partition_seeds <- list(
                    UMAP_snn = list(`2` = 11L, `3` = 211L),
                    UMAP_nn = list(`2` = 11L, `3` = 111L)
                )
                result$partition_rng_kinds <- lapply(result$partition_seeds,
                    function(x) lapply(x, function(seed) RNGkind()))
            }
            result
        },
        .package = "ClustAssess"
    )
    testthat::local_mocked_bindings(
        umap = function(X, ...) {
            replay_arguments <<- list(...)
            matrix(stats::runif(nrow(X) * 2), ncol = 2)
        },
        .package = "uwot"
    )
    expression <- matrix(0, nrow = 5, ncol = 8,
                         dimnames = list(paste0("gene_", seq_len(5)), paste0("cell_", seq_len(8))))
    for (winner in c("PCA_nn", "PCA_snn", "UMAP_nn", "UMAP_snn")) {
        for (prune in c(0.2, -1)) {
            for (keep_umap in c(TRUE, FALSE)) {
                for (keep_idx in c(TRUE, FALSE)) {
                    replay_arguments <- NULL
                    result <- automatic_stability_assessment(
                        expression, 3, c(2, 3), 0.5, list(all = rownames(expression)),
                        list(all = c(3, 5)), seed_sequence = c(11, 111, 211),
                        include_umap_nn_assessment = TRUE, n_top_configs = 1,
                        ranking_criterion = "median",
                        compute_umap_embeddings = keep_umap, return_nn_idx = keep_idx,
                        prune_value = prune, clustering_arguments = list(n.start = 7, n.iter = 3),
                        umap_arguments = list(n_neighbors = 4, min_dist = 0.2, n_threads = 8),
                        save_temp = FALSE, verbose = FALSE
                    )
                    config <- result$all[["3"]]
                    expected_seed <- if (winner == "UMAP_nn") 111L else 211L
                    expected_embedding <- config$pca
                    if (startsWith(winner, "UMAP")) {
                        set.seed(expected_seed)
                        expected_embedding <- matrix(stats::runif(16), ncol = 2,
                            dimnames = list(colnames(expression), c("UMAP_1", "UMAP_2")))
                        expect_identical(config$stable_config$umap_seed, expected_seed)
                        expect_identical(replay_arguments, list(n_neighbors = 4, min_dist = 0.2,
                                                               n_threads = 1, n_sgd_threads = 1))
                        expect_false(identical(RANN::nn2(expected_embedding, k = 3)$nn.idx,
                                               RANN::nn2(config$pca, k = 3)$nn.idx))
                    } else {
                        expect_null(config$stable_config$umap_seed)
                        expect_null(replay_arguments)
                    }
                    expect_equal(config$umap, if (keep_umap) expected_embedding else NULL)
                    expected_idx <- RANN::nn2(expected_embedding, k = 3)$nn.idx
                    expect_equal(config$nn_stability$nn_idx, if (keep_idx) expected_idx else NULL)
                    expect_false("partition_seeds" %in% names(config$nn_stability))
                    expect_false("partition_rng_kinds" %in% names(config$nn_stability))
                    expected_graph <- getNNmatrix(expected_idx, 3, 0, -1)$nn
                    expected_prune <- 0
                    if (endsWith(winner, "snn")) {
                        if (prune >= 0) {
                            expected_graph <- computeSNN(expected_graph, 3, prune)
                            expected_prune <- prune
                        } else {
                            pruned <- get_highest_prune_param(expected_graph, 3)
                            expected_graph <- pruned$adj_matrix
                            expected_prune <- pruned$prune_value
                        }
                    }
                    dimnames(expected_graph) <- list(colnames(expression), colnames(expression))
                    expect_equal(config$adj_matrix, expected_graph)
                    expect_equal(config$stable_config$n_neighbours, 3)
                    expect_equal(config$stable_config$prune_param, expected_prune)
                    expect_identical(config$stable_config$base_embedding,
                                     if (startsWith(winner, "UMAP")) "umap" else "pca")
                    for (args in graph_arguments) {
                        expect_identical(args, list(n.start = 7, n.iter = 3))
                    }
                }
            }
        }
    }
})

test_that("UMAP representative seeds follow merged partitions per neighborhood", {
    testthat::local_mocked_bindings(
        prepare_clustering_functions_for_workers = function() {
            function(object, resolution, seed, ...) {
                a <- rep(c(1L, 2L), each = 4)
                b <- rep(c(1L, 2L), 4)
                k <- sum(object) / nrow(object)
                if (k == 2) {
                    if (seed == 11) a else b
                } else {
                    if (seed == 211) b else a
                }
            }
        },
        .package = "ClustAssess"
    )
    testthat::local_mocked_bindings(
        umap = function(X, ...) cbind(seq_len(nrow(X)), seq_len(nrow(X))^2),
        .package = "uwot"
    )
    embedding <- matrix(seq_len(24), ncol = 3,
                        dimnames = list(paste0("cell_", seq_len(8)), paste0("PC_", seq_len(3))))
    foreach::registerDoSEQ()
    run <- function(keep) assess_nn_stability(
        embedding, c(2, 3), seed_sequence = c(11L, 111L, 211L),
        graph_reduction_type = "UMAP", graph_type = 0, return_partition_seeds = keep
    )
    ordinary <- run(FALSE)
    with_seeds <- run(TRUE)
    expect_identical(with_seeds$partition_seeds, list(UMAP_nn = list(`2` = 111L, `3` = 11L)))
    expect_identical(with_seeds$partition_rng_kinds, list(UMAP_nn = list(`2` = RNGkind(), `3` = RNGkind())))
    with_seeds$partition_seeds <- NULL
    with_seeds$partition_rng_kinds <- NULL
    expect_equal(with_seeds, ordinary)
    expect_named(ordinary, c("n_neigh_k_corresp", "n_neigh_ec_consistency", "n_different_partitions"))
    expect_identical(tail(names(formals(assess_nn_stability)), 2),
                     c("return_nn_idx", "return_partition_seeds"))
    expect_error(assess_nn_stability(embedding, 2, return_partition_seeds = NA),
                 "return_partition_seeds parameter should be a single boolean value", fixed = TRUE)
})

test_that("UMAP merging preserves input seed order with and without provenance", {
    foreach::registerDoSEQ()
    seeds <- c(901L, 7L, 103L, 5L)
    a <- rep(c(1L, 1L, 2L, 2L), 2)
    b <- rep(c(1L, 2L, 1L, 2L), 2)
    third <- rep(c(1L, 2L, 2L, 1L), 2)
    memberships <- list(`901` = a, `7` = b, `103` = third, `5` = a)
    collection_order <- logical()
    original_foreach <- foreach::foreach
    testthat::local_mocked_bindings(
        foreach = function(..., .inorder = TRUE) {
            collection_order <<- c(collection_order, .inorder)
            call <- match.call()
            call[[1]] <- original_foreach
            eval(call, envir = parent.frame())
        },
        .package = "foreach"
    )
    testthat::local_mocked_bindings(
        prepare_clustering_functions_for_workers = function() {
            function(object, resolution, seed, ...) memberships[[as.character(seed)]]
        },
        .package = "ClustAssess"
    )
    testthat::local_mocked_bindings(
        umap = function(X, ...) cbind(seq_len(nrow(X)), seq_len(nrow(X))^2),
        .package = "uwot"
    )
    embedding <- matrix(seq_len(24), ncol = 3,
                        dimnames = list(paste0("cell_", seq_len(8)), paste0("PC_", seq_len(3))))
    # Distinct patterns have tied pairwise agreement, exercising threshold ties.
    tie_score <- element_sim(a, b)
    expect_equal(element_sim(a, third), tie_score)
    expect_equal(element_sim(b, third), tie_score)
    for (threshold in c(1, tie_score - 1e-6)) {
        expected <- merge_partitions(lapply(seeds, function(seed) {
            list(mb = memberships[[as.character(seed)]], freq = 1, seed = seed)
        }), ecs_thresh = threshold)
        outputs <- lapply(c(FALSE, TRUE), function(keep) {
            assess_nn_stability(embedding, 2, seed_sequence = seeds,
                graph_reduction_type = "UMAP", graph_type = 0,
                ecs_thresh = threshold, return_partition_seeds = keep)
        })
        expect_identical(outputs[[2]]$partition_seeds$UMAP_nn[["2"]],
                         expected$partitions[[1]]$seed)
        if (threshold == 1) expect_identical(expected$partitions[[1]]$seed, 901L)
        expect_equal(outputs[[2]]$n_neigh_ec_consistency$UMAP_nn[["2"]], expected$ecc)
        outputs[[2]]$partition_seeds <- NULL
        outputs[[2]]$partition_rng_kinds <- NULL
        expect_equal(outputs[[2]], outputs[[1]])
    }
    expect_identical(collection_order, rep(TRUE, 4))
})

test_that("tiny automatic pipeline supports combined real graph assessment", {
    testthat::skip_if_not_installed("Seurat")
    set.seed(17)
    expression <- matrix(stats::rnorm(10 * 12), nrow = 10,
                         dimnames = list(paste0("gene_", seq_len(10)), paste0("cell_", seq_len(12))))
    result <- suppressWarnings(automatic_stability_assessment(
        expression, 2, c(2, 3), 0.5, list(all = rownames(expression)), list(all = c(6, 10)),
        seed_sequence = c(11, 111), include_umap_nn_assessment = TRUE,
        n_top_configs = 1, matrix_processing = make_pca_matrix_processing(actual_npcs = 3),
        umap_arguments = list(n_neighbors = 4, n_epochs = 20, init = "random"),
        prune_value = 0, algorithms_clustering_assessment = 1,
        save_temp = FALSE, verbose = FALSE
    ))
    config <- result$all[[1]]
    expect_named(config$nn_stability$n_neigh_ec_consistency,
                 c("PCA_snn", "PCA_nn", "UMAP_snn", "UMAP_nn"))
    expected_embedding <- if (config$stable_config$base_embedding == "umap") config$umap else config$pca
    expect_identical(config$nn_stability$nn_idx, RANN::nn2(expected_embedding, k = 3)$nn.idx)
    expect_equal(dim(config$adj_matrix), c(12L, 12L))
    expect_named(config$clustering_stability, c("split_by_resolution", "split_by_k"))
})

test_that("installed PSOCK pipeline replays a winning evaluated UMAP seed", {
    testthat::skip_if_not_installed("ClustAssess")
    testthat::skip_if_not_installed("doParallel")
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if(!file.exists(file.path(system.file(package = "ClustAssess"), "Meta", "package.rds")),
                      "ClustAssess is loaded from source, not an installed package")
    cl <- parallel::makePSOCKcluster(2)
    on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
    }, add = TRUE)
    parallel::clusterCall(cl, function(paths) {
        .libPaths(paths)
        NULL
    }, .libPaths())
    parallel::clusterSetRNGStream(cl, 321)
    doParallel::registerDoParallel(cl)
    assess <- assess_nn_stability
    testthat::local_mocked_bindings(
        assess_nn_stability = function(..., graph_reduction_type) {
            result <- assess(..., graph_reduction_type = graph_reduction_type)
            # Make a real UMAP assessment win without changing its worker results.
            if (graph_reduction_type == "PCA") {
                result$n_neigh_ec_consistency <- lapply(result$n_neigh_ec_consistency,
                    function(config) lapply(config, function(x) rep(-1, length(x))))
            }
            result
        },
        .package = "ClustAssess"
    )
    set.seed(17)
    expression <- matrix(stats::rnorm(10 * 12), nrow = 10,
                         dimnames = list(paste0("gene_", seq_len(10)), paste0("cell_", seq_len(12))))
    umap_args <- list(n_neighbors = 4, n_epochs = 20, init = "random")
    result <- suppressWarnings(automatic_stability_assessment(
        expression, 2, c(2, 3), 0.5, list(all = rownames(expression)), list(all = c(6, 10)),
        seed_sequence = c(11, 111), include_umap_nn_assessment = TRUE,
        n_top_configs = 1, ranking_criterion = "median",
        matrix_processing = make_pca_matrix_processing(actual_npcs = 3),
        umap_arguments = umap_args, prune_value = 0, algorithms_clustering_assessment = 1,
        save_temp = FALSE, verbose = FALSE
    ))
    config <- result$all[[1]]
    expect_identical(config$stable_config$base_embedding, "umap")
    expect_true(config$stable_config$umap_seed %in% c(11L, 111L))
    old_kind <- RNGkind()
    on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(config$stable_config$umap_seed)
    expected_umap <- do.call(uwot::umap, c(list(X = config$pca), process_umap_arguments(umap_args, 12)))
    dimnames(expected_umap) <- dimnames(config$umap)
    expect_equal(config$umap, expected_umap, tolerance = 0)
    expect_identical(config$nn_stability$nn_idx, RANN::nn2(expected_umap, k = 3)$nn.idx)
    expect_false("partition_seeds" %in% names(config$nn_stability))
    expect_false("partition_rng_kinds" %in% names(config$nn_stability))
})

test_that("representative UMAP replay restores main RNG state and kind", {
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- .Random.seed
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    testthat::local_mocked_bindings(
        umap = function(X, ...) matrix(stats::runif(nrow(X) * 2), ncol = 2),
        .package = "uwot"
    )
    input <- matrix(seq_len(24), ncol = 3,
                    dimnames = list(paste0("cell_", seq_len(8)), paste0("PC_", seq_len(3))))
    RNGkind("L'Ecuyer-CMRG")
    worker_kind <- RNGkind()
    set.seed(111)
    expected <- matrix(stats::runif(16), ncol = 2,
                        dimnames = list(rownames(input), c("UMAP_1", "UMAP_2")))
    RNGkind("Mersenne-Twister")
    set.seed(17)
    before_kind <- RNGkind()
    before_seed <- .Random.seed
    observed <- replay_umap_embedding(input, list(n_neighbors = 3), 111L, worker_kind)
    expect_identical(observed, expected)
    expect_identical(RNGkind(), before_kind)
    expect_identical(.Random.seed, before_seed)

    rm(".Random.seed", envir = .GlobalEnv)
    expect_identical(replay_umap_embedding(input, list(), 111L, worker_kind), expected)
    expect_identical(RNGkind(), before_kind)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

    set.seed(17)
    testthat::local_mocked_bindings(umap = function(...) stop("UMAP failed"), .package = "uwot")
    expect_error(replay_umap_embedding(input, list(), 111L, worker_kind), "UMAP failed", fixed = TRUE)
    expect_identical(RNGkind(), before_kind)
    expect_identical(.Random.seed, before_seed)
})

test_that("default PCA processors do not require optional BLAS control", {
    testthat::local_mocked_bindings(
        is_package_installed = function(package_name) FALSE,
        .package = "ClustAssess"
    )
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        testthat::local_mocked_bindings(
            blas_set_num_threads = function(...) stop("Optional BLAS control was used"),
            .package = "RhpcBLASctl"
        )
    }
    input <- matrix(seq_len(40), nrow = 8,
                    dimnames = list(paste0("cell_", seq_len(8)), paste0("gene_", seq_len(5))))
    expected <- stats::prcomp(input, rank. = 2)$x
    colnames(expected) <- paste0("PC_", seq_len(ncol(expected)))
    for (pipeline in list(assess_feature_stability, automatic_stability_assessment)) {
        processor <- eval(formals(pipeline)$matrix_processing, envir = environment(pipeline))
        expect_equal(processor(input), expected)
    }
})

test_that("parallel nearest neighbors use ordinary matrices without SharedObject", {
    testthat::skip_if_not_installed("doParallel")
    testthat::local_mocked_bindings(
        is_package_installed = function(package_name) FALSE,
        .package = "ClustAssess"
    )
    if (requireNamespace("SharedObject", quietly = TRUE)) {
        testthat::local_mocked_bindings(
            share = function(...) stop("Optional sharing was used"),
            .package = "SharedObject"
        )
    }
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
    input <- cbind(c(0, 1, 4, 10, 15, 22, 33), c(2, 3, 2, 1, 4, 2, 3))
    expected <- RANN::nn2(input, k = 3)$nn.idx
    expect_identical(parallel_nn2_idx(input, 3), expected)
    expect_equal(dim(parallel_nn2_idx(input, 1)), c(7L, 1L))
})

test_that("parallel nearest neighbors retain the available SharedObject path", {
    testthat::skip_if_not_installed("SharedObject")
    testthat::skip_if_not_installed("doParallel")
    calls <- 0L
    share <- SharedObject::share
    testthat::local_mocked_bindings(
        share = function(...) {
            calls <<- calls + 1L
            share(...)
        },
        .package = "SharedObject"
    )
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
    input <- matrix(seq_len(21), ncol = 3)
    expect_identical(parallel_nn2_idx(input, 3), RANN::nn2(input, k = 3)$nn.idx)
    expect_gte(calls, 1L)
})
