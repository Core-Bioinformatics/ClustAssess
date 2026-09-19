make_smoke_adjacency <- function() {
    adjacency <- matrix(0, nrow = 8, ncol = 8)
    adjacency[1:4, 1:4] <- 1
    adjacency[5:8, 5:8] <- 1
    diag(adjacency) <- 0
    adjacency <- methods::as(Matrix::Matrix(adjacency, sparse = TRUE), "dgCMatrix")
    cell_names <- paste0("cell_", seq_len(nrow(adjacency)))
    rownames(adjacency) <- cell_names
    colnames(adjacency) <- cell_names
    adjacency
}

reference_louvain_clustering_stability <- function(adjacency, resolution, seed_sequence) {
    result_object <- list(Louvain = list())

    for (res in resolution) {
        different_partitions <- list()

        for (i in seq_along(seed_sequence)) {
            cluster_result <- Seurat::FindClusters(
                adjacency,
                resolution = res,
                random.seed = seed_sequence[i],
                algorithm = 1,
                n.start = 1,
                n.iter = 10,
                verbose = FALSE
            )
            membership <- as.integer(cluster_result[[colnames(cluster_result)[1]]])
            k <- as.character(length(unique(membership)))

            different_partitions[[k]][[length(different_partitions[[k]]) + 1L]] <- list(
                mb = membership,
                freq = 1,
                seed = seed_sequence[i]
            )
        }

        different_partitions <- different_partitions[
            stringr::str_sort(names(different_partitions), numeric = TRUE)
        ]
        unique_partitions <- list()

        for (k in names(different_partitions)) {
            different_partitions[[k]] <- merge_partitions(
                partition_list = different_partitions[[k]],
                ecs_thresh = 1,
                order_logic = "freq"
            )
            unique_partitions <- c(unique_partitions, different_partitions[[k]]$partitions)
        }

        ecc <- if (length(unique_partitions) == 1) {
            rep(1, length(unique_partitions[[1]]$mb))
        } else {
            weighted_element_consistency(
                lapply(unique_partitions, function(x) x$mb),
                sapply(unique_partitions, function(x) x$freq)
            )
        }

        result_object$Louvain[[as.character(res)]] <- list(
            clusters = different_partitions,
            ecc = ecc
        )
    }

    list(
        split_by_resolution = result_object,
        split_by_k = list(
            Louvain = getFromNamespace("merge_resolutions", "ClustAssess")(
                result_object$Louvain
            )
        )
    )
}

test_that("clustering membership summaries reuse exact single-k merge output", {
    summarize <- getFromNamespace("summarize_clustering_memberships", "ClustAssess")
    seed_sequence <- c(1, 101, 201)
    memberships <- list(
        setNames(c(1, 1, 1, 2, 2, 2), paste0("cell_", 1:6)),
        setNames(c(1, 1, 2, 2, 2, 1), paste0("cell_", 1:6)),
        setNames(c("A", "A", "A", "B", "B", "B"), paste0("cell_", 1:6))
    )

    result <- summarize(
        clustering_memberships = memberships,
        seed_sequence = seed_sequence,
        ecs_thresh = 1
    )

    expected_merge <- merge_partitions(
        partition_list = Map(function(mb, seed) {
            list(mb = mb, freq = 1, seed = seed)
        }, memberships, seed_sequence),
        ecs_thresh = 1,
        order_logic = "freq"
    )

    expect_equal(result, list(
        clusters = list(`2` = expected_merge),
        ecc = expected_merge$ecc
    ))
})

test_that("assess_clustering_stability returns deterministic object shape", {
    testthat::skip_if_not_installed("Matrix")
    testthat::skip_if_not_installed("Seurat")

    adjacency <- make_smoke_adjacency()

    result <- suppressWarnings(
        assess_clustering_stability(
            graph_adjacency_matrix = adjacency,
            resolution = 0.5,
            seed_sequence = c(1, 101),
            clustering_algorithm = 1,
            verbose = FALSE
        )
    )

    expect_named(result, c("split_by_resolution", "split_by_k"))
    expect_named(result$split_by_resolution, "Louvain")
    expect_named(result$split_by_resolution$Louvain, "0.5")
    expect_true(is.list(result$split_by_resolution$Louvain[["0.5"]]$clusters))
    expect_true(is.numeric(result$split_by_resolution$Louvain[["0.5"]]$ecc))
    expect_true(is.list(result$split_by_k$Louvain[["2"]]$partitions))
    expect_true(is.numeric(result$split_by_k$Louvain[["2"]]$ecc))
})

test_that("assess_clustering_stability matches resolution-seed reference", {
    testthat::skip_if_not_installed("Matrix")
    testthat::skip_if_not_installed("Seurat")

    adjacency <- make_smoke_adjacency()
    resolution <- c(1, 0.5)
    seed_sequence <- c(1, 101, 201)

    result <- suppressWarnings(
        assess_clustering_stability(
            graph_adjacency_matrix = adjacency,
            resolution = resolution,
            seed_sequence = seed_sequence,
            clustering_algorithm = 1,
            verbose = FALSE
        )
    )
    expected <- reference_louvain_clustering_stability(
        adjacency = adjacency,
        resolution = resolution,
        seed_sequence = seed_sequence
    )

    expect_named(result$split_by_resolution$Louvain, as.character(resolution))
    expect_equal(result, expected)
})

test_that("clustering worker dependencies follow the selected algorithm", {
    testthat::skip_if_not_installed("Matrix")

    # Inspect the dispatch contract without requiring a Python Leiden runtime.
    requested_packages <- NULL
    original_foreach <- foreach::foreach
    testthat::local_mocked_bindings(
        foreach = function(..., .packages = character()) {
            requested_packages <<- union(requested_packages, .packages)
            iterator <- original_foreach(..., .packages = character())
            # foreach must still evaluate loop expressions in its caller.
            iterator$evalenv <- parent.frame()
            iterator
        },
        .package = "foreach"
    )
    testthat::local_mocked_bindings(
        prepare_clustering_functions_for_workers = function() {
            function(object, ...) rep(1:2, each = nrow(object) / 2)
        },
        .package = "ClustAssess"
    )

    for (algorithm in c(1, 4, 5, 6)) {
        requested_packages <- character()
        assess_clustering_stability(
            graph_adjacency_matrix = make_smoke_adjacency(),
            resolution = 0.5,
            seed_sequence = 1,
            clustering_algorithm = algorithm,
            clustering_arguments = list(),
            verbose = FALSE
        )
        expect_identical(
            "leiden" %in% requested_packages,
            algorithm %in% 4:5,
            info = paste("algorithm", algorithm)
        )
    }
})

test_that("assess_clustering_stability works with PSOCK workers when package-installed", {
    testthat::skip_if_not_installed("ClustAssess")
    testthat::skip_if_not_installed("doParallel")
    testthat::skip_if_not_installed("Matrix")
    testthat::skip_if_not_installed("Seurat")

    package_path <- system.file(package = "ClustAssess")
    testthat::skip_if(!nzchar(package_path), "ClustAssess is not installed in a library path")
    testthat::skip_if(
        !file.exists(file.path(package_path, "Meta", "package.rds")),
        "ClustAssess is loaded from source, not an installed package"
    )

    adjacency <- make_smoke_adjacency()

    expected <- suppressWarnings(
        ClustAssess::assess_clustering_stability(
            graph_adjacency_matrix = adjacency,
            resolution = 0.5,
            seed_sequence = c(1, 101),
            clustering_algorithm = 1,
            verbose = FALSE
        )
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

    result <- suppressWarnings(
        ClustAssess::assess_clustering_stability(
            graph_adjacency_matrix = adjacency,
            resolution = 0.5,
            seed_sequence = c(1, 101),
            clustering_algorithm = 1,
            verbose = FALSE
        )
    )

    expect_named(result$split_by_resolution, "Louvain")
    expect_named(result$split_by_resolution$Louvain, "0.5")
    expect_equal(result, expected)
})
