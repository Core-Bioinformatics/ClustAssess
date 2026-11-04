test_that("frequency ordering records the average agreement of every partition", {
    x <- c(1, 1, 2, 2)
    y <- c(1, 2, 1, 2)
    z <- c(1, 1, 1, 2)
    # x is the most frequent partition: no tie at the top
    partitions <- list(x, x, x, y, z)

    for (check_ties in c(TRUE, FALSE)) {
        merged <- merge_partitions(partitions, order_logic = "freq", check_ties = check_ties)
        agreement <- vapply(merged$partitions, `[[`, numeric(1), "avg_agreement")

        expect_length(agreement, 3)
        expect_equal(merged$partitions[[1]]$freq, 3)
        # x agrees fully with its two copies and partly with y and z
        expect_equal(agreement[[1]], (2 + element_sim(x, y) + element_sim(x, z)) / 4)
    }
})

test_that("frequency ordering recomputes supplied average agreement values", {
    x <- c(1, 1, 2, 2)
    y <- c(1, 2, 1, 2)

    merged <- merge_partitions(
        list(
            list(mb = x, freq = 2, avg_agreement = 99),
            list(mb = y, freq = 1, avg_agreement = -99)
        ),
        order_logic = "freq",
        check_ties = FALSE
    )

    expect_equal(
        vapply(merged$partitions, `[[`, numeric(1), "avg_agreement"),
        c((1 + element_sim(x, y)) / 2, element_sim(x, y))
    )
})

test_that("a single partition with missing labels is returned as it is", {
    missing_inputs <- list(
        c(1, NA_real_, 2, 2), c(1, NaN, 2, 2), c(1L, NA_integer_, 2L, 2L),
        c("a", NA_character_, "b", "b"), factor(c("a", NA, "b", "b"))
    )
    for (mb in missing_inputs) {
        for (ecs_thresh in c(1, 0.99)) {
            merged <- merge_partitions(list(mb), ecs_thresh = ecs_thresh)
            expect_length(merged$partitions, 1)
            expect_identical(merged$partitions[[1]]$mb, mb)
            expect_equal(merged$partitions[[1]]$freq, 1)
            expect_equal(merged$ecc, rep(1, 4))
        }
        merged <- merge_partitions(list(list(mb = mb, freq = 2)))
        expect_equal(merged$partitions[[1]]$freq, 2)
        expect_equal(merged$ecc, rep(1, 4))

        expect_equal(element_sim_matrix(list(mb)), matrix(1, nrow = 1, ncol = 1))

        # Comparing partitions with missing labels still stops.
        expect_error(merge_partitions(list(mb, mb)), "missing values")
        expect_error(element_sim_matrix(list(mb, c(1, 1, 2, 2))), "missing values")
    }
})

two_group_adjacency <- function() {
    group <- rep(1:2, each = 10)
    adjacency <- outer(group, group, "==") * 1
    diag(adjacency) <- 0
    adjacency[10, 11] <- adjacency[11, 10] <- 0.1
    names <- paste0("cell_", seq_along(group))
    dimnames(adjacency) <- list(names, names)
    methods::as(adjacency, "dgCMatrix")
}

test_that("algorithm 6 (leidenbase) has a name", {
    expect_identical(ClustAssess:::algorithm_names[[6]], "Leidenbase")
    expect_identical(ClustAssess:::algorithm_names_mapping[["Leidenbase"]], 6)
})

test_that("algorithm 6 (leidenbase) uses the resolution scale of Louvain", {
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if_not_installed("leidenbase")

    adjacency <- two_group_adjacency()
    for (resolution in c(0.5, 1)) {
        louvain <- ClustAssess:::clustering_functions(
            adjacency,
            resolution = resolution,
            seed = 1,
            algorithm = 1,
            verbose = FALSE
        )
        leiden <- ClustAssess:::clustering_functions(
            adjacency,
            resolution = resolution,
            seed = 1,
            algorithm = 6,
            verbose = FALSE
        )
        expect_equal(length(unique(leiden)), 2, info = paste("resolution", resolution))
        expect_equal(length(unique(leiden)), length(unique(louvain)))
    }
})

test_that("assess_clustering_stability names the results of algorithm 6", {
    testthat::skip_if_not_installed("leidenbase")

    result <- suppressWarnings(assess_clustering_stability(
        graph_adjacency_matrix = two_group_adjacency(),
        resolution = c(0.5, 1),
        seed_sequence = c(1, 101),
        clustering_algorithm = 6,
        verbose = FALSE
    ))

    expect_named(result$split_by_resolution, "Leidenbase")
    expect_named(result$split_by_k, "Leidenbase")
})

test_that("feature and graph stability accept algorithm 6 and reject 7", {
    embedding <- matrix(stats::rnorm(40), nrow = 20)
    rownames(embedding) <- paste0("cell_", seq_len(20))

    expect_error(
        assess_nn_stability(embedding, n_neigh_sequence = 5, n_repetitions = 2, clustering_algorithm = 7),
        "between 1 and 6"
    )
    data_matrix <- t(embedding)
    rownames(data_matrix) <- c("gene_1", "gene_2")
    expect_error(
        assess_feature_stability(
            data_matrix = data_matrix,
            feature_set = rownames(data_matrix),
            steps = 2,
            feature_type = "test",
            resolution = 0.5,
            n_repetitions = 2,
            clustering_algorithm = 7
        ),
        "between 1 and 6"
    )
})
