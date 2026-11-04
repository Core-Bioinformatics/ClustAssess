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

test_that("a single partition with missing labels is returned without the label check", {
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

test_that("one partition found by every repetition has agreement 1", {
    x <- c(1, 1, 2, 2)
    relabelled <- c(2, 2, 1, 1)
    records <- list(
        list(mb = x, freq = 2, avg_agreement = 99),
        list(mb = x, freq = 1, avg_agreement = -99)
    )

    for (order_logic in c("freq", "avg_agreement")) {
        merged <- merge_partitions(list(x, relabelled, x), order_logic = order_logic)
        expect_length(merged$partitions, 1)
        expect_equal(merged$partitions[[1]]$freq, 3)
        expect_equal(merged$partitions[[1]]$avg_agreement, 1)

        # A supplied value is replaced, also in pre-merged and nested inputs.
        expect_equal(
            merge_partitions(records, order_logic = order_logic)$partitions[[1]]$avg_agreement,
            1
        )
        expect_equal(
            merge_partitions(list(partitions = records), order_logic = order_logic)$partitions[[1]]$avg_agreement,
            1
        )
        nested <- merge_partitions(list(records, records), order_logic = order_logic)
        expect_equal(nested[[2]]$partitions[[1]]$avg_agreement, 1)

        # A single repetition has no other repetition to agree with.
        single <- merge_partitions(
            list(list(mb = x, freq = 1, avg_agreement = 99)),
            order_logic = order_logic
        )
        expect_null(single$partitions[[1]]$avg_agreement)
        expect_null(merge_partitions(list(x), order_logic = order_logic)$partitions[[1]]$avg_agreement)
    }

    # Without ordering, no agreement is computed.
    expect_null(merge_partitions(list(x, x), order_logic = "none")$partitions[[1]]$avg_agreement)
})

test_that("the method shapes keep those of 1.2.0 and cover every method", {
    expect_equal(ClustAssess:::method_shapes(4), 21:24)
    expect_equal(ClustAssess:::method_shapes(5), 21:25)
    expect_false(anyNA(ClustAssess:::method_shapes(length(ClustAssess:::algorithm_names))))
})

six_configurations <- function() {
    result <- suppressWarnings(assess_clustering_stability(
        graph_adjacency_matrix = two_group_adjacency(),
        resolution = c(0.5, 1),
        seed_sequence = c(1, 101),
        clustering_algorithm = 1,
        verbose = FALSE
    ))
    names <- paste0("configuration_", 1:6)
    list(
        split_by_resolution = stats::setNames(rep(result$split_by_resolution, 6), names),
        split_by_k = stats::setNames(rep(result$split_by_k, 6), names)
    )
}

test_that("the k plots give each of six configurations a shape", {
    object <- six_configurations()
    for (plot in list(plot_k_resolution_corresp(object), plot_k_n_partitions(object))) {
        built <- ggplot2::ggplot_build(plot)
        shapes <- unlist(lapply(built$data, function(layer) layer$shape))
        expect_false(anyNA(shapes))
        expect_setequal(unique(shapes), 21:25)
    }
})

test_that("the app's ggplot k plots give each clustering method a shape", {
    methods <- ClustAssess:::algorithm_names
    summary <- list(first_freq = 8, total_freq = 10, ecc = c(0.7, 0.8))
    object <- list(
        split_by_resolution = stats::setNames(
            lapply(methods, function(method) list("0.5" = list(clusters = list("2" = summary)))),
            methods
        ),
        split_by_k = stats::setNames(
            lapply(methods, function(method) list("2" = c(summary, n_partitions = 2))),
            methods
        )
    )
    for (plot in list(
        ClustAssess:::shiny_ggplot_k_resolution_corresp(object),
        ClustAssess:::shiny_ggplot_k_n_partitions(object)
    )) {
        built <- ggplot2::ggplot_build(plot)
        shapes <- unlist(lapply(built$data, function(layer) layer$shape))
        expect_false(anyNA(shapes))
        expect_setequal(unique(shapes), 21:25)
    }
})

test_that("the app's k plots draw every clustering method", {
    methods <- ClustAssess:::algorithm_names
    summary_df <- expand.grid(k = 2:3, cl_method = methods, stringsAsFactors = FALSE)
    summary_df$res <- rep(c(0.5, 1), length(methods))
    summary_df$first_freq <- 8
    summary_df$total_freq <- 10
    summary_df$ecc <- seq(0.6, 0.9, length.out = nrow(summary_df))
    summary_df$n_partitions <- 2

    grDevices::pdf(NULL, width = 20, height = 20)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_error(ClustAssess:::shiny_plot_k_res(summary_df))
    expect_no_error(ClustAssess:::shiny_plot_k_n_partitions(summary_df, plt_height = 600))
})
