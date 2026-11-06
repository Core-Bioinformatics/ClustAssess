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

        # Comparing different partitions with missing labels still stops.
        expect_error(merge_partitions(list(mb, c(1, 2, 1, 2))), "missing values")
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

test_that("one partition found by every repetition has agreement 1 unless its input gives one", {
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
        # The kept partition is the first; its input has no agreement.
        expect_equal(
            merge_partitions(
                list(list(mb = x, freq = 2), list(mb = x, freq = 1, avg_agreement = -99)),
                order_logic = order_logic
            )$partitions[[1]]$avg_agreement,
            1
        )

        # An agreement given with the kept partition stays, as in 1.2.0, also
        # in pre-merged and nested inputs and for a single repetition.
        expect_equal(
            merge_partitions(records, order_logic = order_logic)$partitions[[1]]$avg_agreement,
            99
        )
        expect_equal(
            merge_partitions(list(partitions = records), order_logic = order_logic)$partitions[[1]]$avg_agreement,
            99
        )
        nested <- merge_partitions(list(records, records), order_logic = order_logic)
        expect_equal(nested[[2]]$partitions[[1]]$avg_agreement, 99)
        single <- merge_partitions(
            list(list(mb = x, freq = 1, avg_agreement = 99)),
            order_logic = order_logic
        )
        expect_equal(single$partitions[[1]]$avg_agreement, 99)

        # A single repetition has no other repetition to agree with.
        expect_null(merge_partitions(list(x), order_logic = order_logic)$partitions[[1]]$avg_agreement)
    }

    # Without ordering, no agreement is computed.
    expect_null(merge_partitions(list(x, x), order_logic = "none")$partitions[[1]]$avg_agreement)
})

test_that("ECS and ECC values keep the arithmetic of 1.2.0", {
    # 1.2.0's arithmetic: corrected_l1_mb() (R code) where either vector has
    # character labels, and the native scores of the labels themselves for
    # numbers and factors; the weighted ECC adds the pairs in this order.
    elscore_1_2_0 <- function(a, b) {
        if (is.character(a) || is.character(b)) corrected_l1_mb(a, b) else disjointECS(a, b)
    }
    ecc_1_2_0 <- function(mbs, weights) {
        total <- sum(weights)
        norm_factor <- total * (total - 1) / 2
        consistency <- 0
        for (i in seq_len(length(mbs) - 1)) {
            consistency <- consistency + weights[i] / norm_factor * (weights[i] - 1) / 2
            for (j in seq(from = i + 1, to = length(mbs))) {
                consistency <- consistency + elscore_1_2_0(mbs[[i]], mbs[[j]]) * weights[i] / norm_factor * weights[j]
            }
        }
        unname(consistency + weights[length(mbs)] / norm_factor * (weights[length(mbs)] - 1) / 2)
    }
    set.seed(7)
    for (i in 1:48) {
        n <- sample(c(11, 60, 300), 1)
        k <- sample(2:6, 1)
        labels <- switch(i %% 4 + 1,
            function() sample(k, n, TRUE),
            function() sample(k, n, TRUE) - 1L,
            function() paste0("c", sample(k, n, TRUE)),
            function() factor(sample(letters[1:k], n, TRUE), levels = letters[6:1])
        )
        mbs <- replicate(sample(2:7, 1), labels(), simplify = FALSE)
        weights <- sample(1:4, length(mbs), TRUE)
        a <- mbs[[1]]
        b <- mbs[[2]]
        expect_identical(unname(weighted_element_consistency(mbs, weights)), ecc_1_2_0(mbs, weights))
        expect_identical(unname(element_sim_elscore(a, b)), elscore_1_2_0(a, b))
        expect_identical(element_sim(a, b), if (is.character(a)) mean(corrected_l1_mb(a, b)) else disjointECSaverage(a, b))
        expect_identical(element_sim_matrix(list(a, b))[1, 2], mean(elscore_1_2_0(a, b)))
        expect_identical(
            element_agreement(a, mbs[-1]),
            Reduce(`+`, lapply(mbs[-1], function(x) corrected_l1_mb(a, x))) / (length(mbs) - 1)
        )
    }

    # The configuration chosen from ECC medians, as the automatic assessment
    # chooses one, is that of 1.2.0.
    set.seed(2)
    p <- replicate(5, sample.int(4, 11, TRUE), simplify = FALSE)
    m <- vapply(list(p, rev(p)), function(x) stats::median(merge_partitions(x)$ecc), numeric(1))
    expect_identical(which(m >= stats::fivenum(m)[4]), 1L)

    # A threshold equal to the 1.2.0 score of a pair merges it, and a threshold
    # equal to a higher score from another arithmetic keeps both partitions.
    find_pair <- function(seed, labels, other, direction) {
        set.seed(seed)
        for (i in 1:1000) {
            pair <- list(labels(), labels())
            if (direction(other(pair[[1]], pair[[2]]), mean(elscore_1_2_0(pair[[1]], pair[[2]])))) {
                return(pair)
            }
        }
        stop("no pair found")
    }
    survivors <- function(pair, threshold) {
        length(merge_partitions(pair, ecs_thresh = threshold, order_logic = "none")$partitions)
    }
    cases <- list(
        numeric = list(seed = 981, labels = function() sample(1:4, 31, TRUE), other = element_sim),
        character = list(seed = 52, labels = function() paste0("c", sample(1:5, 40, TRUE)),
                         other = function(a, b) disjointECSaverage(match(a, unique(a)), match(b, unique(b))))
    )
    for (case in cases) {
        below <- find_pair(case$seed, case$labels, case$other, `<`)
        expect_equal(survivors(below, mean(elscore_1_2_0(below[[1]], below[[2]]))), 1)
        above <- find_pair(case$seed, case$labels, case$other, `>`)
        expect_equal(survivors(above, case$other(above[[1]], above[[2]])), 2)
    }

    # corrected_l1_mb_at_once() gives the values of corrected_l1_mb(), with
    # labels equal only when they are equal.
    set.seed(11)
    for (i in 1:60) {
        n <- sample(c(5, 31, 200), 1)
        first <- paste0("c", sample(sample(2:9, 1), n, TRUE))
        second <- switch(i %% 3 + 1,
            paste0("d", sample(sample(2:9, 1), n, TRUE)),
            sample(sample(2:9, 1), n, TRUE) * 1.5,
            factor(sample(letters[1:5], n, TRUE), levels = letters[6:1])
        )
        expect_identical(corrected_l1_mb_at_once(first, second), corrected_l1_mb(first, second))
        expect_identical(corrected_l1_mb_at_once(second, first), corrected_l1_mb(second, first))
    }
    by_column <- c("a", "a", "b", "b")
    expect_identical(
        corrected_l1_mb_at_once(c(1, 1 + 1e-15, 1, 1 + 1e-15), by_column),
        corrected_l1_mb(c(1, 2, 1, 2), by_column)
    )
    expect_identical(
        corrected_l1_mb_at_once(c("", "x", "", "x"), by_column),
        corrected_l1_mb(c("e", "x", "e", "x"), by_column)
    )
})

test_that("feature and graph stability take the repetitions in seed order", {
    # A backend that returns unordered results in reverse order, as a backend
    # that returns results as the workers finish can.
    reverse_completion <- function() {
        foreach::setDoPar(function(obj, expr, envir, data) {
            results <- foreach:::doSEQ(obj, expr, envir, data)
            if (!obj$combineInfo$in.order && is.list(results)) rev(results) else results
        }, data = NULL, info = function(data, item) {
            switch(item, workers = 1L, name = "reverse completion", version = "1")
        })
    }
    on.exit(foreach::registerDoSEQ(), add = TRUE)
    seeds <- c(901L, 7L, 103L, 5L)
    patterns <- list(c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L), c(1L, 2L, 2L, 1L), c(1L, 1L, 1L, 2L))
    memberships <- stats::setNames(lapply(patterns, function(p) rep(p, 6)), seeds)
    seed_of <- function(x) {
        mb <- if (is.list(x)) x$mb else x
        names(memberships)[vapply(memberships, function(m) identical(as.integer(mb), m), logical(1))]
    }
    set.seed(3)
    genes <- matrix(stats::rpois(40 * 24, 5) + 0, nrow = 40,
                    dimnames = list(paste0("g", 1:40), paste0("c", 1:24)))
    embedding <- matrix(stats::rnorm(24 * 4), nrow = 24,
                        dimnames = list(paste0("c", 1:24), paste0("PC_", 1:4)))

    # The partitions reach merge_partitions() in seed order.
    received <- character()
    real_merge <- merge_partitions
    testthat::local_mocked_bindings(
        prepare_clustering_functions_for_workers = function() {
            function(object, resolution, seed, ...) memberships[[as.character(seed)]]
        },
        merge_partitions = function(partition_list, ...) {
            received <<- c(received, paste(vapply(partition_list, seed_of, ""), collapse = ","))
            real_merge(partition_list, ...)
        },
        .package = "ClustAssess"
    )
    reverse_completion()
    assess_feature_stability(genes, rownames(genes), steps = 20, feature_type = "genes",
        resolution = 0.8, n_repetitions = 4, seed_sequence = seeds, prune_value = 0.05,
        matrix_processing = function(dt_mtx, actual_npcs = 4, ...) {
            z <- stats::prcomp(dt_mtx, rank. = actual_npcs)$x
            colnames(z) <- paste0("PC_", seq_len(ncol(z)))
            z
        }, compute_umap_embeddings = FALSE)
    assess_nn_stability(embedding, c(4, 6), n_repetitions = 4, seed_sequence = seeds,
        graph_reduction_type = "PCA", graph_type = 2, prune_value = 0.05)
    expect_gt(length(received), 0)
    expect_true(all(received == paste(seeds, collapse = ",")))

    # The counts of connected components of the UMAP graphs follow the seeds.
    components <- lapply(c(FALSE, TRUE), function(reversed) {
        if (reversed) reverse_completion() else foreach::registerDoSEQ()
        get_nn_conn_comps(embedding, c(3, 5), n_repetitions = 4, seed_sequence = seeds, include_umap = TRUE,
            umap_arguments = list(n_neighbors = 5, n_epochs = 30, init = "random"))
    })
    expect_identical(components[[2]], components[[1]])
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

# The arguments of each plot() call of an app plot function.
capture_plot_calls <- function(fn, ...) {
    env <- new.env(parent = environment(fn))
    env$calls <- list()
    env$plot <- function(...) {
        env$calls[[length(env$calls) + 1L]] <- list(...)
        graphics::plot(...)
    }
    environment(fn) <- env
    fn(...)
    env$calls
}

test_that("the app's k plots give each clustering method its shape", {
    methods <- ClustAssess:::algorithm_names
    summary_df <- expand.grid(k = 2:3, cl_method = methods, stringsAsFactors = FALSE)
    summary_df$res <- rep(c(0.5, 1), length(methods))
    summary_df$first_freq <- 8
    summary_df$total_freq <- 10
    summary_df$ecc <- seq(0.6, 0.9, length.out = nrow(summary_df))
    summary_df$n_partitions <- 2
    # Two points of each of the six methods, in their order.
    expected <- rep(c(21:25, 21), each = 2)

    grDevices::pdf(NULL, width = 10, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)
    calls <- capture_plot_calls(ClustAssess:::shiny_plot_k_res, summary_df)
    expect_equal(calls[[1]]$pch, expected)
    # With its legend, at a small text size and at the default of the app.
    for (text_size in c(1, 1.5)) {
        calls <- capture_plot_calls(
            ClustAssess:::shiny_plot_k_n_partitions,
            summary_df,
            plt_height = 576,
            display_legend = TRUE,
            text_size = text_size
        )
        expect_equal(calls[[1]]$pch, expected)
        expect_equal(calls[[2]]$pch, c(21:25, 21))
    }
})

# The k summary of four clustering methods, as the app plots it.
k_summary <- function() {
    methods <- ClustAssess:::algorithm_names[c(1, 2, 3, 6)]
    summary_df <- expand.grid(k = 2:9, cl_method = methods, stringsAsFactors = FALSE)
    summary_df$first_freq <- 8
    summary_df$total_freq <- rep(c(5, 10, 20, 15), length.out = nrow(summary_df))
    summary_df$ecc <- seq(0.6, 0.95, length.out = nrow(summary_df))
    summary_df$n_partitions <- rep(1:4, length.out = nrow(summary_df))
    summary_df
}

test_that("a legend row keeps its text size while its plot has room", {
    fit <- ClustAssess:::fit_legend_size
    # A legend of one inch per unit of text size, beside a plot with 5 inches
    # for the legend and itself.
    height <- function(size) size
    region <- function(legend) 5 - legend
    expect_equal(fit(2, height, region), 2)
    smaller <- fit(10, height, region)
    expect_lt(smaller, 10)
    expect_gte(region(height(smaller)), 0.1)
    expect_null(fit(10, function(size) 5, region))
})

test_that("the app's plots with a legend row draw at every size of their sliders", {
    ns <- asNamespace("ClustAssess")
    summary_df <- k_summary()
    groups <- rep(c("Louvain", "Louvain.refined", "SLM"), times = 6)
    resolutions <- rep(seq(0.6, 1.1, by = 0.1), each = 3)
    boxes <- lapply(seq_along(groups), function(i) seq(0.5, 1, length.out = 20))
    # The app sets the colours of each number of groups when it starts.
    had_colours <- exists("discrete_colors", envir = ns$pkg_env, inherits = FALSE)
    old_colours <- if (had_colours) get("discrete_colors", envir = ns$pkg_env)
    ns$add_env_variable("discrete_colors", list("3" = c("#1B9E77", "#D95F02", "#7570B3")))
    on.exit(
        if (had_colours) {
            assign("discrete_colors", old_colours, envir = ns$pkg_env)
        } else {
            rm("discrete_colors", envir = ns$pkg_env)
        },
        add = TRUE
    )
    embedding <- cbind(seq(-1, 1, length.out = 200), sin(seq(0, 6, length.out = 200)))
    labels <- paste("cluster", 1:8)
    discrete <- factor(rep(labels, length.out = 200), levels = labels)
    continuous <- seq(0, 1, length.out = 200)

    draw <- list(
        k = function(h, w, size) {
            ns$shiny_plot_k_n_partitions(summary_df, plt_height = h * ns$ppi, display_legend = TRUE,
                text_size = size, pt_size_range = c(1.35, 2.35))
        },
        k_large_points = function(h, w, size) {
            ns$shiny_plot_k_n_partitions(summary_df, plt_height = h * ns$ppi, display_legend = TRUE,
                text_size = size, pt_size_range = c(0.1, 5))
        },
        boxplot = function(h, w, size) {
            ns$grouped_boxplot_list(boxes, groups, resolutions, plt_height = h * ns$ppi,
                plt_width = w * ns$ppi, xlab_text = "resolution", display_legend = TRUE, text_size = size)
        },
        umap = function(h, w, size) {
            ns$color_plot2(embedding, discrete, grDevices::rainbow(8), plt_height = h * ns$ppi,
                plt_width = w * ns$ppi, unique_values = labels, display_legend = TRUE, legend_text_size = size)
        },
        colour_bar = function(h, w, size) {
            ns$color_plot2(embedding, continuous, NULL, plt_height = h * ns$ppi, plt_width = w * ns$ppi,
                display_legend = TRUE, legend_text_size = size)
        },
        # A legend height given by the caller, as the sandbox gives it.
        given_height = function(h, w, size) {
            ns$color_plot2(embedding, continuous, NULL, plt_height = h * ns$ppi, plt_width = w * ns$ppi,
                display_legend = TRUE, predicted_height = size / 2)
        }
    )
    # Downloads are from 3 inches high and wide; the sliders go up to 10.
    for (name in names(draw)) {
        for (h in c(3, 7)) {
            for (w in c(3, 7)) {
                for (size in c(0.1, 0.5, 1.5, 2, 4.5, 7.5, 10)) {
                    grDevices::pdf(NULL, width = w, height = h)
                    result <- tryCatch(
                        {
                            draw[[name]](h, w, size)
                            "drawn"
                        },
                        error = conditionMessage
                    )
                    grDevices::dev.off()
                    expect_identical(result, "drawn", label = sprintf("%s, %g x %g inches, size %g", name, w, h, size))
                }
            }
        }
    }
})

test_that("the app's k plot legend is drawn smaller only when the plot has no room", {
    ns <- asNamespace("ClustAssess")
    # Each plot of the app has its own device, as here.
    k_plot_calls <- function(height, text_size) {
        grDevices::pdf(NULL, width = 10, height = height)
        on.exit(grDevices::dev.off())
        capture_plot_calls(ns$shiny_plot_k_n_partitions, k_summary(),
            plt_height = height * ns$ppi, display_legend = TRUE, text_size = text_size)
    }
    calls <- k_plot_calls(8, 1.5)
    expect_equal(calls[[1]]$cex.axis, 1.5)
    expect_equal(calls[[2]]$cex.main, 1.5)
    calls <- k_plot_calls(7, 10)
    expect_equal(calls[[1]]$cex.axis, 10)
    expect_lt(calls[[2]]$cex.main, 10)
})

test_that("the app's plots keep their legend size wherever it leaves the plot room", {
    ns <- asNamespace("ClustAssess")
    had_colours <- exists("discrete_colors", envir = ns$pkg_env, inherits = FALSE)
    old_colours <- if (had_colours) get("discrete_colors", envir = ns$pkg_env)
    ns$add_env_variable("discrete_colors", list("3" = c("#1B9E77", "#D95F02", "#7570B3")))
    on.exit(
        if (had_colours) {
            assign("discrete_colors", old_colours, envir = ns$pkg_env)
        } else {
            rm("discrete_colors", envir = ns$pkg_env)
        },
        add = TRUE
    )
    groups <- rep(c("Louvain", "Louvain.refined", "SLM"), times = 6)
    resolutions <- rep(seq(0.6, 1.1, by = 0.1), each = 3)
    boxes <- lapply(seq_along(groups), function(i) seq(0.5, 1, length.out = 20))
    embedding <- cbind(seq(-1, 1, length.out = 200), sin(seq(0, 6, length.out = 200)))
    labels <- paste("cluster", 1:8)
    discrete <- factor(rep(labels, length.out = 200), levels = labels)

    # The margins (inches) above and below a plot, at a character size.
    grDevices::pdf(NULL)
    margins <- sum(graphics::par("mar")[c(1, 3)]) * graphics::par("cin")[2]
    grDevices::dev.off()
    # The height (inches) of the legend row, the last row of the layout.
    legend_row <- function(draw, height) {
        grDevices::pdf(NULL, width = 7, height = height)
        on.exit(grDevices::dev.off())
        draw(height)
        graphics::par("fin")[2]
    }
    plots <- list(
        # A layout of three columns draws the margins at 0.66 of their size;
        # the plot is 0.1 inch shorter than the room above the legend.
        k = list(
            draw = function(height) {
                ns$shiny_plot_k_n_partitions(k_summary(), plt_height = height * ns$ppi, display_legend = TRUE,
                    text_size = 3, pt_size_range = c(1.35, 2.35))
            },
            height = function(row, room) row + 0.1 + margins * 0.66 + room
        ),
        boxplot = list(
            draw = function(height) {
                ns$grouped_boxplot_list(boxes, groups, resolutions, plt_height = height * ns$ppi,
                    plt_width = 7 * ns$ppi, xlab_text = "resolution", display_legend = TRUE, text_size = 1)
            },
            height = function(row, room) row + 0.05 + margins + room
        ),
        # The rows of the UMAP are drawn at 2.5 cm per inch.
        umap = list(
            draw = function(height) {
                ns$color_plot2(embedding, discrete, grDevices::rainbow(8), plt_height = height * ns$ppi,
                    plt_width = 7 * ns$ppi, unique_values = labels, display_legend = TRUE, legend_text_size = 1)
            },
            height = function(row, room) (row + margins + room) * 2.54 / 2.5
        )
    )
    for (name in names(plots)) {
        plot <- plots[[name]]
        # The legend row at its size, with room to spare.
        row <- legend_row(plot$draw, 12)
        # 0.05 inch of room inside the margins of the plot: the legend keeps
        # its size. Without room, the plot stopped: the legend is smaller.
        expect_equal(legend_row(plot$draw, plot$height(row, 0.05)), row, label = name)
        expect_lt(legend_row(plot$draw, plot$height(row, -0.05)), row, label = name)
    }
})

# A stability.h5 with the overall feature stability of two feature sets of
# 100 and 200 features, as the app reads it from its folder.
write_feature_stability <- function(folder) {
    file <- file.path(folder, "stability.h5")
    rhdf5::h5createFile(file)
    rhdf5::h5createGroup(file, "feature_stability")
    rhdf5::h5createGroup(file, "feature_stability/overall")
    scores <- expand.grid(ftype = c("A", "B"), fsize = c("100", "200"), rep = 1:20, stringsAsFactors = FALSE)
    scores$ecc <- seq(0.5, 1, length.out = nrow(scores))
    rhdf5::h5write(scores, file, "feature_stability/overall/by_step")
    rhdf5::h5write(scores, file, "feature_stability/overall/incremental")
}

test_that("the app's feature stability boxplots draw at every size of their sliders", {
    testthat::skip_if_not_installed("rhdf5")
    ns <- asNamespace("ClustAssess")
    had_colours <- exists("discrete_colors", envir = ns$pkg_env, inherits = FALSE)
    old_colours <- if (had_colours) get("discrete_colors", envir = ns$pkg_env)
    ns$add_env_variable("discrete_colors", list("2" = c("#1B9E77", "#D95F02")))
    folder <- tempfile("feature-stability")
    dir.create(folder)
    old_folder <- setwd(folder)
    on.exit(
        {
            setwd(old_folder)
            unlink(folder, recursive = TRUE)
            if (had_colours) {
                assign("discrete_colors", old_colours, envir = ns$pkg_env)
            } else {
                rm("discrete_colors", envir = ns$pkg_env)
            }
        },
        add = TRUE
    )
    write_feature_stability(folder)
    ordering <- list(
        original = list(A = c("100", "200"), B = c("100", "200")),
        original_incremental = list(A = c("100", "200"), B = c("100", "200"))
    )
    draw <- list(
        by_step = ns$shiny_plot_feature_stability_boxplot,
        incremental = ns$shiny_plot_feature_stability_incremental
    )
    for (name in names(draw)) {
        for (h in c(3, 7)) {
            for (w in c(3, 7)) {
                for (size in c(0.1, 0.5, 1.5, 2, 4.5, 7.5, 10)) {
                    grDevices::pdf(NULL, width = w, height = h)
                    result <- tryCatch(
                        {
                            draw[[name]]("overall", ordering, plt_height = h * ns$ppi, plt_width = w * ns$ppi,
                                text_size = size)
                            "drawn"
                        },
                        error = conditionMessage
                    )
                    grDevices::dev.off()
                    expect_identical(result, "drawn", label = sprintf("%s, %g x %g inches, size %g", name, w, h, size))
                }
            }
        }
    }
})

test_that("the app's feature stability boxplots keep their legend size wherever it leaves room", {
    testthat::skip_if_not_installed("rhdf5")
    ns <- asNamespace("ClustAssess")
    had_colours <- exists("discrete_colors", envir = ns$pkg_env, inherits = FALSE)
    old_colours <- if (had_colours) get("discrete_colors", envir = ns$pkg_env)
    ns$add_env_variable("discrete_colors", list("2" = c("#1B9E77", "#D95F02")))
    folder <- tempfile("feature-stability")
    dir.create(folder)
    old_folder <- setwd(folder)
    on.exit(
        {
            setwd(old_folder)
            unlink(folder, recursive = TRUE)
            if (had_colours) {
                assign("discrete_colors", old_colours, envir = ns$pkg_env)
            } else {
                rm("discrete_colors", envir = ns$pkg_env)
            }
        },
        add = TRUE
    )
    write_feature_stability(folder)
    ordering <- list(
        original = list(A = c("100", "200"), B = c("100", "200")),
        original_incremental = list(A = c("100", "200"), B = c("100", "200"))
    )
    draw <- list(
        by_step = ns$shiny_plot_feature_stability_boxplot,
        incremental = ns$shiny_plot_feature_stability_incremental
    )
    # The legend is drawn in the bottom margin (inches) of the boxplots.
    bottom_margin <- function(plot, height) {
        grDevices::pdf(NULL, width = 7, height = height)
        on.exit(grDevices::dev.off())
        plot("overall", ordering, plt_height = height * ns$ppi, plt_width = 7 * ns$ppi, text_size = 3)
        graphics::par("mai")[1]
    }
    grDevices::pdf(NULL)
    top_margin <- graphics::par("mai")[3]
    grDevices::dev.off()
    for (name in names(draw)) {
        # The bottom margin with the legend at its size, with room to spare.
        margin <- bottom_margin(draw[[name]], 12)
        # 0.05 inch of room inside the margins of the boxplots: the legend
        # keeps its size. Without room, the plot stopped: the legend is smaller.
        expect_equal(bottom_margin(draw[[name]], margin + top_margin + 0.05), margin, label = name)
        expect_lt(bottom_margin(draw[[name]], margin + top_margin - 0.05), margin, label = name)
    }
})

test_that("the app's separate colour bar draws at every legend size", {
    ns <- asNamespace("ClustAssess")
    values <- seq(0, 1, length.out = 500)
    for (size in c(0.1, 0.2, 0.3, 0.5, 1, 2, 4.5, 10)) {
        # The app makes the legend as high as three lines of its text.
        grDevices::pdf(NULL)
        text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = size)
        grDevices::dev.off()
        for (w in c(2, 3, 7)) {
            grDevices::pdf(NULL, width = w, height = text_height)
            result <- tryCatch(
                {
                    ns$only_legend_plot(NULL, NULL, values, plt_width = w * ns$ppi, text_size = size)
                    "drawn"
                },
                error = conditionMessage
            )
            grDevices::dev.off()
            expect_identical(result, "drawn", label = sprintf("width %g, size %g", w, size))
        }
    }
})

test_that("the app's k plot draws on a narrow screen", {
    ns <- asNamespace("ClustAssess")
    for (w in c(2, 2.4, 2.5)) {
        for (size in c(0.5, 1.5, 4.5)) {
            grDevices::pdf(NULL, width = w, height = 7)
            result <- tryCatch(
                {
                    ns$shiny_plot_k_n_partitions(k_summary(), plt_height = 7 * ns$ppi, display_legend = TRUE,
                        text_size = size, pt_size_range = c(1.35, 2.35))
                    "drawn"
                },
                error = conditionMessage
            )
            grDevices::dev.off()
            expect_identical(result, "drawn", label = sprintf("width %g, size %g", w, size))
        }
    }
})

test_that("a UMAP colour bar keeps the height given by the caller where it fits", {
    ns <- asNamespace("ClustAssess")
    embedding <- cbind(seq(-1, 1, length.out = 200), sin(seq(0, 6, length.out = 200)))
    values <- seq(0, 1, length.out = 200)
    colour_bar_height <- function(predicted_height) {
        grDevices::pdf(NULL, width = 7, height = 7)
        on.exit(grDevices::dev.off())
        ns$color_plot2(embedding, values, NULL, plt_height = 7 * ns$ppi, plt_width = 7 * ns$ppi,
            display_legend = TRUE, predicted_height = predicted_height)
        # The colour bar is drawn last, in its row of the layout.
        graphics::par("fin")[2]
    }
    # The rows are drawn at 2.5 cm per inch.
    expect_equal(colour_bar_height(0.205), 0.205 * 2.5 / 2.54)
    expect_equal(colour_bar_height(1), 1 * 2.5 / 2.54)
    # A row too low for the 0.1 inch above and below the colour bar is raised.
    expect_equal(colour_bar_height(0.1), 0.21)
    # A row that leaves the UMAP 0.05 inch inside its margins is kept; a
    # higher one, which leaves it no room, is lowered.
    grDevices::pdf(NULL)
    margins <- sum(graphics::par("mar")[c(1, 3)]) * graphics::par("cin")[2]
    grDevices::dev.off()
    tight <- 7 - (margins + 0.05) * 2.54 / 2.5
    expect_equal(colour_bar_height(tight), tight * 2.5 / 2.54)
    expect_lt(colour_bar_height(tight + 0.1), (tight + 0.1) * 2.5 / 2.54)
})

# Whether R draws a plot in every panel of a layout (FALSE: it stops).
draws_every_panel <- function(mat, heights, last_par = NULL) {
    tryCatch(
        {
            graphics::layout(mat, heights = heights)
            for (panel in seq_len(max(mat))) {
                if (panel == max(mat) && !is.null(last_par)) {
                    graphics::par(last_par)
                }
                graphics::plot.new()
            }
            TRUE
        },
        error = function(e) FALSE
    )
}

test_that("the check of a legend layout gives R's answer at the limits of its margins", {
    ns <- asNamespace("ClustAssess")
    grDevices::pdf(NULL)
    margins <- sum(graphics::par("mar")[c(1, 3)]) * graphics::par("cin")[2]
    grDevices::dev.off()
    two_rows <- matrix(c(1, 2), nrow = 2)
    k_rows <- matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = TRUE)
    # Each layout at the limit of one panel, as the app's plots make them: a
    # colour bar of 0.2 inch with 0.1 inch above and below it, a UMAP or a k
    # plot with no room inside its margins, and k legend plots in a row as
    # high as their margins (a layout of three columns draws them at 0.66).
    layouts <- function(d) {
        bar <- 0.2 * 2.54 / 2.5 + d
        umap <- 7 - (margins + d) * 2.54 / 2.5
        k_plot <- 7 - 0.1 - margins * 0.66 - d
        k_legend <- margins * 0.66 + d
        list(
            bar = list(two_rows, c(graphics::lcm((7 - bar) * 2.5), graphics::lcm(bar * 2.5)),
                list(mai = c(0.1, 0, 0.1, 0))),
            umap = list(two_rows, c(graphics::lcm((7 - umap) * 2.5), graphics::lcm(umap * 2.5)),
                list(mar = c(0, 0, 0, 0))),
            k_plot = list(k_rows, c(graphics::lcm((7 - k_plot - 0.1) * 2.54), graphics::lcm(k_plot * 2.54))),
            k_legend = list(k_rows, c(graphics::lcm((7 - k_legend - 0.1) * 2.54), graphics::lcm(k_legend * 2.54)))
        )
    }
    answers <- character()
    for (d in c(-1e-12, -4e-16, -2e-16, 0, 2e-16, 4e-16, 1e-12)) {
        cases <- layouts(d)
        for (name in names(cases)) {
            case <- cases[[name]]
            last_par <- if (length(case) > 2) case[[3]]
            grDevices::pdf(NULL, width = 7, height = 7)
            expected <- draws_every_panel(case[[1]], case[[2]], last_par)
            grDevices::dev.off()
            grDevices::pdf(NULL, width = 7, height = 7)
            checked <- ns$layout_draws(case[[1]], case[[2]], last_par)
            grDevices::dev.off()
            expect_identical(checked, expected, label = sprintf("%s at %g", name, d))
            answers <- c(answers, expected)
        }
    }
    # The limits are crossed: some layouts draw and some do not.
    expect_true(any(answers == "TRUE") && any(answers == "FALSE"))
})

test_that("a UMAP colour bar keeps the height given by the caller wherever R draws it", {
    ns <- asNamespace("ClustAssess")
    embedding <- cbind(seq(-1, 1, length.out = 200), sin(seq(0, 6, length.out = 200)))
    values <- seq(0, 1, length.out = 200)
    for (d in c(-4e-16, -2e-16, 0, 2e-16, 4e-16)) {
        # A colour bar of 0.2 inch (the rows are drawn at 2.5 cm per inch).
        given <- 0.2 * 2.54 / 2.5 + d
        grDevices::pdf(NULL, width = 7, height = 7)
        drawn <- draws_every_panel(matrix(c(1, 2), nrow = 2),
            c(graphics::lcm((7 - given) * 2.5), graphics::lcm(given * 2.5)),
            list(mai = c(0.1, 0, 0.1, 0)))
        grDevices::dev.off()
        grDevices::pdf(NULL, width = 7, height = 7)
        result <- tryCatch(
            {
                ns$color_plot2(embedding, values, NULL, plt_height = 7 * ns$ppi, plt_width = 7 * ns$ppi,
                    display_legend = TRUE, predicted_height = given)
                "drawn"
            },
            error = conditionMessage
        )
        bar <- graphics::par("fin")[2]
        grDevices::dev.off()
        expect_identical(result, "drawn")
        # Kept where R draws it; raised to 0.21 inch where it does not.
        expect_identical(abs(bar - given * 2.5 / 2.54) < 1e-9, drawn, label = sprintf("bar kept at %g", d))
    }
})

test_that("the app's plots keep the caller's text size, margins and layout", {
    ns <- asNamespace("ClustAssess")
    # the app keeps its colours in this environment; leave them as they were
    had_colours <- exists("discrete_colors", envir = ns$pkg_env, inherits = FALSE)
    old_colours <- if (had_colours) get("discrete_colors", envir = ns$pkg_env)
    on.exit(
        if (had_colours) {
            assign("discrete_colors", old_colours, envir = ns$pkg_env)
        } else {
            rm("discrete_colors", envir = ns$pkg_env)
        },
        add = TRUE
    )
    # Values from 1.2.0: a legend drawn at the caller's text size, and a
    # legend row measured at it.
    grDevices::pdf(NULL, width = 7, height = 7)
    graphics::par(cex = 2)
    ns$only_legend_plot(NULL, grDevices::rainbow(50), seq(0, 1, length.out = 20),
        plt_width = 7 * ns$ppi, text_size = 1.5)
    expect_equal(graphics::par("cex"), 2)
    grDevices::dev.off()
    grDevices::pdf(NULL, width = 7, height = 7)
    graphics::par(cex = 2)
    ns$add_env_variable("discrete_colors", list("2" = c("red", "blue")))
    ns$grouped_boxplot_list(list(c(0.2, 0.3, 0.4), c(0.4, 0.6, 0.8)), c("A", "B"), c(1, 1),
        plt_height = 7 * ns$ppi, plt_width = 7 * ns$ppi, xlab_text = "resolution", text_size = 1.5)
    expect_equal(graphics::par("fin")[2], 1.059, tolerance = 1e-9)
    grDevices::dev.off()

    # The check of a plot of one panel keeps the caller's layout; the check of
    # a layout keeps the caller's text size, margin line height and margins.
    grDevices::pdf(NULL, width = 7, height = 7)
    graphics::par(mfrow = c(2, 2))
    graphics::par(cex = 1.5, mex = 2, mar = c(1, 2, 3, 4))
    expect_true(ns$layout_draws(matrix(1), 1, list(mai = c(0.1, 0, 0.1, 0))))
    expect_equal(graphics::par("mfrow"), c(2, 2))
    expect_equal(graphics::par(c("cex", "mex", "mar")), list(cex = 1.5, mex = 2, mar = c(1, 2, 3, 4)))
    expect_true(ns$layout_draws(matrix(c(1, 2), nrow = 2), c(graphics::lcm(10), graphics::lcm(5))))
    expect_equal(graphics::par(c("cex", "mex", "mar")), list(cex = 1.5, mex = 2, mar = c(1, 2, 3, 4)))
    grDevices::dev.off()

    # In a figure of the caller's own layout, the next one, a plot of one
    # panel is drawn as in 1.2.0, also where the current figure is too small.
    for (advance in c(FALSE, TRUE)) {
        grDevices::pdf(NULL, width = 7, height = 7)
        graphics::layout(matrix(1:2, nrow = 2), heights = if (advance) c(0.1, 6.9) else c(6.9, 0.1))
        graphics::par(mar = rep(0, 4))
        if (advance) {
            graphics::plot.new()
        }
        expect_error(ns$only_legend_plot(NULL, grDevices::rainbow(50), seq(0, 1, length.out = 20),
            plt_width = 7 * ns$ppi, text_size = 1.5), NA)
        expect_equal(graphics::par("mai"), c(0.1, 0, 0.1, 0))
        expect_equal(graphics::par("fin")[2], 6.9)
        grDevices::dev.off()
    }
})

# The graphical parameters of each panel as plot.new() sets it up, and of the
# device after `draw()`.
drawn_states <- function(draw) {
    hooks <- getHook("plot.new")
    on.exit(setHook("plot.new", hooks, "replace"))
    states <- list()
    setHook("plot.new", function() {
        states[[length(states) + 1]] <<- graphics::par(c("fig", "plt", "pin", "mai", "mar", "cex", "mex"))
    })
    draw()
    c(states, list(graphics::par(no.readonly = TRUE)))
}

test_that("the plots that 1.2.0 drew keep the caller's plot region and margins", {
    ns <- asNamespace("ClustAssess")
    had_colours <- exists("discrete_colors", envir = ns$pkg_env, inherits = FALSE)
    old_colours <- if (had_colours) get("discrete_colors", envir = ns$pkg_env)
    on.exit(
        if (had_colours) {
            assign("discrete_colors", old_colours, envir = ns$pkg_env)
        } else {
            rm("discrete_colors", envir = ns$pkg_env)
        },
        add = TRUE
    )
    ns$add_env_variable("discrete_colors", list("2" = c("red", "blue")))
    k_df <- expand.grid(k = 2:5, cl_method = c("Louvain", "Leiden"), stringsAsFactors = FALSE)
    k_df$first_freq <- 8
    k_df$total_freq <- rep(c(5, 10, 20, 15), length.out = nrow(k_df))
    k_df$ecc <- seq(0.6, 0.95, length.out = nrow(k_df))
    k_df$n_partitions <- rep(1:4, length.out = nrow(k_df))
    embedding <- cbind(seq(-1, 1, length.out = 200), sin(seq(0, 6, length.out = 200)))
    values <- seq(0, 1, length.out = 200)
    labels <- rep(c("a", "b"), 100)
    plots <- list(
        k = function() {
            ns$shiny_plot_k_n_partitions(k_df, plt_height = 7 * ns$ppi, display_legend = TRUE,
                pt_size_range = c(1.35, 2.35))
        },
        boxplot = function() {
            ns$grouped_boxplot_list(list(c(0.2, 0.3, 0.4), c(0.4, 0.6, 0.8)), c("A", "B"), c(1, 1),
                plt_height = 7 * ns$ppi, plt_width = 7 * ns$ppi, xlab_text = "resolution")
        },
        umap = function() {
            ns$color_plot2(embedding, labels, color_values = c("red", "blue"), plt_height = 7 * ns$ppi,
                plt_width = 7 * ns$ppi, unique_values = c("a", "b"), display_legend = TRUE)
        },
        colour_bar = function() {
            ns$color_plot2(embedding, values, NULL, plt_height = 7 * ns$ppi, plt_width = 7 * ns$ppi,
                display_legend = TRUE)
        },
        # Two rows of 3.5 inches: with `pty = "s"`, R maps squares of the same
        # size from the margins in both panels and in the right half of the
        # device.
        colour_bar_halves = function() {
            ns$color_plot2(embedding, values, NULL, plt_height = 7 * 2.54 / 2.5 * ns$ppi,
                plt_width = 7 * ns$ppi, display_legend = TRUE, predicted_height = 7 * 2.54 / 2.5 / 2)
        },
        legend = function() {
            ns$only_legend_plot(NULL, grDevices::rainbow(50), values, plt_width = 7 * ns$ppi)
        }
    )
    # A plot region or margins that the caller set, which 1.2.0 kept: margins
    # in inches with a text size that the layout of a plot changes, a region,
    # a size, and a square region, also in margins of the same size on every
    # side.
    settings <- list(
        mai = function() graphics::par(cex = 1.5, mai = c(0.1, 0.3, 0.2, 0.1)),
        plt = function() graphics::par(plt = c(0.3, 0.6, 0.4, 0.6)),
        pin = function() graphics::par(pin = c(1.5, 0.4)),
        square = function() graphics::par(pty = "s"),
        square_margins = function() graphics::par(mar = rep(1, 4), pty = "s")
    )
    draw_with <- function(plot, setting) {
        grDevices::pdf(NULL, width = 7, height = 7)
        on.exit(grDevices::dev.off())
        graphics::par(mar = c(2, 2, 1, 1))
        settings[[setting]]()
        drawn_states(plots[[plot]])
    }
    for (plot in names(plots)) {
        for (setting in names(settings)) {
            # 1.2.0 drew these plots with their legend: the check finds that
            # R draws them, and every panel is drawn as without the check.
            unchecked <- with_mocked_bindings(
                draw_with(plot, setting),
                layout_draws = function(...) TRUE,
                .package = "ClustAssess"
            )
            expect_identical(draw_with(plot, setting), unchecked, label = paste(plot, setting))
        }
    }
    # Margins too small to change the regions in figures of the device's size,
    # but not in a plot row of 0.3 inch, or of 2e-15 (times 2.5 / 2.54) inch,
    # and margins too large to count in lines.
    small_margins <- list(
        lines = list(margins = function() graphics::par(mar = c(0, 0, 1e-16, 1e-16)), row = 0.3),
        inches = list(margins = function() graphics::par(mai = c(0, 0, 3e-17, 3e-17)), row = 0.3),
        tiny_inches = list(margins = function() graphics::par(mai = c(0, 0, 1e-31, 1e-31)), row = 2e-15),
        # Margins of 1e308 inches are more lines than R can count; with a
        # region or a size that the caller set, the plot is drawn.
        huge_inches_region = list(margins = function() {
            graphics::par(mai = rep(1e308, 4))
            graphics::par(plt = c(0.1, 0.9, 0.1, 0.9))
        }, row = 3.5),
        huge_inches_size = list(margins = function() {
            graphics::par(mai = rep(1e308, 4))
            graphics::par(pin = c(2, 1))
        }, row = 3.5),
        huge_lines_region = list(margins = function() {
            graphics::par(mar = rep(1e308, 4))
            graphics::par(plt = c(0.1, 0.9, 0.1, 0.9))
        }, row = 3.5)
    )
    for (margins in names(small_margins)) {
        draw_row <- function() {
            grDevices::pdf(NULL, width = 7, height = 7)
            on.exit(grDevices::dev.off())
            small_margins[[margins]]$margins()
            drawn_states(function() {
                ns$color_plot2(embedding, values, NULL, plt_height = 7 * 2.54 / 2.5 * ns$ppi,
                    plt_width = 7 * ns$ppi, display_legend = TRUE,
                    predicted_height = 7 * 2.54 / 2.5 - small_margins[[margins]]$row)
            })
        }
        unchecked <- with_mocked_bindings(
            draw_row(),
            layout_draws = function(...) TRUE,
            .package = "ClustAssess"
        )
        expect_identical(draw_row(), unchecked, label = paste("small plot row, margins", margins))
    }

    # The check of a colour bar row sets the caller's margins and plot region
    # back as the caller set them, where R draws the row and where it does
    # not: a plot drawn next at a larger text size has the margins and plot
    # region that it would have had, also where the region or size that the
    # caller set is the one that R maps from the margins in the first panel,
    # and where R maps squares of the same size in figures of other shapes
    # (two rows of 3.5 inches).
    bar_rows <- list(drawn = c(6.7, 0.3), not_drawn = c(6.9, 0.1), halves = c(3.5, 3.5) * 2.54 / 2.5)
    first_panel <- function(rows, what) {
        grDevices::pdf(NULL, width = 7, height = 7)
        on.exit(grDevices::dev.off())
        graphics::par(mar = c(2, 2, 1, 1))
        graphics::layout(matrix(c(2, 1), nrow = 2), heights = graphics::lcm(rows * 2.5))
        graphics::par(what)
    }
    callers <- list(
        mai = function(rows) graphics::par(cex = 1.5, mai = c(0.1, 0.3, 0.2, 0.1)),
        plt = function(rows) graphics::par(plt = c(0.3, 0.6, 0.4, 0.6)),
        pin = function(rows) graphics::par(pin = c(1.5, 0.4)),
        square = function(rows) graphics::par(pty = "s"),
        square_margins = function(rows) graphics::par(mar = rep(1, 4), pty = "s"),
        first_plt = function(rows) graphics::par(plt = first_panel(rows, "plt")),
        first_pin = function(rows) graphics::par(pin = first_panel(rows, "pin")),
        figure_pin = function(rows) graphics::par(mar = rep(0, 4), pin = first_panel(rows, "fin")),
        # `pty` set after a region gives the region up.
        region_then_pty = function(rows) {
            graphics::par(mar = rep(0, 4), plt = c(0.15, 0.85, 0.2, 0.7))
            graphics::par(pty = "m")
        }
    )
    regions <- function(caller, rows, check) {
        grDevices::pdf(NULL, width = 7, height = 7)
        on.exit(grDevices::dev.off())
        graphics::par(mar = c(2, 2, 1, 1))
        callers[[caller]](rows)
        if (check) {
            drawn <- ns$layout_draws(matrix(c(1, 2), nrow = 2), graphics::lcm(rows * 2.5),
                list(mai = c(0.1, 0, 0.1, 0)))
            expect_identical(drawn, rows[2] > 0.2 * 2.54 / 2.5, label = paste("check", caller, rows[2]))
        }
        graphics::par(cex = 3)
        graphics::plot.new()
        graphics::par(c("mar", "mai", "plt", "pin"))
    }
    for (caller in names(callers)) {
        for (rows in bar_rows) {
            expect_identical(regions(caller, rows, TRUE), regions(caller, rows, FALSE),
                label = paste("colour bar row", rows[2], caller))
        }
    }
    # Margins too small to change the regions of figures of these sizes show
    # in a plot drawn next with margin lines 1e16 times as high, or in a figure
    # 1e-30 of the device wide and high: margins in lines (also so small that
    # they are 0 inches) or in inches, and such margins with a region that the
    # caller set. Margins of 0 in lines and in inches differ with margin lines
    # of no height (5e-324 times as high), set before or after the check: 0
    # inches are then no finite number of lines. Margins of 1e308 inches are
    # no finite number of lines at any height.
    small_callers <- list(
        lines = function() graphics::par(mar = c(0, 0, 1e-16, 1e-16)),
        lines_and_region = function() {
            graphics::par(mar = c(0, 0, 1e-16, 1e-16))
            graphics::par(plt = c(0, 1, 0, 1))
        },
        lines_of_no_inches = function() graphics::par(mar = c(0, 0, 1e-323, 1e-323)),
        inches = function() graphics::par(mai = c(0, 0, 1e-31, 1e-31)),
        inches_and_region = function() {
            graphics::par(mai = c(0, 0, 1e-31, 1e-31))
            graphics::par(plt = c(0, 1, 0, 1))
        },
        zero_lines = function() graphics::par(mar = rep(0, 4)),
        zero_inches = function() graphics::par(mai = rep(0, 4)),
        zero_lines_no_height = function() graphics::par(mar = rep(0, 4), mex = 5e-324),
        zero_inches_no_height = function() graphics::par(mai = rep(0, 4), mex = 5e-324),
        huge_inches = function() graphics::par(mai = rep(1e308, 4)),
        huge_inches_region = function() {
            graphics::par(mai = rep(1e308, 4))
            graphics::par(plt = c(0.1, 0.9, 0.1, 0.9))
        },
        huge_inches_size = function() {
            graphics::par(mai = rep(1e308, 4))
            graphics::par(pin = c(2, 1))
        },
        huge_lines_region = function() {
            graphics::par(mar = rep(1e308, 4))
            graphics::par(plt = c(0.1, 0.9, 0.1, 0.9))
        }
    )
    next_plots <- list(
        tall_lines = function() graphics::par(mex = 1e16),
        tiny_figure = function() graphics::par(fig = c(0, 1e-30, 0, 1e-30)),
        lines_of_no_height = function() graphics::par(mex = 5e-324)
    )
    small_regions <- function(caller, next_plot, rows, check) {
        grDevices::pdf(NULL, width = 7, height = 7)
        on.exit(grDevices::dev.off())
        small_callers[[caller]]()
        if (check) {
            ns$layout_draws(matrix(c(1, 2), nrow = 2), graphics::lcm(rows * 2.5),
                list(mai = c(0.1, 0, 0.1, 0)))
        }
        next_plots[[next_plot]]()
        tryCatch(
            {
                graphics::plot.new()
                graphics::par(c("mar", "mai", "plt", "pin"))
            },
            error = conditionMessage
        )
    }
    for (caller in names(small_callers)) {
        for (next_plot in names(next_plots)) {
            for (rows in bar_rows) {
                expect_identical(small_regions(caller, next_plot, rows, TRUE),
                    small_regions(caller, next_plot, rows, FALSE),
                    label = paste("colour bar row", rows[2], "small margins", caller, next_plot))
            }
        }
    }
})

# The nearest-neighbour and shared-nearest-neighbour graphs of 300 cells in
# three groups.
three_group_graphs <- function() {
    set.seed(3)
    centers <- matrix(stats::rnorm(3 * 10, sd = 3), 3)
    embedding <- centers[rep(1:3, each = 100), ] + matrix(stats::rnorm(300 * 10), 300)
    rownames(embedding) <- paste0("cell_", seq_len(300))
    ClustAssess:::build_neigh_matrix(
        nn_idx = RANN::nn2(embedding, k = 15)$nn.idx,
        n_neigh = 15,
        cell_names = rownames(embedding),
        graph_type = 2,
        prune_value = -1
    )
}

test_that("the workers build the igraph of a graph once for its consecutive clusterings", {
    ns <- asNamespace("ClustAssess")
    graphs <- three_group_graphs()
    builds <- 0
    counting_env <- new.env(parent = ns)
    counting_env$leiden_graph <- function(object) {
        builds <<- builds + 1
        ns$leiden_graph(object)
    }
    remembering_leiden_graph <- ns$remembering_leiden_graph
    environment(remembering_leiden_graph) <- counting_env
    graph_of <- remembering_leiden_graph()

    snn <- graph_of(graphs$snn)
    expect_false(igraph::is_directed(snn$graph))
    expect_identical(snn$weights, igraph::E(snn$graph)$weight)
    expect_identical(graph_of(graphs$snn), snn)
    expect_equal(builds, 1)
    # A matrix of ones is a directed graph without weights.
    nn <- graph_of(graphs$nn)
    expect_true(igraph::is_directed(nn$graph))
    expect_null(nn$weights)
    expect_equal(builds, 2)
    graph_of(graphs$snn)
    expect_equal(builds, 3)
    # A copy with the same content uses the graph already built.
    copy <- graphs$snn
    copy@x <- copy@x + 0
    expect_identical(graph_of(copy), graph_of(graphs$snn))
    expect_equal(builds, 3)
    # An igraph is used as it is, each time.
    graph_of(snn$graph)
    graph_of(snn$graph)
    expect_equal(builds, 5)
})

test_that("the workers' Leidenbase graph follows a shared matrix that changes in place", {
    testthat::skip_if_not_installed("SharedObject")
    testthat::skip_if_not_installed("leidenbase")
    # A graph of 300 cells without groups, whose partitions depend on its weights.
    set.seed(5)
    embedding <- matrix(stats::rnorm(300 * 10), 300)
    rownames(embedding) <- paste0("cell_", seq_len(300))
    shared <- SharedObject::share(ClustAssess:::build_neigh_matrix(
        nn_idx = RANN::nn2(embedding, k = 15)$nn.idx,
        n_neigh = 15,
        cell_names = rownames(embedding),
        graph_type = 2,
        prune_value = -1
    )$snn)
    SharedObject::setCopyOnWrite(shared@x, FALSE)
    worker_clustering <- ClustAssess:::prepare_clustering_functions_for_workers()
    cluster <- function(clustering) {
        clustering(shared, resolution = 1, seed = 1, algorithm = 6, num_iters = 10, num_starts = 1,
            verbose = FALSE)
    }
    weighted <- cluster(worker_clustering)
    # Weights of one, written in place: a directed graph without weights.
    shared@x[] <- 1
    expect_identical(cluster(worker_clustering), cluster(ClustAssess:::clustering_functions))
    expect_false(identical(cluster(worker_clustering), weighted))
})

test_that("the copy that keys the Leidenbase graph keeps slot attributes and has its own data", {
    ns <- asNamespace("ClustAssess")
    x <- methods::as(Matrix::sparseMatrix(i = 1:4, j = 1:4, x = 1, dims = c(4, 4)), "generalMatrix")
    # A slot with an attribute: its copy is identical, so the graph is reused.
    attr(x@x, "source") <- "probe"
    expect_identical(ns$unshared_copy(x), x)
    graph_of <- ns$remembering_leiden_graph()
    first <- graph_of(x)
    expect_identical(graph_of(x), first)

    testthat::skip_if_not_installed("SharedObject")
    # A slot with a class, shared without copy-on-write and changed in place
    # through another object on the same memory.
    x@p <- SharedObject::share(structure(x@p, levels = as.character(0:4), class = "factor"))
    SharedObject::setCopyOnWrite(x@p, FALSE)
    same_memory <- unclass(x@p)
    SharedObject::setCopyOnWrite(same_memory, FALSE)
    graph_of <- ns$remembering_leiden_graph()
    graph_of(x)
    same_memory[] <- c(0L, 2L, 2L, 3L, 4L)
    edges <- function(graph) igraph::as_edgelist(graph$graph, names = FALSE)
    expect_identical(edges(graph_of(x)), edges(ns$leiden_graph(x)))
})

test_that("algorithm 6 gives the partitions of leidenbase on a new igraph of each graph", {
    testthat::skip_if_not_installed("leidenbase")
    graphs <- three_group_graphs()
    worker_clustering <- ClustAssess:::prepare_clustering_functions_for_workers()
    leidenbase_partition <- function(adjacency, directed, resolution, seed) {
        g <- igraph::graph_from_adjacency_matrix(
            adjacency,
            mode = if (directed) "directed" else "undirected",
            weighted = if (directed) NULL else TRUE
        )
        leidenbase::leiden_find_partition(
            igraph = g,
            partition_type = "RBConfigurationVertexPartition",
            edge_weights = igraph::E(g)$weight,
            resolution_parameter = resolution,
            num_iter = 10,
            seed = seed
        )$membership
    }
    # The worker alternates between the graphs, as the graph step does.
    for (resolution in c(0.5, 1.5)) {
        for (seed in c(1, 42)) {
            for (graph in c("snn", "nn")) {
                expect_identical(
                    worker_clustering(
                        graphs[[graph]],
                        resolution = resolution,
                        seed = seed,
                        algorithm = 6,
                        num_iters = 10,
                        num_starts = 1,
                        verbose = FALSE
                    ),
                    leidenbase_partition(graphs[[graph]], graph == "nn", resolution, seed),
                    label = sprintf("%s graph, resolution %g, seed %d", graph, resolution, seed)
                )
            }
        }
    }
})

test_that("algorithm 6 with several starts gives the same partition for the same seed", {
    testthat::skip_if_not_installed("leidenbase")
    # A graph of 300 cells without groups, whose partitions depend on the seed.
    set.seed(5)
    embedding <- matrix(stats::rnorm(300 * 10), 300)
    rownames(embedding) <- paste0("cell_", seq_len(300))
    snn <- ClustAssess:::build_neigh_matrix(
        nn_idx = RANN::nn2(embedding, k = 15)$nn.idx,
        n_neigh = 15,
        cell_names = rownames(embedding),
        graph_type = 2,
        prune_value = -1
    )$snn
    cluster <- function(seed) {
        ClustAssess:::clustering_functions(snn, resolution = 1, seed = seed, algorithm = 6,
            num_iters = 1, num_starts = 3, verbose = FALSE)
    }
    first <- cluster(1)
    # leidenbase without a seed seeds itself from the clock, in seconds.
    Sys.sleep(1.1)
    expect_identical(cluster(1), first)
    expect_false(identical(cluster(2), first))
    # The random number state is as set.seed() left it.
    set.seed(1)
    seeded_state <- .Random.seed
    cluster(1)
    expect_identical(.Random.seed, seeded_state)
})

test_that("algorithm 6 with one start and no seed still clusters", {
    testthat::skip_if_not_installed("leidenbase")
    membership <- ClustAssess:::leiden_clustering(igraph::make_ring(10), resolution = 1, seed = NULL,
        num_starts = 1)
    expect_length(membership, 10)
})

test_that("feature stability seeds the UMAP of each repetition from its seed", {
    set.seed(7)
    expression <- matrix(
        stats::rpois(40 * 60, 3), 40,
        dimnames = list(paste0("gene_", 1:40), paste0("cell_", 1:60))
    )
    run <- function(outer_seed, umap_arguments = list(n_epochs = 20)) {
        set.seed(outer_seed)
        result <- suppressWarnings(assess_feature_stability(
            data_matrix = expression,
            feature_set = rownames(expression),
            steps = 40,
            feature_type = "genes",
            resolution = c(0.5, 1),
            seed_sequence = c(11, 23),
            graph_reduction_type = "UMAP",
            umap_arguments = umap_arguments,
            compute_umap_embeddings = FALSE,
            verbose = FALSE
        ))
        result$by_steps$genes[["40"]]
    }

    first <- run(10)
    # R's random number state before the call does not change the result
    expect_identical(run(20), first)
    # a seed given in umap_arguments is replaced by the seed of the repetition
    expect_identical(run(30, list(n_epochs = 20, seed = 999)), first)
})

test_that("graph stability assesses a repeated number of neighbours once", {
    set.seed(123)
    embedding <- matrix(
        stats::rnorm(80 * 12), 80,
        dimnames = list(paste0("c", 1:80), paste0("PC", 1:12))
    )
    run <- function(n_neigh_sequence, reduction) {
        suppressWarnings(assess_nn_stability(
            embedding,
            n_neigh_sequence = n_neigh_sequence,
            n_repetitions = 3,
            graph_reduction_type = reduction,
            prune_value = 0.05,
            umap_arguments = list(n_epochs = 25)
        ))
    }
    for (reduction in c("PCA", "UMAP")) {
        expect_identical(run(c(5, 5, 10), reduction), run(c(5, 10), reduction))
    }
})

test_that("comparisons create R's random number state only where 1.2.0 did", {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) {
        old_seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    }
    creates_state <- function(compare) {
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
        }
        compare()
        exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    }
    a <- c(1, 1, 2, 2, 3, 3)
    b <- c(1, 2, 2, 3, 3, 1)
    text <- function(x) letters[x]
    five <- lapply(1:5, function(i) text(if (i %% 2 == 1) a else b))

    # 1.2.0 compared character labels, and ran element_agreement(), in R code
    expect_false(creates_state(function() element_sim(text(a), text(b))))
    expect_false(creates_state(function() element_sim(text(a), b)))
    expect_false(creates_state(function() element_sim_elscore(text(a), text(b))))
    expect_false(creates_state(function() element_consistency(five)))
    expect_false(creates_state(function() weighted_element_consistency(five, c(2, 1, 1, 1, 1))))
    expect_false(creates_state(function() element_sim_matrix(list(text(a), text(b), text(a)))))
    expect_false(creates_state(function() merge_partitions(list(text(a), text(a), text(b)))))
    expect_false(creates_state(function() element_agreement(a, list(b, a))))
    expect_false(creates_state(function() element_agreement(factor(a), list(factor(b)))))
    # and numeric or factor labels in native code
    expect_true(creates_state(function() element_sim(a, b)))
    expect_true(creates_state(function() element_sim(factor(a), factor(b))))
    expect_true(creates_state(function() element_consistency(lapply(1:5, function(i) if (i %% 2 == 1) a else b))))
})

test_that("partitions with missing labels give the results of 1.2.0 where it gave one", {
    valid <- c(1, 1, 2, 2)
    other <- c(1, 2, 1, 2)
    for (mb in list(c("a", NA, "b", "b"), factor(c("a", NA, "b", "b")))) {
        # identical partitions are merged before anything is compared
        merged <- merge_partitions(list(mb, mb), order_logic = "none")
        expect_length(merged$partitions, 1)
        expect_identical(merged$partitions[[1]]$mb, mb)
        expect_equal(merged$partitions[[1]]$freq, 2)
        expect_equal(merged$ecc, rep(1, 4))
        expect_equal(element_consistency(list(mb, mb, mb)), rep(1, 4))
        # character and factor labels are compared where both partitions
        # have a label
        expect_length(merge_partitions(list(mb, valid))$partitions, 1)
        # a partition without missing labels found first is compared as usual
        merged <- merge_partitions(list(valid, mb, other), order_logic = "none")
        expect_equal(vapply(merged$partitions, `[[`, numeric(1), "freq"), c(2, 1))
        expect_equal(merged$ecc, element_consistency(list(valid, valid, other)))
    }
    # 1.2.0 compared numeric labels in native code, which failed on them
    expect_error(merge_partitions(list(c(1, NA, 2, 2), c(1, NA, 2, 2))), "missing values")

    # A hierarchical clustering is compared with the element of a missing label
    # in no cluster; values from 1.2.0.
    hc <- stats::hclust(stats::dist(c(1, 2, 4, 8)))
    scores <- c(0.368878280201, 0.444444444444, 0.446497648926, 0.541449078797)
    for (mb in list(c("a", NA, "b", "b"), c(1, NA, 2, 2))) {
        expect_equal(unname(element_sim_elscore(mb, hc)), scores, tolerance = 1e-10)
        expect_equal(element_sim(mb, hc), 0.450317363092, tolerance = 1e-10)
        expect_equal(suppressWarnings(element_agreement(mb, list(hc))), scores, tolerance = 1e-10)
    }
    for (mb in list(c("a", NA, "b", "b"), c(1, NA, 2, 2), factor(c("a", NA, "b", "b")))) {
        expect_equal(unname(weighted_element_consistency(list(mb, hc), c(1, 1))), scores, tolerance = 1e-10)
    }

    # 1.2.0 compared two vectors of numbers or factors in native code, where a
    # vector whose labels were all missing was one cluster; values from 1.2.0.
    for (all_missing in list(rep(NA_integer_, 4), rep(NaN, 4), factor(rep(NA, 4), levels = "a"))) {
        expect_equal(element_sim(all_missing, valid), 0.5)
        expect_equal(element_sim_elscore(all_missing, valid), rep(0.5, 4))
        expect_equal(weighted_element_consistency(list(all_missing, valid), c(1, 1)), rep(0.5, 4))
        expect_equal(weighted_element_consistency(c(list(all_missing), rep(list(valid), 4))), rep(0.8, 4))
        expect_equal(element_sim_matrix(list(all_missing, valid)), matrix(c(1, 0.5, 0.5, 1), 2))
    }
    expect_length(merge_partitions(list(rep(NA_integer_, 4), rep(NA_integer_, 4)))$partitions, 1)
    # Any missing label stopped where character labels were compared.
    expect_error(element_sim(rep(NA_integer_, 4), c("a", "a", "b", "b")), "missing values")

    # The identity check of 1.2.0 counted the unused levels of a factor: two
    # such factors with a missing label were left to compare, which stopped R.
    unused <- factor(c("a", NA, "b", "b"), levels = c("a", "b", "unused"))
    expect_error(merge_partitions(list(unused, unused)), "missing values")
    # A factor without levels, all its labels missing, was the same partition
    # as itself; one with a level was not, and was compared as one cluster.
    expect_length(merge_partitions(list(factor(rep(NA, 4)), factor(rep(NA, 4))))$partitions, 1)
    with_level <- factor(rep(NA, 4), levels = "a")
    merged <- merge_partitions(list(with_level, with_level))
    expect_equal(vapply(merged$partitions, `[[`, numeric(1), "freq"), c(1, 1))
    expect_equal(merged$ecc, rep(1, 4))

    # 1.2.0 compared each pair of a partition list on its own: a vector with
    # missing labels was compared with a tree as above; values from 1.2.0.
    hc6 <- stats::hclust(stats::dist(c(1, 2, 4, 8, 16, 32)))
    merged <- merge_partitions(list(hc6, c(1, NA, 2, 3, 3, 1), hc6), ecs_thresh = 0.9)
    expect_equal(vapply(merged$partitions, `[[`, numeric(1), "freq"), c(2, 1))
    expect_equal(unname(merged$ecc), c(0.561324735177, 0.62962962963, 0.490709564862, 0.546025826201,
        0.572292678288, 0.648829859017), tolerance = 1e-10)
    # A vector whose labels are all missing is not the same partition as a
    # matrix of another length.
    overlap <- matrix(c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1), nrow = 6)
    merged <- merge_partitions(list(rep(NA_real_, 6), overlap))
    expect_equal(vapply(merged$partitions, `[[`, numeric(1), "freq"), c(1, 1))
    expect_equal(unname(merged$ecc), rep(4 / 9, 6))

    # 1.2.0 merged a list with missing labels by its own check, in its order:
    # a factor with an unused level and a numeric vector of the same partition
    # stayed apart, and a factor with a missing label could join the second.
    # That merge comes first; the partitions without missing labels are then
    # merged. The consistency is that of 1.2.0.
    u <- factor(c(1, 1, 2, 2), levels = 1:3)
    v <- c(1, 1, 2, 2)
    m <- factor(c("a", NA, "b", "b"))
    for (order in list(c(1, 2, 3), c(2, 1, 3), c(2, 3, 1))) {
        merged <- merge_partitions(list(u, v, m)[order])
        expect_equal(vapply(merged$partitions, `[[`, numeric(1), "freq"), 3)
        expect_equal(merged$ecc, rep(1, 4))
        expect_equal(element_consistency(list(u, v, m)[order]), rep(1, 4))
    }
    for (freq in list(c(0, 1, 3), c(2, 2, 4))) {
        merged <- merge_partitions(Map(function(mb, f) list(mb = mb, freq = f), list(u, v, m), freq))
        expect_equal(vapply(merged$partitions, `[[`, numeric(1), "freq"), sum(freq))
        expect_equal(merged$ecc, rep(1, 4))
    }
    # Where the factor with a missing label stayed apart, 1.2.0 failed.
    for (order in list(c(1, 3, 2), c(3, 1, 2), c(3, 2, 1))) {
        expect_error(merge_partitions(list(u, v, m)[order]), "missing values")
    }
    # Vectors of different lengths are not the same partition.
    expect_error(merge_partitions(list(NA_real_, rep(NA_real_, 4))), "same length")
})
