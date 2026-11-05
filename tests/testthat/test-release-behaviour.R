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
