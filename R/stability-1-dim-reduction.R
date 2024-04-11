# TODO close the socket connections in case of interruptions
# TODO use the marker JSI as indicator for stabiliity instead of / together with the ECC

#' Assess the stability for configurations of feature types and sizes
#'
#' @description Evaluate the stability of clusterings obtained
#' based on incremental subsets of a given feature set.
#'
#' @param data_matrix A data matrix having the features on the rows
#' and the observations on the columns.
#' @param feature_set A set of feature names that can be found on the rownames
#' of the data matrix.
#' @param steps Vector containing the sizes of the subsets; negative values will
#' be interpreted as using all features.
#' @param feature_type A name associated to the feature_set.
#' @param resolution A vector containing the resolution values used for
#' clustering.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user. Defaults
#' to `100`.
#' @param seed_sequence A custom seed sequence; if the value is NULL, the
#' sequence will be built starting from 1 with a step of 100. Defaults to
#' `NULL`.
#' @param graph_reduction_type The graph reduction type, denoting if the graph
#' should be built on either the PCA or the UMAP embedding. Defaults to `PCA`.
#' @param ecs_thresh The ECS threshold used for merging similar clusterings. We
#' recommend using the 1 value. Defaults to `1`.
#' @param matrix_processing A function that will be used to process the data
#' matrix
#' by using a dimensionality reduction technique. The function should have
#' one parameter, the data matrix, and should return an embedding describing the
#' reduced space. By default, the function will use the precise PCA method with
#' `prcomp`.
#' @param umap_arguments A list containing the arguments that will be passed
#' to the UMAP function. Refer to the `uwot::umap` function for more details.
#' @param prune_value Argument indicating whether to prune the SNN graph. If the
#' value is 0, the graph won't be pruned. If the value is between 0 and 1, the
#' edges with weight under the pruning value will be removed. If the value is
#' -1, the highest pruning value will be calculated automatically and used.
#' @param clustering_algorithm An index indicating which community detection
#' algorithm will be used: Louvain (1), Louvain refined (2), SLM (3) or
#' Leiden (4). More details can be found in the Seurat's
#' `FindClusters` function.
#' @param clustering_arguments A list containing the arguments that will be
#' passed to the community detection algorithm, such as the number of iterations
#' and the number of starts. Refer to the Seurat's `FindClusters` function
#' for more details.
#' @param verbose A boolean indicating if the intermediate progress will
#' be printed or not.
#'
#' @return A list having one field associated with a step value. Each step
#' contains a list with three fields:
#'
#' * ecc - the EC-Consistency of the partitions obtained on all repetitions
#' * embedding - one UMAP embedding generated on the feature subset
#' * most_frequent_partition - the most common partition obtained across
#' repetitions
#'
#' @note The algorithm assumes that the feature_set is already sorted when
#' performing the subsetting based on the steps values. For example, if the
#' user wants to analyze highly variable feature set, they should provide them
#' sorted by their variability.
#'
#'
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(100 * 10), runif(100 * 10, min = 3, max = 4)),
#'     nrow = 200, byrow = TRUE
#' )
#' rownames(expr_matrix) <- as.character(1:200)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- assess_feature_stability(
#'     data_matrix = t(expr_matrix),
#'     feature_set = colnames(expr_matrix),
#'     steps = 5,
#'     feature_type = "feature_name",
#'     resolution = c(0.1, 0.5, 1),
#'     n_repetitions = 10,
#'     umap_arguments = list(
#'         # the following parameters are used by the umap function
#'         # and are not mandatory
#'         n_neighbors = 3,
#'         approx_pow = TRUE,
#'         n_epochs = 0,
#'         init = "random",
#'         min_dist = 0.3
#'     ),
#'     clustering_algorithm = 1
#' )
#' plot_feature_overall_stability_boxplot(feature_stability_result)
assess_feature_stability <- function(data_matrix,
                                     feature_set,
                                     steps,
                                     feature_type,
                                     resolution,
                                     n_repetitions = 100,
                                     seed_sequence = NULL,
                                     graph_reduction_type = "PCA",
                                     ecs_thresh = 1,
                                     matrix_processing = function(dt_mtx, actual_npcs = 30, ...) {
                                         # WARNING using more threads will lead to slightly different results
                                         # check using RhpcBLASctl::blas_get_num_procs()
                                         actual_npcs <- min(actual_npcs, ncol(dt_mtx) %/% 2)
                                         RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
                                         embedding <- stats::prcomp(x = dt_mtx, rank. = actual_npcs)$x
                                         RhpcBLASctl::blas_set_num_threads(1)

                                         rownames(embedding) <- rownames(dt_mtx)
                                         colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))

                                         return(embedding)
                                     },
                                     umap_arguments = list(),
                                     prune_value = -1,
                                     clustering_algorithm = 1,
                                     clustering_arguments = list(),
                                     verbose = FALSE) {
    # TODO create a function that does checks of the parameters
    # NOTE irbla using all cores; fix with RhpcBLASctl::blas_set_num_threads(1)? seems to be fixed, keep under observation
    # check parameters
    if (!is.matrix(data_matrix) && !inherits(data_matrix, "Matrix")) {
        stop("the data matrix parameter should be a matrix")
    }

    # transpose the matrix to have observations on rows and features on columns
    if (inherits(data_matrix, "dgCMatrix")) {
        data_matrix <- Matrix::t(data_matrix)
    } else {
        data_matrix <- t(data_matrix)
    }
    ncells <- nrow(data_matrix)

    if (!is.character(feature_set)) {
        stop("feature_set parameter should be a vector of strings")
    }

    if (!all(feature_set %in% colnames(data_matrix))) {
        stop("all features from the feature_set should be found in the data_matrix")
    }

    if (!is.numeric(steps)) {
        stop("steps parameter should be numeric")
    }
    # convert steps to integers
    steps <- as.integer(steps)

    if (!is.character(feature_type)) {
        stop("feature_type parameter should be a string")
    }

    if (!is.numeric(n_repetitions) || length(n_repetitions) > 1) {
        stop("n_repetitions parameter should be numeric")
    }
    # convert n_repetitions to integers
    n_repetitions <- as.integer(n_repetitions)

    if (!(graph_reduction_type %in% c("PCA", "UMAP"))) {
        stop("graph_reduction_type parameter should take one of these values: 'PCA' or 'UMAP'")
    }

    if (!is.numeric(ecs_thresh)) {
        stop("ecs_thresh parameter should be numeric")
    }

    if (!is.numeric(clustering_algorithm) || length(clustering_algorithm) > 1 || !(clustering_algorithm %in% 1:4)) {
        stop("algorithm should be a number between 1 and 4")
    }

    ncores <- foreach::getDoParWorkers()

    partitions_list <- list()

    object_name <- feature_type
    steps_ecc_list <- list()
    steps_ecc_list[[object_name]] <- list()
    embedding_list <- list()
    embedding_list[[object_name]] <- list()
    pca_list <- list()
    pca_list[[object_name]] <- list()

    # create a seed sequence if it's not provided
    if (is.null(seed_sequence)) {
        seed_sequence <- seq(
            from = 1,
            by = 100,
            length.out = n_repetitions
        )
    } else {
        n_repetitions <- length(seed_sequence)
    }

    # convert the negative steps to the size of the feature
    steps[steps <= 0 |
        steps > length(feature_set)] <- length(feature_set)
    steps <- sort(unique(steps))

    umap_arguments <- process_umap_arguments(umap_arguments, ncells)

    if (verbose) {
        pb <- progress::progress_bar$new(
            format = glue::glue("{feature_type} - :featuresize [:bar] eta: :eta  total elapsed: :elapsed"),
            total = length(steps),
            show_after = 0,
            clear = FALSE,
            width = 80
        )
    }

    clustering_arguments <- process_clustering_arguments(clustering_arguments, clustering_algorithm)

    for (step in steps) {
        if (verbose) {
            pb$tick(0, tokens = list(featuresize = step))
        }
        used_features <- feature_set[1:step]
        step <- as.character(step)

        partitions_list[[step]] <- list()
        steps_ecc_list[[object_name]][[step]] <- list()

        embedding <- matrix_processing(data_matrix[, used_features])

        if (is.null(rownames(embedding))) {
            rownames(embedding) <- rownames(data_matrix)
        }

        if (is.null(colnames(embedding))) {
            colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))
        }

        pca_list[[object_name]][[step]] <- embedding # nolint

        # build the SNN graph before if PCA
        if (graph_reduction_type == "PCA") {
            used_n_neigh <- min(20, ncells - 1)

            if (prune_value >= 0) {
                neigh_matrix <- getNNmatrix(
                    RANN::nn2(embedding, k = used_n_neigh)$nn.idx,
                    used_n_neigh,
                    0,
                    prune_value
                )$snn
            } else {
                neigh_matrix <- get_highest_prune_param(
                    nn_matrix = getNNmatrix(
                        RANN::nn2(embedding, k = used_n_neigh)$nn.idx,
                        used_n_neigh,
                        0,
                        -1
                    )$nn,
                    n_neigh = used_n_neigh
                )$adj_matrix
            }
            rownames(neigh_matrix) <- rownames(embedding)
            colnames(neigh_matrix) <- rownames(embedding)
        }

        # the variables needed in each PSOCK process
        needed_vars <- c(
            "graph_reduction_type",
            "umap_arguments",
            "clustering_arguments",
            "resolution",
            "prune_value",
            "ncells"
        )

        if (graph_reduction_type == "PCA") {
            needed_vars <- c(needed_vars, "shared_neigh_matrix")
            if (ncores > 1) {
                shared_neigh_matrix <- SharedObject::share(neigh_matrix)
            } else {
                shared_neigh_matrix <- neigh_matrix
            }
        } else {
            needed_vars <- c(needed_vars, "shared_embedding")
            if (ncores > 1) {
                shared_embedding <- SharedObject::share(embedding)
            } else {
                shared_embedding <- embedding
            }
        }

        all_vars <- ls()
        seed <- 0

        # send the name of the umap arguments
        partitions_list[[step]] <- foreach::foreach(
            seed = seed_sequence,
            .inorder = FALSE,
            .noexport = all_vars[!(all_vars %in% needed_vars)]
        ) %dopar% {
            if (graph_reduction_type == "UMAP") {
                row_names <- rownames(shared_embedding)
                # calculate the umap embedding
                if ("seed" %in% names(umap_arguments)) {
                    umap_arguments$seed <- seed
                }

                ump_embedding <- do.call(
                    uwot::umap,
                    c(
                        list(X = shared_embedding),
                        umap_arguments
                    )
                )
                colnames(ump_embedding) <- paste0("UMAP_", seq_len(ncol(ump_embedding)))
                rownames(ump_embedding) <- row_names
                gc()

                # build the nn and snn graphs and apply graph clustering
                used_n_neigh <- min(20, ncells - 1)
                if (prune_value >= 0) {
                    umap_neigh_matrix <- getNNmatrix(
                        RANN::nn2(ump_embedding, k = used_n_neigh)$nn.idx,
                        used_n_neigh,
                        0,
                        prune_value
                    )$snn
                } else {
                    umap_neigh_matrix <- get_highest_prune_param(
                        nn_matrix = getNNmatrix(
                            RANN::nn2(ump_embedding, k = used_n_neigh)$nn.idx,
                            used_n_neigh,
                            0,
                            -1
                        )$nn,
                        n_neigh = used_n_neigh
                    )$adj_matrix
                }
                rownames(umap_neigh_matrix) <- row_names
                colnames(umap_neigh_matrix) <- row_names

                # apply graph clustering to the snn graph
                cluster_results <- lapply(resolution, function(r) {
                    cl_res <- do.call(
                        clustering_functions,
                        c(
                            list(
                                object = umap_neigh_matrix,
                                resolution = r,
                                seed = seed
                            ),
                            clustering_arguments
                        )
                    )

                    # return the clustering
                    list(
                        mb = cl_res,
                        freq = 1,
                        seed = seed
                    )
                })

                names(cluster_results) <- as.character(resolution)
                return(cluster_results)
            }

            # apply graph clustering to the snn graph
            cluster_results <- lapply(resolution, function(r) {
                cl_res <- do.call(
                    clustering_functions,
                    c(
                        list(
                            object = shared_neigh_matrix,
                            resolution = r,
                            seed = seed
                        ),
                        clustering_arguments
                    )
                )

                # return the clustering
                list(
                    mb = cl_res,
                    freq = 1,
                    seed = seed
                )
            })

            names(cluster_results) <- as.character(resolution)
            cluster_results
        }

        if (graph_reduction_type == "PCA") {
            shared_neigh_matrix <- SharedObject::unshare(neigh_matrix)
        } else {
            shared_embedding <- SharedObject::unshare(embedding)
        }

        set.seed(42) # to match Seurat's default
        embedding <- base::do.call(
            uwot::umap,
            c(
                list(X = embedding),
                umap_arguments
            )
        )

        colnames(embedding) <- paste0("UMAP_", seq_len(ncol(embedding)))
        rownames(embedding) <- rownames(data_matrix)
        gc()
        embedding_list[[object_name]][[step]] <- embedding

        # merge the identical partitions into the same object
        partitions_list[[step]] <- lapply(as.character(resolution), function(r) {
            merge_partitions(
                lapply(partitions_list[[step]], function(x) {
                    x[[r]]
                }),
                ecs_thresh = ecs_thresh,
                order = TRUE
            )
        })

        names(partitions_list[[step]]) <- as.character(resolution)

        # update the returned object with a list containing the ecc,
        # the most frequent partition and the number of different partitions
        steps_ecc_list[[object_name]][[step]] <- lapply(as.character(resolution), function(r) {
            list(
                ecc = weighted_element_consistency(
                    lapply(partitions_list[[step]][[r]], function(x) {
                        x$mb
                    }),
                    sapply(partitions_list[[step]][[r]], function(x) {
                        x$freq
                    }) # NOTE ECS should be fast enough without parallelization
                ),
                most_frequent_partition = partitions_list[[step]][[r]][[1]],
                n_different_partitions = length(partitions_list[[step]][[r]])
            )
        })

        names(steps_ecc_list[[object_name]][[step]]) <- as.character(resolution)

        if (verbose) {
            pb$tick(tokens = list(featuresize = step))
        }
    }

    steps <- as.character(steps)
    resolution <- as.character(resolution)

    incremental_ecs_list <- list()
    incremental_ecs_list[[object_name]] <- list()

    if (length(steps) > 1) {
        for (i in 2:length(steps)) {
            temp_list <- lapply(resolution, function(r) {
                element_sim_elscore(
                    steps_ecc_list[[object_name]][[steps[i - 1]]][[r]]$most_frequent_partition$mb,
                    steps_ecc_list[[object_name]][[steps[i]]][[r]]$most_frequent_partition$mb
                )
            })
            names(temp_list) <- resolution
            incremental_ecs_list[[object_name]][[paste(steps[i - 1], steps[i], sep = "-")]] <- temp_list
        }
    }

    list(
        by_steps = steps_ecc_list,
        incremental = incremental_ecs_list,
        embedding_list = embedding_list,
        pca_list = pca_list
    )
}

#' Per resolution Feature Stability Boxplot
#'
#' @description Display EC consistency for each feature set and for each step.
#' Above each boxplot there is a number representing
#' the step (or the size of the subset). The ECC values are extracted depdening
#' on the resolution value provided by the user.
#'
#' @param feature_object_list An object or a concatenation of objects returned
#' by the `assess_feature_stability` method
#' @param resolution The resolution value for which the ECC will be extracted.
#' @param violin_plot If TRUE, the function will return a violin plot instead
#' of a boxplot. Default is `FALSE`.
#' @param text_size The size of the labels above boxplots
#' @param boxplot_width Used for adjusting the width of the boxplots; the value
#' will be passed to the `width` argument of the `ggplot2::geom_boxplot` method.
#' @param dodge_width Used for adjusting the horizontal position of the boxplot;
#' the value will be passed to the `width` argument of the
#' `ggplot2::position_dodge` method.
#' @param return_df If TRUE, the function will return the ECS values as a
#' dataframe. Default is `FALSE`.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(100 * 10), runif(100 * 10, min = 3, max = 4)),
#'     nrow = 200, byrow = TRUE
#' )
#' rownames(expr_matrix) <- as.character(1:200)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- assess_feature_stability(
#'     data_matrix = t(expr_matrix),
#'     feature_set = colnames(expr_matrix),
#'     steps = 5,
#'     feature_type = "feature_name",
#'     resolution = c(0.1, 0.5, 1),
#'     n_repetitions = 10,
#'     umap_arguments = list(
#'         # the following parameters are used by the umap function
#'         # and are not mandatory
#'         n_neighbors = 3,
#'         approx_pow = TRUE,
#'         n_epochs = 0,
#'         init = "random",
#'         min_dist = 0.3
#'     ),
#'     clustering_algorithm = 1
#' )
#' plot_feature_per_resolution_stability_boxplot(feature_stability_result, 0.5)
plot_feature_per_resolution_stability_boxplot <- function(feature_object_list,
                                                          resolution,
                                                          violin_plot = FALSE,
                                                          text_size = 4,
                                                          boxplot_width = 0.4,
                                                          dodge_width = 0.7,
                                                          return_df = FALSE) {
    resolution <- as.character(resolution)
    min_index <- -1 # number of steps that will be displayed on the plot
    feature_object_list <- feature_object_list$by_steps
    final_melt_df <- NULL

    # create a dataframe based on the object returned by `get_feature_stability`
    for (config_name in names(feature_object_list)) {
        list_ecc <- lapply(feature_object_list[[config_name]], function(x) {
            # NOTE as.numeric is necessary for the shiny app, which creates arrays, not numeric vectors
            as.numeric(x[[resolution]]$ecc)
        })

        melt_object <- reshape2::melt(list_ecc)
        melt_object$L1 <- factor(
            melt_object$L1,
            levels = stringr::str_sort(unique(melt_object$L1), numeric = TRUE)
        )
        melt_object[["feature_set"]] <- rep(config_name, nrow(melt_object))

        temp_steps_df <- data.frame(
            step = levels(melt_object$L1),
            index = 1:nlevels(melt_object$L1),
            stringsAsFactors = FALSE
        )
        temp_steps_df$feature_set <- rep(config_name, nrow(temp_steps_df))

        levels(melt_object$L1) <- 1:nlevels(melt_object$L1)

        if (min_index == -1) {
            final_steps_df <- temp_steps_df
            final_melt_df <- melt_object

            min_index <- nrow(temp_steps_df)
        } else {
            if (nrow(temp_steps_df) < min_index) {
                min_index <- nrow(temp_steps_df)
            }

            final_steps_df <- rbind(final_steps_df, temp_steps_df)
            final_melt_df <- rbind(final_melt_df, melt_object)
        }
    }

    # given that the input object can have multiple configurations with different
    # number of steps, we will use only the first `min_index` steps
    final_melt_df <- final_melt_df %>% dplyr::filter(as.numeric(.data$L1) <= min_index)
    final_steps_df <- final_steps_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)

    final_melt_df$feature_set <- factor(final_melt_df$feature_set, levels = names(feature_object_list))
    final_steps_df$feature_set <- factor(final_steps_df$feature_set, levels = names(feature_object_list))

    names(final_melt_df) <- c("ecc", "step_index", "feature_set")

    if (return_df) {
        return(final_melt_df)
    }

    # generate the coordinates where the sizes of the steps will be displayed
    text_position <-
        stats::aggregate(ecc ~ step_index + feature_set, final_melt_df, max)
    text_position$ecc <- max(text_position$ecc) + (max(final_melt_df$ecc) - min(final_melt_df$ecc)) / 100

    geom_function <- ifelse(violin_plot, ggplot2::geom_violin, ggplot2::geom_boxplot)
    # return the ggplot object
    ggplot2::ggplot(
        final_melt_df,
        ggplot2::aes(
            x = .data$step_index,
            y = .data$ecc,
            fill = .data$feature_set
        )
    ) +
        geom_function(
            position = ggplot2::position_dodge(width = dodge_width),
            width = boxplot_width,
            outlier.shape = NA
        ) +
        ggplot2::geom_text(
            data = text_position,
            ggplot2::aes(label = final_steps_df$step),
            position = ggplot2::position_dodge(width = dodge_width),
            size = text_size
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::xlab("# features") +
        ggplot2::ylab("EC consistency") +
        ggplot2::ggtitle(glue::glue("Stability for resolution = {resolution}"))
}

#' Overall Feature Stability Boxplot
#'
#' @description Display EC consistency for each feature set and for each step.
#' Above each boxplot there is a number representing
#' the step (or the size of the subset). The ECC values are extracted for each
#' resolution value and summarized using the `summary_function` parameter.
#'
#' @param feature_object_list An object or a concatenation of objects returned
#' by the `assess_feature_stability` method
#' @param summary_function The function that will be used to summarize the ECC
#' values. Defaults to `median`.
#' @param text_size The size of the labels above boxplots
#' @param boxplot_width Used for adjusting the width of the boxplots; the value
#' will be passed to the `width` argument of the `ggplot2::geom_boxplot` method.
#' @param dodge_width Used for adjusting the horizontal position of the boxplot;
#' the value will be passed to the `width` argument of the
#' `ggplot2::position_dodge` method.
#' @param return_df If TRUE, the function will return the ECS values as a
#' dataframe. Default is `FALSE`.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(100 * 10), runif(100 * 10, min = 3, max = 4)),
#'     nrow = 200, byrow = TRUE
#' )
#' rownames(expr_matrix) <- as.character(1:200)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- assess_feature_stability(
#'     data_matrix = t(expr_matrix),
#'     feature_set = colnames(expr_matrix),
#'     steps = 5,
#'     feature_type = "feature_name",
#'     resolution = c(0.1, 0.5, 1),
#'     n_repetitions = 10,
#'     umap_arguments = list(
#'         # the following parameters are used by the umap function
#'         # and are not mandatory
#'         n_neighbors = 3,
#'         approx_pow = TRUE,
#'         n_epochs = 0,
#'         init = "random",
#'         min_dist = 0.3
#'     ),
#'     clustering_algorithm = 1
#' )
#' plot_feature_overall_stability_boxplot(feature_stability_result)
plot_feature_overall_stability_boxplot <- function(feature_object_list,
                                                   summary_function = stats::median,
                                                   text_size = 4,
                                                   boxplot_width = 0.4,
                                                   dodge_width = 0.7,
                                                   return_df = FALSE) {
    min_index <- -1 # number of steps that will be displayed on the plot
    feature_object_list <- feature_object_list$by_steps
    final_melt_df <- NULL

    # create a dataframe based on the object returned by `get_feature_stability`
    for (config_name in names(feature_object_list)) {
        melt_object <- reshape2::melt(lapply(feature_object_list[[config_name]], function(by_step) {
            sapply(by_step, function(by_res) {
                summary_function(by_res$ecc)
            })
        }))
        melt_object$L1 <- factor(
            melt_object$L1,
            levels = stringr::str_sort(unique(melt_object$L1), numeric = TRUE)
        )
        melt_object[["feature_set"]] <- rep(config_name, nrow(melt_object))

        temp_steps_df <- data.frame(
            step = levels(melt_object$L1),
            index = 1:nlevels(melt_object$L1),
            stringsAsFactors = FALSE
        )
        temp_steps_df$feature_set <- rep(config_name, nrow(temp_steps_df))

        levels(melt_object$L1) <- 1:nlevels(melt_object$L1)

        if (min_index == -1) {
            final_steps_df <- temp_steps_df
            final_melt_df <- melt_object

            min_index <- nrow(temp_steps_df)
        } else {
            if (nrow(temp_steps_df) < min_index) {
                min_index <- nrow(temp_steps_df)
            }

            final_steps_df <- rbind(final_steps_df, temp_steps_df)
            final_melt_df <- rbind(final_melt_df, melt_object)
        }
    }

    # given that the input object can have multiple configurations with different
    # number of steps, we will use only the first `min_index` steps
    final_melt_df <- final_melt_df %>% dplyr::filter(as.numeric(.data$L1) <= min_index)
    final_steps_df <- final_steps_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)

    final_melt_df$feature_set <- factor(final_melt_df$feature_set, levels = names(feature_object_list))
    final_steps_df$feature_set <- factor(final_steps_df$feature_set, levels = names(feature_object_list))

    names(final_melt_df) <- c("ecc", "step_index", "feature_set")

    if (return_df) {
        return(final_melt_df)
    }

    # generate the coordinates where the sizes of the steps will be displayed
    text_position <-
        stats::aggregate(ecc ~ step_index + feature_set, final_melt_df, max)
    text_position$ecc <- max(text_position$ecc) + (max(final_melt_df$ecc) - min(final_melt_df$ecc)) / 100

    # return the ggplot object
    ggplot2::ggplot(
        final_melt_df,
        ggplot2::aes(
            x = .data$step_index,
            y = .data$ecc,
            fill = .data$feature_set
        )
    ) +
        ggplot2::geom_boxplot(
            position = ggplot2::position_dodge(width = dodge_width),
            width = boxplot_width,
            outlier.shape = NA
        ) +
        ggplot2::geom_text(
            data = text_position,
            ggplot2::aes(label = final_steps_df$step),
            position = ggplot2::position_dodge(width = dodge_width),
            size = text_size
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::xlab("# features") +
        ggplot2::ylab("EC consistency") +
        ggplot2::ggtitle("Overall stability")
}

#' Feature Stability - Cluster Membership Facet Plot
#'
#' @description Display a facet of plots where each subpanel is associated with
#' a feature set and illustrates the distribution of the most frequent partition
#' over the UMAP embedding.
#'
#' @param feature_object_list An object or a concatenation of objects returned
#' by the `assess_feature_stability` method
#' @param resolution The resolution value for which the ECS will be extracted.
#' @param text_size The size of the cluster label
#' @param n_facet_cols The number of facet's columns.
#' @param point_size The size of the points displayed on the plot.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(100 * 10), runif(50 * 10, min = 3, max = 4)),
#'     nrow = 150, byrow = TRUE
#' )
#' rownames(expr_matrix) <- as.character(1:150)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- assess_feature_stability(
#'     data_matrix = t(expr_matrix),
#'     feature_set = colnames(expr_matrix),
#'     steps = 5,
#'     feature_type = "feature_name",
#'     resolution = c(0.1, 0.5, 1),
#'     n_repetitions = 10,
#'     clustering_algorithm = 1
#' )
#' plot_feature_stability_mb_facet(
#'     feature_stability_result,
#'     0.5,
#'     point_size = 2
#' )
plot_feature_stability_mb_facet <- function(feature_object_list,
                                            resolution,
                                            text_size = 5,
                                            n_facet_cols = 3,
                                            point_size = 0.3) {
    resolution <- as.character(resolution)
    if (!is.numeric(text_size) || length(text_size) > 1) {
        stop("text_size parameter should be numeric")
    }

    first_temp <- TRUE
    umap_object_list <- feature_object_list$embedding_list
    feature_object_list <- feature_object_list$by_steps

    for (config_name in names(feature_object_list)) {
        for (steps in names(feature_object_list[[config_name]])) {
            temp_df <- data.frame(
                x = umap_object_list[[config_name]][[steps]][, 1],
                y = umap_object_list[[config_name]][[steps]][, 2],
                mb = feature_object_list[[config_name]][[steps]][[resolution]]$most_frequent_partition$mb
            )

            temp_df[["config_name"]] <- rep(config_name, nrow(temp_df))
            temp_df[["steps"]] <- rep(steps, nrow(temp_df))

            if (first_temp) {
                first_temp <- FALSE
                final_df <- temp_df
            } else {
                final_df <- rbind(final_df, temp_df)
            }
        }
    }

    final_df$steps <- factor(final_df$steps, levels = stringr::str_sort(unique(final_df$steps), numeric = TRUE))
    final_df$mb <- factor(final_df$mb)
    text.labs <- final_df %>%
        dplyr::group_by(.data$config_name, .data$steps, .data$mb) %>%
        dplyr::summarise(
            mean_x = stats::median(.data$x),
            mean_y = stats::median(.data$y)
        )

    ggplot2::ggplot(
        final_df,
        ggplot2::aes(
            x = .data$x,
            y = .data$y,
            color = .data$mb
        )
    ) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::theme_classic() +
        ggplot2::geom_text(
            data = text.labs,
            ggplot2::aes(
                x = .data$mean_x,
                y = .data$mean_y,
                label = .data$mb
            ),
            size = text_size,
            color = "black"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "none",
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::labs(x = "UMAP_1", y = "UMAP_2") +
        ggplot2::facet_wrap(~ config_name + steps, ncol = n_facet_cols)
}

#' Feature Stability - EC Consistency Facet Plot
#'
#' @description Display a facet of plots where each subpanel is associated with
#' a feature set and illustrates the distribution of the EC consistency score
#' over the UMAP embedding.
#'
#' @param feature_object_list An object or a concatenation of objects returned
#' by the `assess_feature_stability` method
#' @param resolution The resolution value for which the ECS will be extracted.
#' @param n_facet_cols The number of facet's columns.
#' @param point_size The size of the points displayed on the plot.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(100 * 10), runif(50 * 10, min = 3, max = 4)),
#'     nrow = 150, byrow = TRUE
#' )
#' rownames(expr_matrix) <- as.character(1:150)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- assess_feature_stability(
#'     data_matrix = t(expr_matrix),
#'     feature_set = colnames(expr_matrix),
#'     steps = 5,
#'     feature_type = "feature_name",
#'     resolution = c(0.1, 0.5, 1),
#'     n_repetitions = 10,
#'     clustering_algorithm = 1
#' )
#' plot_feature_stability_ecs_facet(
#'     feature_stability_result,
#'     0.5,
#'     point_size = 2
#' )
plot_feature_stability_ecs_facet <- function(feature_object_list,
                                             resolution,
                                             n_facet_cols = 3,
                                             point_size = 0.3) {
    resolution <- as.character(resolution)
    first_temp <- TRUE
    umap_object_list <- feature_object_list$embedding_list
    feature_object_list <- feature_object_list$by_steps

    for (config_name in names(feature_object_list)) {
        for (steps in names(feature_object_list[[config_name]])) {
            temp_df <- data.frame(
                x = umap_object_list[[config_name]][[steps]][, 1],
                y = umap_object_list[[config_name]][[steps]][, 2],
                ecc = feature_object_list[[config_name]][[steps]][[resolution]]$ecc
            )

            temp_df[["config_name"]] <- rep(config_name, nrow(temp_df))
            temp_df[["steps"]] <- rep(steps, nrow(temp_df))

            if (first_temp) {
                first_temp <- FALSE
                final_df <- temp_df
            } else {
                final_df <- rbind(final_df, temp_df)
            }
        }
    }

    final_df$steps <- factor(final_df$steps, levels = stringr::str_sort(unique(final_df$steps), numeric = TRUE))

    ggplot2::ggplot(
        final_df,
        ggplot2::aes(
            x = .data$x,
            y = .data$y,
            color = .data$ecc
        )
    ) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::labs(
            color = "ECC",
            x = "UMAP_1",
            y = "UMAP_2"
        ) +
        ggplot2::facet_wrap(~ config_name + steps, ncol = n_facet_cols)
}

#' Per resolution - Feature Stability Incremental Boxplot
#'
#' @description Perform an incremental ECS between two consecutive feature
#' steps. The ECS values are extracted only for a specified resolution value.
#'
#'
#' @param feature_object_list An object or a concatenation of objects returned
#' by the `assess_feature_stability` method.
#' @param resolution The resolution value for which the ECS will be extracted.
#' @param dodge_width Used for adjusting the horizontal position of the boxplot;
#' the value will be passed to the `width` argument of the
#' `ggplot2::position_dodge` method.
#' @param text_size The size of the labels above boxplots.
#' @param boxplot_width Used for adjusting the width of the boxplots; the value
#' will be passed to the `width` argument of the `ggplot2::geom_boxplot` method.
#' @param return_df If TRUE, the function will return the ECS values as a
#' dataframe. Default is `FALSE`.
#'
#' @return A ggplot2 object with ECS distribution will be displayed as a
#' boxplot. Above each boxplot there will be a pair of numbers representing the
#' two steps that are compared.
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(50 * 10), runif(50 * 10, min = 3, max = 4)),
#'     nrow = 100, byrow = TRUE
#' )
#' rownames(expr_matrix) <- as.character(1:100)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- assess_feature_stability(
#'     data_matrix = t(expr_matrix),
#'     feature_set = colnames(expr_matrix),
#'     steps = c(5, 10),
#'     feature_type = "feature_name",
#'     resolution = c(0.1, 0.5),
#'     n_repetitions = 3,
#'     umap_arguments = list(
#'         # the following parameters are used by the umap function
#'         # and are not mandatory
#'         n_neighbors = 3,
#'         approx_pow = TRUE,
#'         n_epochs = 0,
#'         init = "random",
#'         min_dist = 0.3
#'     ),
#'     clustering_algorithm = 1
#' )
#' plot_feature_per_resolution_stability_incremental(feature_stability_result, 0.1)
plot_feature_per_resolution_stability_incremental <- function(feature_object_list,
                                                              resolution,
                                                              dodge_width = 0.7,
                                                              text_size = 4,
                                                              boxplot_width = 0.4,
                                                              return_df = FALSE) {
    resolution <- as.character(resolution)
    if (!is.numeric(dodge_width) || length(dodge_width) > 1) {
        stop("dodge_width parameter should be numeric")
    }

    min_index <- -1
    feature_object_list <- feature_object_list$incremental
    first_df <- TRUE
    ecs_df <- NULL

    for (config_name in names(feature_object_list)) {
        n_pairs <- length(feature_object_list[[config_name]])
        pairs.names <- names(feature_object_list[[config_name]])
        first_from_pair <- sapply(pairs.names, function(x) {
            strsplit(x, "-")[[1]][1]
        })
        second_from_pair <- sapply(pairs.names, function(x) {
            strsplit(x, "-")[[1]][2]
        })

        # treat the case with only one step
        index <- 0
        for (i in stringr::str_order(first_from_pair, numeric = TRUE)) {
            index <- index + 1
            first_n_steps <- first_from_pair[i]
            second_n_steps <- second_from_pair[i]

            temp_df <- data.frame(
                ecs = feature_object_list[[config_name]][[i]][[resolution]],
                index = index,
                feature_set = config_name
            )

            if (first_df) {
                first_df <- FALSE

                steps_df <- data.frame(
                    step = paste(
                        first_n_steps,
                        second_n_steps,
                        sep = "-\n"
                    ),
                    index = index,
                    feature_set = config_name,
                    stringsAsFactors = FALSE
                )

                ecs_df <- temp_df
            } else {
                ecs_df <- rbind(ecs_df, temp_df)

                steps_df <- rbind(steps_df, c(
                    paste(
                        first_n_steps,
                        second_n_steps,
                        sep = "-\n"
                    ),
                    index,
                    config_name
                ))
            }
        }

        if (min_index == -1 || n_pairs < min_index) {
            min_index <- n_pairs
        }
    }

    if (is.null(ecs_df)) {
        return(ggplot2::ggplot() +
            ggplot2::theme_void())
    }

    ecs_df <- ecs_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)
    steps_df <- steps_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)

    ecs_df$feature_set <- factor(ecs_df$feature_set, levels = names(feature_object_list))
    steps_df$feature_set <- factor(steps_df$feature_set, levels = names(feature_object_list))

    ecs_df$index <- factor(ecs_df$index)

    text_position <-
        stats::aggregate(ecs ~ index + feature_set, ecs_df, max)
    text_position$ecs <- max(text_position$ecs) + (max(ecs_df$ecs) - min(ecs_df$ecs)) / 25

    if (return_df) {
        return(ecs_df)
    }

    ggplot2::ggplot(
        ecs_df,
        ggplot2::aes(
            x = .data$index,
            y = .data$ecs,
            fill = .data$feature_set
        )
    ) +
        ggplot2::geom_boxplot(
            position = ggplot2::position_dodge(width = dodge_width),
            width = boxplot_width,
            outlier.shape = NA
        ) +
        ggplot2::geom_text(
            data = text_position,
            ggplot2::aes(label = steps_df$step),
            position = ggplot2::position_dodge(width = dodge_width),
            size = text_size
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::xlab("# features") +
        ggplot2::ylab("EC similiarity") +
        ggplot2::ggtitle(glue::glue("Incremental stability for resolution = {resolution}"))
}

#' Overall Feature Stability Incremental Boxplot
#'
#' @description Perform an incremental ECS between two consecutive feature
#' steps. The ECS values are extracted for every resolution value and summarized
#' using a function (e.g. median, mean, etc.).
#'
#'
#' @param feature_object_list An object or a concatenation of objects returned
#' by the `assess_feature_stability` method.
#' @param summary_function The function used to summarize the ECS values.
#' Default is `median`.
#' @param dodge_width Used for adjusting the horizontal position of the boxplot;
#' the value will be passed to the `width` argument of the
#' `ggplot2::position_dodge` method.
#' @param text_size The size of the labels above boxplots.
#' @param boxplot_width Used for adjusting the width of the boxplots; the value
#' will be passed to the `width` argument of the `ggplot2::geom_boxplot` method.
#' @param return_df If TRUE, the function will return the ECS values as
#' a dataframe. Default is `FALSE`.
#'
#' @return A ggplot2 object with ECS distribution will be displayed as a
#' boxplot. Above each boxplot there will be a pair of numbers representing the
#' two steps that are compared.
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(50 * 10), runif(50 * 10, min = 3, max = 4)),
#'     nrow = 100, byrow = TRUE
#' )
#' rownames(expr_matrix) <- as.character(1:100)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- assess_feature_stability(
#'     data_matrix = t(expr_matrix),
#'     feature_set = colnames(expr_matrix),
#'     steps = c(5, 10),
#'     feature_type = "feature_name",
#'     resolution = c(0.1, 0.5),
#'     n_repetitions = 3,
#'     umap_arguments = list(
#'         # the following parameters are used by the umap function
#'         # and are not mandatory
#'         n_neighbors = 3,
#'         approx_pow = TRUE,
#'         n_epochs = 0,
#'         init = "random",
#'         min_dist = 0.3
#'     ),
#'     clustering_algorithm = 1
#' )
#' plot_feature_overall_stability_incremental(feature_stability_result)
plot_feature_overall_stability_incremental <- function(feature_object_list,
                                                       summary_function = stats::median,
                                                       dodge_width = 0.7,
                                                       text_size = 4,
                                                       boxplot_width = 0.4,
                                                       return_df = FALSE) {
    if (!is.numeric(dodge_width) || length(dodge_width) > 1) {
        stop("dodge_width parameter should be numeric")
    }

    min_index <- -1
    feature_object_list <- feature_object_list$incremental
    first_df <- TRUE
    ecs_df <- NULL

    for (config_name in names(feature_object_list)) {
        n_pairs <- length(feature_object_list[[config_name]])
        pairs.names <- names(feature_object_list[[config_name]])
        first_from_pair <- sapply(pairs.names, function(x) {
            strsplit(x, "-")[[1]][1]
        })
        second_from_pair <- sapply(pairs.names, function(x) {
            strsplit(x, "-")[[1]][2]
        })

        # treat the case with only one step
        index <- 0
        for (i in stringr::str_order(first_from_pair, numeric = TRUE)) {
            index <- index + 1
            first_n_steps <- first_from_pair[i]
            second_n_steps <- second_from_pair[i]

            temp_df <- data.frame(
                ecs = sapply(feature_object_list[[config_name]][[i]], summary_function),
                index = index,
                feature_set = config_name
            )

            if (first_df) {
                first_df <- FALSE

                steps_df <- data.frame(
                    step = paste(
                        first_n_steps,
                        second_n_steps,
                        sep = "-\n"
                    ),
                    index = index,
                    feature_set = config_name,
                    stringsAsFactors = FALSE
                )

                ecs_df <- temp_df
            } else {
                ecs_df <- rbind(ecs_df, temp_df)

                steps_df <- rbind(steps_df, c(
                    paste(
                        first_n_steps,
                        second_n_steps,
                        sep = "-\n"
                    ),
                    index,
                    config_name
                ))
            }
        }

        if (min_index == -1 || n_pairs < min_index) {
            min_index <- n_pairs
        }
    }

    if (is.null(ecs_df)) {
        return(ggplot2::ggplot() +
            ggplot2::theme_void())
    }

    ecs_df <- ecs_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)
    steps_df <- steps_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)

    ecs_df$feature_set <- factor(ecs_df$feature_set, levels = names(feature_object_list))
    steps_df$feature_set <- factor(steps_df$feature_set, levels = names(feature_object_list))

    ecs_df$index <- factor(ecs_df$index)

    text_position <-
        stats::aggregate(ecs ~ index + feature_set, ecs_df, max)
    text_position$ecs <- max(text_position$ecs) + (max(ecs_df$ecs) - min(ecs_df$ecs)) / 25

    if (return_df) {
        return(ecs_df)
    }

    ggplot2::ggplot(
        ecs_df,
        ggplot2::aes(
            x = .data$index,
            y = .data$ecs,
            fill = .data$feature_set
        )
    ) +
        ggplot2::geom_boxplot(
            position = ggplot2::position_dodge(width = dodge_width),
            width = boxplot_width,
            outlier.shape = NA
        ) +
        ggplot2::geom_text(
            data = text_position,
            ggplot2::aes(label = steps_df$step),
            position = ggplot2::position_dodge(width = dodge_width),
            size = text_size
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::xlab("# features") +
        ggplot2::ylab("EC similiarity") +
        ggplot2::ggtitle("Overall incremental stability")
}
