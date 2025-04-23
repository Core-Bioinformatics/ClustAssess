# 3. The assessment of the stability for clustering paramters
# (such as clustering method, number of clusters)

algorithm_names <- c("Louvain", "Louvain.refined", "SLM", "Leiden")
algorithm_names_mapping <- list(
    Louvain = 1,
    Louvain.refined = 2,
    SLM = 3,
    Leiden = 4
)

#' Merge Partitions from different Resolutions
#'
#' @description Merge partitions obtained with different resolution values. The
#' partitions will be grouped based on the number of clusters. The identical
#' partitions will be merged into a single partition by updating the frequency
#' using the `merge_partitions` method.
#'
#' @param res_obj A list associated to a configuration field from the object
#' returned by the `assess_clustering_importance` method.
#'
#' @return A list having one field assigned to each number of clusters. A number
#' of cluster will contain a list of all merged partitions. To avoid duplicates,
#' `merged_partitions` with threshold 1 is applied.
merge_resolutions <- function(res_obj) {
    clusters_obj <- list()
    ncores <- foreach::getDoParWorkers()

    # concatenate all partitions having the same number of clusters into a list
    for (res.val in names(res_obj)) {
        for (k in names(res_obj[[res.val]]$clusters)) {
            if (!(k %in% names(clusters_obj))) {
                clusters_obj[[k]] <- res_obj[[res.val]]$clusters[[k]]$partitions
            } else {
                clusters_obj[[k]] <- c(
                    clusters_obj[[k]],
                    res_obj[[res.val]]$clusters[[k]]$partitions
                )
            }
        }
    }

    k_vals <- names(clusters_obj)

    # if (ncores > 1) {
    #     shared_clusters_obj <- SharedObject::share(clusters_obj)
    # } else {
    #     shared_clusters_obj <- clusters_obj
    # }

    all_vars <- ls()
    needed_vars <- c("shared_clusters_obj")

    i <- 1

    # merge identical partitions from each list for a fixed k
    clusters_obj <- foreach::foreach(
        i = seq_along(clusters_obj),
        .inorder = TRUE,
        .noexport = all_vars[!(all_vars %in% needed_vars)]
    ) %do% {
        ClustAssess::merge_partitions(clusters_obj[[i]])
    }

    names(clusters_obj) <- k_vals

    clusters_obj
}

#' Assessment of Stability for Graph Clustering
#' @description Evaluates the stability of different graph clustering methods
#' in the clustering pipeline. The method will iterate through different values of
#' the resolution parameter and compare, using the EC Consistency score, the
#' partitions obtained at different seeds.
#'
#' @param graph_adjacency_matrix A square adjacency matrix based on which an igraph
#' object will be built. The matrix should have rownames and colnames that correspond
#' to the names of the cells.
#' @param resolution A sequence of resolution values. The resolution parameter
#' controls the coarseness of the clustering. The higher the resolution, the more
#' clusters will be obtained. The resolution parameter is used in the community
#' detection algorithms.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user. Defaults to `100`.
#' @param seed_sequence A custom seed sequence; if the value is NULL, the
#' sequence will be built starting from 1 with a step of 100.
#' @param ecs_thresh The ECS threshold used for merging similar clusterings.
#' @param clustering_algorithm An index or a list of indexes indicating which community detection
#' algorithm will be used: Louvain (1), Louvain refined (2), SLM (3) or Leiden (4).
#' More details can be found in the Seurat's `FindClusters` function. Defaults to `1:3`.
#' @param clustering_arguments A list of additional arguments that will be passed to the
#' clustering method. More details can be found in the Seurat's `FindClusters` function.
#' @param verbose Boolean value used for displaying the progress bar.
#'
#' @return A list having two fields:
#'
#' * `all` - a list that contains, for each clustering method and each resolution
#' value, the EC consistency between the partitions obtained by changing the seed
#' * `filtered` - similar to `all`, but for each configuration, we determine the
#' number of clusters that appears the most and use only the partitions with this
#' size
#'
#' @export
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_embedding <- matrix(runif(100 * 30), nrow = 100)
#' rownames(pca_embedding) <- paste0("cell_", seq_len(nrow(pca_embedding)))
#' colnames(pca_embedding) <- paste0("PC_", 1:30)
#'
#'
#' adj_matrix <- getNNmatrix(
#'     RANN::nn2(pca_embedding, k = 10)$nn.idx,
#'     10,
#'     0,
#'     -1
#' )$nn
#' rownames(adj_matrix) <- paste0("cell_", seq_len(nrow(adj_matrix)))
#' colnames(adj_matrix) <- paste0("cell_", seq_len(ncol(adj_matrix)))
#'
#' # alternatively, the adj_matrix can be calculated
#' # using the `Seurat::FindNeighbors` function.
#'
#' clust_diff_obj <- assess_clustering_stability(
#'     graph_adjacency_matrix = adj_matrix,
#'     resolution = c(0.5, 1),
#'     n_repetitions = 10,
#'     clustering_algorithm = 1:2,
#'     verbose = TRUE
#' )
#' plot_clustering_overall_stability(clust_diff_obj)
assess_clustering_stability <- function(graph_adjacency_matrix,
                                        resolution,
                                        n_repetitions = 100,
                                        seed_sequence = NULL,
                                        ecs_thresh = 1,
                                        clustering_algorithm = 1:3,
                                        clustering_arguments = list(),
                                        verbose = TRUE) {
    # BUG there are some performance issues, make sure the main runtime is caused by the clustering method
    # TODO add leidenbase implementation for the algorithm = 4
    ncores <- foreach::getDoParWorkers()
    # check the parameters
    if (!is.numeric(resolution)) {
        stop("resolution parameter should be numeric")
    }

    if (!is.numeric(ecs_thresh) || length(ecs_thresh) > 1) {
        stop("ecs_thresh pa              ld be numer,
              algorithm = alg_index                        ic")
    }

    if (!is.numeric(n_repetitions) || length(n_repetitions) > 1) {
        stop("n_repetitions parameter should be numeric")
    }
    # convert n_repetitions to an integer
    n_repetitions <- as.integer(n_repetitions)

    if (!is.numeric(clustering_algorithm) || !(all(clustering_algorithm %in% 1:4))) {
        stop("algorithm should be a vector of numbers between 1 and 4")
    }

    if (!is.logical(verbose)) {
        stop("verbose parameter should be logical")
    }

    # create a seed sequence if it's not provided
    if (is.null(seed_sequence)) {
        seed_sequence <- seq(
            from = 1,
            by = 100,
            length.out = n_repetitions
        )
    } else {
        if (!is.numeric(seed_sequence)) {
            stop("seed_sequence parameter should be numeric")
        }

        seed_sequence <- as.integer(seed_sequence)
    }

    result_object <- list()

    # additional arguments used by the clustering method

    # if (min(algorithm) < 4) {
    if (ncores > 1 && is_package_installed("SharedObject")) {
        graph_adjacency_matrix_shared <- SharedObject::share(graph_adjacency_matrix)
    } else {
        graph_adjacency_matrix_shared <- graph_adjacency_matrix
    }
    # }

    # the variables needed in the PSOCK processes
    seurat_clustering <- seurat_clustering
    needed_vars <- c(
        "res",
        "seurat_clustering",
        "graph_adjacency_matrix_shared",
        "alg_index",
        "current_clustering_arguments"
    )
    package_needed <- c()
    if (4 %in% clustering_algorithm) {
        package_needed <- c(package_needed, "leiden")
    }
    all_vars <- ls()
    seed <- 0



    for (alg_index in clustering_algorithm) {
        # FIXME edit this once you will be able to share an igraph
        # if (algorithm == 4) {
        #     g <- igraph::graph_from_adjacency_matrix(
        #         adjmatrix = graph_adjacency_matrix,
        #         mode = ifelse(Matrix::isSymmetric(graph_adjacency_matrix), "undirected", "directed"),
        #         weighted = TRUE
        #     )

        #     if (ncores > 1) {
        #         graph_adjacency_matrix_shared <- SharedObject::share(g)
        #     } else {
        #         graph_adjacency_matrix_shared <- g
        #     }
        # }

        alg_name <- algorithm_names[alg_index]
        current_clustering_arguments <- process_clustering_arguments(clustering_arguments, alg_index)

        if (verbose) {
            pb <- progress::progress_bar$new(
                format = glue::glue("{alg_name} - res :res [:bar] eta: :eta  total elapsed: :elapsed"),
                total = length(resolution),
                clear = FALSE,
                width = 80
            )
            pb$tick(0)
        }

        result_object[[alg_name]] <- list()
        for (res in resolution) {
            if (verbose) {
                pb$tick(0, tokens = list(res = res))
            }

            different_partitions_temp <- foreach::foreach(
                seed = seed_sequence,
                .inorder = FALSE,
                .noexport = all_vars[!(all_vars %in% needed_vars)],
                .export = needed_vars,
                .packages = package_needed
            ) %dopar% {
                # apply the clustering, which should return a membership vector
                do.call(
                    clustering_functions,
                    c(
                        list(
                            object = graph_adjacency_matrix_shared,
                            resolution = res,
                            seed = seed
                        ),
                        current_clustering_arguments
                    )
                )
            }
            different_partitions <- list()

            # group the partitions by the number of clusters
            for (i in seq_along(seed_sequence)) {
                k <- as.character(length(unique(different_partitions_temp[[i]])))

                if (!(k %in% names(different_partitions))) {
                    different_partitions[[k]] <- list()
                    different_partitions[[k]][[1]] <- list(
                        mb = different_partitions_temp[[i]],
                        freq = 1,
                        seed = seed_sequence[i]
                    )
                } else {
                    index <- length(different_partitions[[k]])
                    different_partitions[[k]][[index + 1]] <- list(
                        mb = different_partitions_temp[[i]],
                        freq = 1,
                        seed = seed_sequence[i]
                    )
                }
            }

            different_partitions <- different_partitions[stringr::str_sort(names(different_partitions), numeric = TRUE)]

            # merge the partitions using the ecs threshold
            unique_partitions <- list()
            for (k in names(different_partitions)) {
                different_partitions[[k]] <- merge_partitions(
                    partition_list = different_partitions[[k]],
                    ecs_thresh = ecs_thresh,
                    order_logic = "freq" # TODO change to avg agreement
                )

                unique_partitions <- c(
                    unique_partitions,
                    different_partitions[[k]]$partitions
                )

                # compute the EC-consistency of the partition list
                # ec_consistency <- weighted_element_consistency(
                #   clustering_list = lapply(different_partitions[[k]]$partitions, function(x) {
                #     x$mb
                #   }),
                #   weights = sapply(different_partitions[[k]]$partitions, function(x) {
                #     x$freq
                #   }),
                #   calculate_sim_matrix = TRUE
                # )

                # different_partitions[[k]][["ecc"]] <- ec_consistency
            }

            if (verbose) {
                pb$tick(tokens = list(res = res))
            }

            if (length(unique_partitions) == 1) {
                result_object[[alg_name]][[as.character(res)]] <- list(
                    clusters = different_partitions,
                    ecc = rep(1, length(unique_partitions[[1]]$mb))
                )
                next
            }

            result_object[[alg_name]][[as.character(res)]] <- list(
                clusters = different_partitions,
                ecc = weighted_element_consistency(
                    lapply(unique_partitions, function(x) {
                        x$mb
                    }),
                    sapply(unique_partitions, function(x) {
                        x$freq
                    })
                )
            )

        }
    }

    on.exit(
        {
            if (ncores > 1 && is_package_installed("SharedObject")) {
                graph_adjacency_matrix_shared <- SharedObject::unshare(graph_adjacency_matrix)
                rm(graph_adjacency_matrix_shared)
            }
            rm(graph_adjacency_matrix)
            gc()
        },
        add = TRUE
    )

    split_by_k <- lapply(algorithm_names[clustering_algorithm], function(alg_name) {
        merge_resolutions(result_object[[alg_name]])
    })

    names(split_by_k) <- algorithm_names[clustering_algorithm]

    list(
        split_by_resolution = result_object,
        split_by_k = split_by_k
    )
}

#' Clustering Method per value Stability Boxplot
#'
#' @description Display EC consistency across clustering methods, calculated
#' for each value of the resolution parameter or the number of clusters.
#'
#' @param clust_object An object returned by the
#' `assess_clustering_stability` method.
#' @param value_type A string that specifies the type of value that was used
#' for grouping the partitions and calculating the ECC score. It can be either
#' `k` or `resolution`. Defaults to `k`.
#'
#' @return A ggplot2 object with the EC consistency distributions grouped by
#' the clustering methods. Higher consistency indicates a more stable clustering.
#' The X axis is decided by the `value_type` parameter.
#' @export
#'
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_embedding <- matrix(runif(100 * 30), nrow = 100)
#' rownames(pca_embedding) <- paste0("cell_", seq_len(nrow(pca_embedding)))
#' colnames(pca_embedding) <- paste0("PC_", 1:30)
#'
#'
#' adj_matrix <- getNNmatrix(
#'     RANN::nn2(pca_embedding, k = 10)$nn.idx,
#'     10,
#'     0,
#'     -1
#' )$nn
#' rownames(adj_matrix) <- paste0("cell_", seq_len(nrow(adj_matrix)))
#' colnames(adj_matrix) <- paste0("cell_", seq_len(ncol(adj_matrix)))
#'
#' # alternatively, the adj_matrix can be calculated
#' # using the `Seurat::FindNeighbors` function.
#'
#' clust_diff_obj <- assess_clustering_stability(
#'     graph_adjacency_matrix = adj_matrix,
#'     resolution = c(0.5, 1),
#'     n_repetitions = 10,
#'     clustering_algorithm = 1:2,
#'     verbose = FALSE
#' )
#' plot_clustering_per_value_stability(clust_diff_obj)
plot_clustering_per_value_stability <- function(clust_object,
                                                value_type = c("k", "resolution")) {
    value_type <- value_type[value_type %in% c("k", "resolution")]
    # TODO add empty boxplots for the missing values for k at least to help differentiate

    if (length(value_type) > 1) {
        value_type <- value_type[1]
    }

    if (length(value_type) == 0) {
        stop("`value_type` should contain either `k` or `resolution`")
    }

    ecc_vals <- lapply(clust_object[[paste0("split_by_", value_type)]], function(by_alg) {
        lapply(by_alg, function(by_value) {
            by_value$ecc
        })
    })

    melted_df <- reshape2::melt(ecc_vals)
    colnames(melted_df) <- c("ecc", value_type, "method")
    melted_df$method <- factor(melted_df$method)
    unique_vals <- stringr::str_sort(unique(melted_df[[value_type]]), numeric = TRUE)
    melted_df[[value_type]] <- factor(melted_df[[value_type]], levels = unique_vals)


    ggplot2::ggplot(
        melted_df,
        ggplot2::aes(x = .data[[value_type]], y = .data$ecc, fill = .data$method)
    ) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste0("Clustering stability per ", value_type))
}



#' Clustering Method Overall Stability Boxplot
#'
#' @description Display EC consistency across clustering methods by summarising
#' the distribution of the EC consistency for each number of clusters.
#'
#' @param clust_object An object returned by the
#' `assess_clustering_stability` method.
#' @param value_type A string that specifies the type of value that was used
#' for grouping the partitions and calculating the ECC score. It can be either
#' `k` or `resolution`. Defaults to `k`.
#' @param summary_function The function that will be used to summarize the
#' distribution of the ECC values obtained for each number of clusters. Defaults
#' to `median`.
#'
#' @return A ggplot2 object with the EC consistency distributions grouped by
#' the clustering methods. Higher consistency indicates
#' a more stable clustering.
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_embedding <- matrix(runif(100 * 30), nrow = 100)
#' rownames(pca_embedding) <- paste0("cell_", seq_len(nrow(pca_embedding)))
#' colnames(pca_embedding) <- paste0("PC_", 1:30)
#'
#'
#' adj_matrix <- getNNmatrix(
#'     RANN::nn2(pca_embedding, k = 10)$nn.idx,
#'     10,
#'     0,
#'     -1
#' )$nn
#' rownames(adj_matrix) <- paste0("cell_", seq_len(nrow(adj_matrix)))
#' colnames(adj_matrix) <- paste0("cell_", seq_len(ncol(adj_matrix)))
#'
#' # alternatively, the adj_matrix can be calculated
#' # using the `Seurat::FindNeighbors` function.
#'
#' clust_diff_obj <- assess_clustering_stability(
#'     graph_adjacency_matrix = adj_matrix,
#'     resolution = c(0.5, 1),
#'     n_repetitions = 10,
#'     clustering_algorithm = 1:2,
#'     verbose = FALSE
#' )
#' plot_clustering_overall_stability(clust_diff_obj)
plot_clustering_overall_stability <- function(clust_object,
                                              value_type = c("k", "resolution"),
                                              summary_function = stats::median) {
    value_type <- value_type[value_type %in% c("k", "resolution")]

    if (length(value_type) > 1) {
        value_type <- value_type[1]
    }

    if (length(value_type) == 0) {
        stop("`value_type` should contain either `k` or `resolution`")
    }

    ecc_vals <- lapply(
        clust_object[[paste0("split_by_", value_type)]],
        function(by_alg) {
            lapply(by_alg, function(by_value) {
                summary_function(by_value$ecc)
            })
        }
    )

    melted_df <- reshape2::melt(ecc_vals)
    colnames(melted_df) <- c("ecc", value_type, "method")
    melted_df$method <- factor(melted_df$method)
    melted_df[[value_type]] <- factor(melted_df[[value_type]])

    ggplot2::ggplot(
        melted_df,
        ggplot2::aes(x = .data$method, y = .data$ecc, fill = .data$method)
    ) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste0("Overall clustering stability grouped by ", value_type))
}

#' Clustering Method Stability Facet Plot
#'
#' @description Display the distribution of the EC consistency for each
#' clustering method and each resolution value on a given embedding The `all`
#' field of the object returned by the `get_clustering_difference_object` method is used.
#'
#' @param clust_object An object returned by the
#' `assess_clustering_stability` method.
#' @param embedding An embedding (only the first two dimensions will be used for
#' visualization).
#' @param low_limit The lowest value of ECC that will be displayed on the embedding.
#' @param high_limit The highest value of ECC that will be displayed on the embedding.
#' @param grid Boolean value indicating whether the facet should be a grid (where each
#' row is associated with a resolution value and each column with a clustering method) or
#' a wrap.
#'
#' @return A ggplot2 object.
#' #TODO should export
#'
#'
#' @examples
#' # FIXME fix the examples
#' # set.seed(2021)
#' # # create an artificial PCA embedding
#' # pca_embedding <- matrix(runif(100 * 30), nrow = 100)
#' # rownames(pca_embedding) <- as.character(1:100)
#' # colnames(pca_embedding) <- paste0("PCA_", 1:30)
#'
#' # adj_matrix <- Seurat::FindNeighbors(pca_embedding,
#' #     k.param = 10,
#' #     nn.method = "rann",
#' #     verbose = FALSE,
#' #     compute.SNN = FALSE
#' # )$nn
#' # clust_diff_obj <- assess_clustering_stability(
#' #     graph_adjacency_matrix = adj_matrix,
#' #     resolution = c(0.5, 1),
#' #     n_repetitions = 10,
#' #     algorithm = 1:2,
#' #     verbose = FALSE
#' # )
#' # plot_clustering_difference_facet(clust_diff_obj, pca_embedding)
plot_clustering_difference_facet <- function(clust_object,
                                             embedding,
                                             low_limit = 0,
                                             high_limit = 1,
                                             grid = TRUE) {
    # check parameters

    # TODO adapt to the new object
    return(1)
    if (!is.logical(grid)) {
        stop("grid parameter should be logical")
    }

    if (!is.numeric(low_limit)) {
        stop("low_limit parameter should be numeric")
    }

    if (!is.numeric(high_limit)) {
        stop("high_limit parameter should be numeric")
    }

    npoints <- nrow(embedding)
    if (length(clust_object$all[[1]][[1]]) != npoints) {
        stop("The provided embedding and the consistency arrays must have the same number of elements!")
    }

    if (ncol(embedding) < 2) {
        stop("The embedding should have at least two dimensions!")
    }

    melt_obj <- reshape2::melt(clust_object$all)

    n_embeeding_repetitions <- nrow(melt_obj) / npoints

    melt_obj["x"] <- rep(embedding[, 1], n_embeeding_repetitions)
    melt_obj["y"] <- rep(embedding[, 2], n_embeeding_repetitions)

    melt_obj["value"][melt_obj["value"] < low_limit] <- low_limit
    melt_obj["value"][melt_obj["value"] > high_limit] <- high_limit

    return_plot <- ggplot2::ggplot(
        melt_obj,
        ggplot2::aes(
            x = .data$x,
            y = .data$y,
            color = .data$value
        )
    ) +
        ggplot2::geom_point() +
        ggplot2::scale_color_viridis_c() +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::labs(
            color = "ECC",
            x = colnames(embedding)[1],
            y = colnames(embedding)[2]
        )

    if (grid) {
        return(return_plot + ggplot2::facet_grid(L2 ~ L1))
    }

    return_plot + ggplot2::facet_wrap(~ L2 + L1)
}



#' Correspondence Between Resolution and the Number of Clusters
#'
#' @description For each configuration provided in the clust_object, display
#' what number of clusters appear for different values of the
#' resolution parameters.
#'
#'
#' @param clust_object An object returned by the
#' `assess_clustering_stability` method.
#' @param colour_information String that specifies the information type that
#' will be illustrated using gradient colour: either `freq_part` for the
#' frequency of the most common partition or `ecc` for the
#' Element-Centric Consistency of the partitions obtained when the the number
#' of clusters is fixed. Defaults to `ecc`.
#' @param dodge_width Used for adjusting the distance between the boxplots
#' representing a clustering method. Defaults to `0.3`.
#' @param pt_size_range Indicates the minimum and the maximum size a point
#' on the plot can have. Defaults to `c(1.5, 4)`.
#' @param summary_function The function that will be used to summarize the
#' distribution of the ECC values obtained for each number of clusters. Defaults
#' to `median`.
#'
#' @return A ggplot2 object. Different shapes of points indicate different
#' parameter configuration, while the color
#' illustrates the frequency of the most common partition or the Element-Centric Consistency
#' of the partitions. The frequency is calculated
#' as the fraction between the number of total appearances of partitions with a specific
#' number of clusters and resolution value and the number of
#' runs. The size
#' illustrates the frequency of the most common partition with *k* clusters relative to the
#' partitions obtained with the same resolution value and have *k* clusters.
#' @export
#'
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_embedding <- matrix(runif(100 * 30), nrow = 100)
#' rownames(pca_embedding) <- paste0("cell_", seq_len(nrow(pca_embedding)))
#' colnames(pca_embedding) <- paste0("PC_", 1:30)
#'
#'
#' adj_matrix <- getNNmatrix(
#'     RANN::nn2(pca_embedding, k = 10)$nn.idx,
#'     10,
#'     0,
#'     -1
#' )$nn
#' rownames(adj_matrix) <- paste0("cell_", seq_len(nrow(adj_matrix)))
#' colnames(adj_matrix) <- paste0("cell_", seq_len(ncol(adj_matrix)))
#'
#' # alternatively, the adj_matrix can be calculated
#' # using the `Seurat::FindNeighbors` function.
#'
#' clust_diff_obj <- assess_clustering_stability(
#'     graph_adjacency_matrix = adj_matrix,
#'     resolution = c(0.5, 1),
#'     n_repetitions = 10,
#'     clustering_algorithm = 1:2,
#'     verbose = FALSE
#' )
#' plot_k_resolution_corresp(clust_diff_obj)
plot_k_resolution_corresp <- function(clust_object,
                                      colour_information = c("ecc", "freq_k"),
                                      dodge_width = 0.3,
                                      pt_size_range = c(1.5, 4),
                                      summary_function = stats::median) {
    # TODO check the colors and the vertical lines, try to help the user
    if (length(colour_information) > 1) {
        colour_information <- colour_information[1]
    }

    if (!(colour_information %in% c("ecc", "freq_k"))) {
        stop("colour_information can be either `ecc` or `freq_k`")
    }

    clust_object <- clust_object$split_by_resolution

    # use the names of the fields from the list
    res_object_names <- names(clust_object)

    # create a dataframe that contains the number of cases when,
    # for a given resolution, a number of clusters was obtained
    for (i in seq_along(clust_object)) {
        res_object <- clust_object[[i]]

        n_runs <- sum(sapply(res_object[[names(res_object)[1]]]$clusters, function(x) {
            sum(sapply(x$partitions, function(y) {
                y$freq
            }))
        }))

        list_appereances <- lapply(res_object, function(x) {
            lapply(x$clusters, function(y) {
                y$partitions[[1]]$freq
            })
        })
        temp_appereances <- reshape2::melt(list_appereances)
        colnames(temp_appereances) <- c(
            "freq_partition",
            "number_clusters",
            "resolution_value"
        )

        temp_appereances[["freq_k"]] <- unlist(lapply(res_object, function(x) {
            lapply(x$clusters, function(y) {
                sum(sapply(y$partitions, function(z) {
                    z$freq
                }))
            })
        }))

        temp_appereances[["configuration"]] <- rep(res_object_names[i], nrow(temp_appereances))
        temp_appereances$freq_partition <- temp_appereances$freq_partition / temp_appereances$freq_k
        temp_appereances$freq_k <- temp_appereances$freq_k / n_runs
        temp_appereances$ecc <- unlist(lapply(res_object, function(x) {
            sapply(x$clusters, function(k) {
                summary_function(k$ecc)
            })
        }))

        if (i == 1) {
            appearances_df <- temp_appereances
        } else {
            appearances_df <- rbind(appearances_df, temp_appereances)
        }
    }

    appearances_df[["configuration"]] <- factor(appearances_df[["configuration"]])
    appearances_df[["number_clusters"]] <- factor(as.numeric(appearances_df[["number_clusters"]]))

    main_plot <- ggplot2::ggplot(
        appearances_df,
        ggplot2::aes(
            y = .data$number_clusters,
            x = .data$resolution_value,
            size = .data$freq_partition,
            fill = .data[[colour_information]],
            shape = .data$configuration,
            group = .data$configuration
        )
    ) +
        ggplot2::geom_hline(
            yintercept = unique(appearances_df$number_clusters),
            linetype = "dashed",
            color = "#e3e3e3"
        ) +
        ggplot2::geom_vline(
            xintercept = unique(appearances_df$resolution_value),
            linetype = "dashed",
            color = "#e3e3e3"
        ) +
        ggplot2::geom_point(position = ggplot2::position_dodge(width = dodge_width)) +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_viridis_c(guide = "colorbar") +
        ggplot2::scale_shape_manual(
            name = "Clustering method",
            values = 21:24,
            guide = "legend"
        ) +
        ggplot2::labs(
            x = "resolution",
            y = "k"
        ) +
        ggplot2::scale_size_continuous(range = pt_size_range, guide = "legend") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
        )) +
        ggplot2::guides(
            shape = ggplot2::guide_legend(override.aes = list(size = max(pt_size_range)))
        )

    return(main_plot)
}

#' Relationship Between the Number of Clusters and the Number of Unique Partitions
#'
#' @description For each configuration provided in clust_object, display how
#' many different partitions with the same number of clusters can be obtained
#' by changing the seed.
#'
#' @param clust_object An object returned by the
#' `assess_clustering_stability` method.
#' @param colour_information String that specifies the information type that will be
#' illustrated using gradient colour: either `freq_part` for the frequency of the
#' most common partition or `ecc` for the Element-Centric Consistency
#' of the partitions obtained when the the number of clusters is fixed. Defaults to `ecc`.
#' @param dodge_width Used for adjusting the distance between the boxplots representing
#' a clustering method. Defaults to `0.3`.
#' @param pt_size_range Indicates the minimum and the maximum size a point on the plot can have.
#' Defaults to `c(1.5, 4)`.
#' @param summary_function The function that will be used to summarize the
#' distribution of the ECC values obtained for each number of clusters. Defaults
#' to `median`.
#' @param y_step The step used for the y-axis. Defaults to `5`.
#'
#' @return A ggplot2 object. The color gradient suggests the frequency of the most
#' common partition relative to the total number of appearances of that specific
#' number of clusters or the Element-Centric Consistency of the partitions. The size
#' illustrates the frequency of the partitions with *k* clusters relative to the
#' total number of partitions. The shape of the points indicates the clustering method.
#' @export
#'
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_embedding <- matrix(runif(100 * 30), nrow = 100)
#' rownames(pca_embedding) <- paste0("cell_", seq_len(nrow(pca_embedding)))
#' colnames(pca_embedding) <- paste0("PC_", 1:30)
#'
#'
#' adj_matrix <- getNNmatrix(
#'     RANN::nn2(pca_embedding, k = 10)$nn.idx,
#'     10,
#'     0,
#'     -1
#' )$nn
#' rownames(adj_matrix) <- paste0("cell_", seq_len(nrow(adj_matrix)))
#' colnames(adj_matrix) <- paste0("cell_", seq_len(ncol(adj_matrix)))
#'
#' # alternatively, the adj_matrix can be calculated
#' # using the `Seurat::FindNeighbors` function.
#'
#' clust_diff_obj <- assess_clustering_stability(
#'     graph_adjacency_matrix = adj_matrix,
#'     resolution = c(0.5, 1),
#'     n_repetitions = 10,
#'     clustering_algorithm = 1:2,
#'     verbose = FALSE
#' )
#' plot_k_n_partitions(clust_diff_obj)
plot_k_n_partitions <- function(clust_object,
                                colour_information = c("ecc", "freq_part"),
                                dodge_width = 0.3,
                                pt_size_range = c(1.5, 4),
                                summary_function = stats::median,
                                y_step = 5) {
    if (length(colour_information) > 1) {
        colour_information <- colour_information[1]
    }

    if ("split_by_k" %in% names(clust_object)) {
        clust_object <- clust_object$split_by_k
    }

    if (!(colour_information %in% c("ecc", "freq_part"))) {
        stop("colour_information can be either `ecc` or `freq_part`")
    }

    # use the names of the fields from the list
    object_names <- names(clust_object)

    max_n_part <- 0

    # creates a dataframe that contains, for each configuration
    # the number of different partitions with a given number of clusters
    for (i in seq_along(clust_object)) {
        partition_object <- clust_object[[i]]

        unique_parts_temp <- reshape2::melt(lapply(partition_object, function(x) {
            length(x$partitions)
        }))
        colnames(unique_parts_temp) <- c("n.partitions", "n.clusters")

        unique_parts_temp[["configuration"]] <- rep(object_names[i], nrow(unique_parts_temp))

        unique_parts_temp[["first.occ"]] <- as.numeric(lapply(partition_object, function(x) {
            max(sapply(x$partitions, function(y) {
                y$freq
            }))
        }))
        unique_parts_temp[["total.occ"]] <- as.numeric(lapply(partition_object, function(x) {
            sum(sapply(x$partitions, function(y) {
                y$freq
            }))
        }))
        unique_parts_temp[["freq_part"]] <- unique_parts_temp$first.occ / unique_parts_temp$total.occ
        unique_parts_temp[["ecc"]] <- sapply(partition_object, function(x) {
            summary_function(x$ecc)
        })
        overall_total_occ <- sum(unique_parts_temp$total.occ)
        unique_parts_temp[["frequency_k"]] <- unique_parts_temp$total.occ / overall_total_occ

        max_n_part <- max(c(max(
            unique_parts_temp$n.partitions
        ), max_n_part))

        if (i == 1) {
            unique_parts <- unique_parts_temp
        } else {
            unique_parts <- rbind(unique_parts, unique_parts_temp)
        }
    }

    unique_parts$configuration <- factor(unique_parts$configuration)
    unique_parts$n.clusters <- factor(unique_parts$n.clusters,
        levels = stringr::str_sort(unique(
            unique_parts$n.clusters
        ), numeric = TRUE)
    )

    main_plot <- ggplot2::ggplot(
        unique_parts,
        ggplot2::aes(
            x = .data$n.clusters,
            y = .data$n.partitions,
            shape = .data$configuration,
            size = .data$frequency_k,
            fill = .data[[colour_information]],
            group = .data$configuration
        )
    ) +
        ggplot2::scale_y_continuous(breaks = seq(
            from = 0,
            to = max_n_part,
            by = y_step
        )) +
        ggplot2::scale_size_continuous(range = pt_size_range, guide = "legend") +
        ggplot2::geom_hline(
            yintercept = seq(from = 0, to = max_n_part, by = y_step),
            linetype = "dashed",
            color = "#C3C3d3"
        ) +
        ggplot2::geom_vline(
            xintercept = unique(unique_parts$n.clusters),
            linetype = "dashed",
            color = "#c3c3c3d3"
        ) +
        ggplot2::geom_point(position = ggplot2::position_dodge(dodge_width)) +
        ggplot2::theme_classic() +
        ggplot2::scale_shape_manual(name = "Clustering method", values = 21:24, guide = "legend") +
        ggplot2::scale_fill_viridis_c(
            name = ifelse(colour_information == "ecc",
                "ECC",
                "partition\nfrequency"
            ),
            guide = "colorbar"
        ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("k") +
        ggplot2::ylab("# partitions") +
        ggplot2::guides(
            shape = ggplot2::guide_legend(override.aes = list(size = max(pt_size_range)))
        )

    return(main_plot)
}
