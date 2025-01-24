#### connected components ####
get_nn_conn_comps_umap <- function(embedding,
                                   n_neigh_sequence,
                                   n_repetitions = 100,
                                   seed_sequence = NULL,
                                   umap_arguments = list()) {
    ncores <- foreach::getDoParWorkers()

    umap_arguments[["n_threads"]] <- 1
    umap_arguments[["n_sgd_threads"]] <- 1

    ncells <- nrow(embedding)
    if (!("n_neighbors" %in% umap_arguments)) {
        umap_arguments[["n_neighbors"]] <- 15 # uwot's default
    }

    umap_arguments[["n_neighbors"]] <- min(umap_arguments[["n_neighbors"]], ncells - 1)

    nn_conn_comps_list <- list()
    if (ncores > 1 && is_package_installed("SharedObject")) {
        shared_embedding <- SharedObject::share(embedding)
    } else {
        shared_embedding <- embedding
    }

    if (length(n_neigh_sequence) > 1) {
        n_neighs <- c(n_neigh_sequence[1], sapply(2:length(n_neigh_sequence), function(i) {
            n_neigh_sequence[i] - n_neigh_sequence[i - 1]
        }))
        start_k <- n_neigh_sequence - n_neighs
    } else {
        start_k <- 0
        n_neighs <- n_neigh_sequence
    }

    all_vars <- ls()
    # the variables needed in each PSOCK process
    needed_vars <- c(
        "shared_embedding",
        "n_neigh_sequence",
        "graph_reduction_type",
        "start_k",
        "n_neighs",
        "umap_arguments"
    )
    seed <- NA
    # getNNmatrix <- getNNmatrix

    # send the name of the dim reduction arguments
    nn_conn_comps_list_temp <- foreach::foreach(
        seed = seed_sequence,
        .inorder = FALSE,
        .noexport = all_vars[!(all_vars %in% needed_vars)],
        .export = c("getNNmatrix"),
        .packages = c("ClustAssess")
    ) %dopar% {
        # perform the UMAP dimensionality reduction
        set.seed(seed)
        umap_embedding <- do.call(uwot::umap, c(list(X = shared_embedding), umap_arguments))
        colnames(umap_embedding) <- paste0("UMAP_", seq_len(ncol(umap_embedding)))
        rownames(umap_embedding) <- rownames(shared_embedding)

        nn2_res <- RANN::nn2(
            umap_embedding,
            k = max(n_neigh_sequence)
        )$nn.idx

        # for each neighbour, return the number of connected components
        # of the generated graph
        # TODO maybe add the evolution of the conn comps by varying the pruning param
        g <- NULL
        n_comps <- rep(0, length(n_neigh_sequence))
        for (neigh_index in seq_along(n_neigh_sequence)) {
            current_g <- igraph::graph_from_adjacency_matrix(
                getNNmatrix(nn2_res, n_neighs[neigh_index], start_k[neigh_index], -1)$nn
            )

            if (neigh_index == 1) {
                g <- current_g
            } else {
                g <- igraph::union(g, current_g)
            }

            n_comps[neigh_index] <- igraph::clusters(g)$no
        }

        n_comps
    }

    if (ncores > 1 && is_package_installed("SharedObject")) {
        shared_embedding <- SharedObject::unshare(embedding)
    }

    # store the results obtained for each number of neighbours in different lists
    for (i in seq_along(n_neigh_sequence)) {
        n_neigh <- n_neigh_sequence[i]
        nn_conn_comps_list[[as.character(n_neigh)]] <- sapply(nn_conn_comps_list_temp, function(x) {
            x[i]
        })
    }

    nn_conn_comps_list <- list(nn_conn_comps_list)
    names(nn_conn_comps_list) <- "UMAP"

    nn_conn_comps_list
}

get_nn_conn_comps_pca <- function(embedding,
                                  n_neigh_sequence) {
    ncores <- foreach::getDoParWorkers()
    nn_conn_comps_list <- list()


    nn2_res <- RANN::nn2(
        embedding,
        k = max(n_neigh_sequence)
    )$nn.idx

    if (ncores > 1 && is_package_installed("SharedObject")) {
        shared_nn2_res <- SharedObject::share(nn2_res)
    } else {
        shared_nn2_res <- nn2_res
    }

    if (length(n_neigh_sequence) > 1) {
        n_neighs <- c(n_neigh_sequence[1], sapply(2:length(n_neigh_sequence), function(i) {
            n_neigh_sequence[i] - n_neigh_sequence[i - 1]
        }))
        start_k <- n_neigh_sequence - n_neighs
    } else {
        start_k <- 0
        n_neighs <- n_neigh_sequence
    }

    all_vars <- ls()
    # the variables needed in each PSOCK process
    needed_vars <- c("shared_nn2_res", "start_k", "n_neighs")

    # send the name of the dim reduction arguments
    g_list <- foreach::foreach(
        i = seq_along(n_neigh_sequence),
        .inorder = TRUE,
        .noexport = all_vars[!(all_vars %in% needed_vars)],
        .export = c("getNNmatrix"),
        .packages = c("ClustAssess")
    ) %dopar% {
        # for each neighbour, return the number of connected components
        # of the generated graph
        igraph::graph_from_adjacency_matrix(
            getNNmatrix(shared_nn2_res, n_neighs[i], start_k[i], -1)$nn
        )
    }

    if (ncores > 1 && is_package_installed("SharedObject")) {
        shared_nn2_res <- SharedObject::unshare(nn2_res)
    }

    g <- NULL
    # store obtained results for each number of neighbours in different lists
    for (i in seq_along(n_neigh_sequence)) {
        n_neigh <- n_neigh_sequence[i]
        if (is.null(g)) {
            g <- g_list[[i]]
        } else {
            g <- igraph::union(g, g_list[[i]])
        }
        nn_conn_comps_list[[as.character(n_neigh)]] <- igraph::clusters(g)$no
    }

    nn_conn_comps_list <- list(nn_conn_comps_list)
    names(nn_conn_comps_list) <- "PCA"

    nn_conn_comps_list
}
#' Relationship Between Nearest Neighbours and Connected Components
#'
#' @description One of the steps in the clustering pipeline is building a
#' k-nearest neighbour graph on a reduced-space embedding. This method assesses
#' the relationship between different number of nearest
#' neighbours and the connectivity of the graph. In the context of graph clustering,
#' the number of connected components can be used as a
#' lower bound for the number of clusters. The calculations are performed multiple
#' times by changing the seed at each repetition.
#'
#' @param embedding A matrix associated with a PCA embedding. Embeddings from
#' other dimensionality reduction techniques (such as LSI) can be used.
#' @param n_neigh_sequence A sequence of the number of nearest neighbours.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user. Defaults to `100``.
#' @param seed_sequence A custom seed sequence; if the value is NULL, the
#' sequence will be built starting from 1 with a step of 100.
#' @param include_umap A boolean value indicating whether to calculate the number
#' of connected components for the UMAP embedding. Defaults to `FALSE`.
#' @param umap_arguments Additional arguments passed to the the `uwot::umap` method.
#'
#'
#' @return A list having one field associated with a number of nearest neighbours.
#' Each value contains an array of the number of connected components
#' obtained on the specified number of repetitions.
#'
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_emb <- matrix(runif(100 * 30), nrow = 100, byrow = TRUE)
#' rownames(pca_emb) <- as.character(1:100)
#' colnames(pca_emb) <- paste0("PCA_", 1:30)
#'
#' nn_conn_comps_obj <- get_nn_conn_comps(
#'     embedding = pca_emb,
#'     n_neigh_sequence = c(2, 5),
#'     n_repetitions = 3,
#'     # arguments that are passed to the uwot function
#'     umap_arguments = list(
#'         min_dist = 0.3,
#'         metric = "cosine"
#'     )
#' )
#' plot_connected_comps_evolution(nn_conn_comps_obj)
get_nn_conn_comps <- function(embedding,
                              n_neigh_sequence,
                              n_repetitions = 100,
                              seed_sequence = NULL,
                              include_umap = FALSE,
                              umap_arguments = list()) {
    # check parameters
    if (!is.numeric(n_neigh_sequence)) {
        stop("n_neigh_sequence parameter should be numeric")
    }
    # convert number of neighbours to integers
    n_neigh_sequence <- sort(as.integer(n_neigh_sequence))

    if (!is.numeric(n_repetitions) || length(n_repetitions) > 1) {
        stop("n_repetitions parameter should be numeric")
    }
    # convert n_repetitions to integers
    n_repetitions <- as.integer(n_repetitions)

    if (!is.matrix(embedding) && !inherits(embedding, "Matrix")) {
        stop("embedding parameter should be a matrix")
    }

    ncells <- nrow(embedding)
    n_neigh_sequence <- n_neigh_sequence[which(n_neigh_sequence < ncells)]

    if (length(n_neigh_sequence) == 0) {
        warning(glue::glue("The provided values for the `n_neigh_sequence` are greater than the number of cells ({ncells}). For the downstream analysis, we will set `n_neigh` to {ncells}."))
        n_neigh_sequence <- c(ncells - 1)
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

    nn_comps_result <- get_nn_conn_comps_pca(
        embedding = embedding,
        n_neigh_sequence = n_neigh_sequence
    )

    if (!include_umap) {
        return(nn_comps_result)
    }

    return(c(
        nn_comps_result,
        get_nn_conn_comps_umap(
            embedding = embedding,
            n_neigh_sequence = n_neigh_sequence,
            n_repetitions = n_repetitions,
            seed_sequence = seed_sequence,
            umap_arguments = umap_arguments
        )
    ))
}

#' Relationship Between Number of Nearest Neighbours and Graph Connectivity
#'
#' @description Display the distribution of the number connected components
#' obtained for each number of neighbours across random seeds.
#'
#' @param nn_conn_comps_object An object or a concatenation of objects returned
#' by the `get_nn_conn_comps` method.
#'
#'
#' @return A ggplot2 object with boxplots for the connected component distributions.
#' @export
#'
#' @note The number of connected components is displayed on a logarithmic scale.
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_emb <- matrix(runif(100 * 30), nrow = 100, byrow = TRUE)
#' rownames(pca_emb) <- as.character(1:100)
#' colnames(pca_emb) <- paste0("PCA_", 1:30)
#'
#' nn_conn_comps_obj <- get_nn_conn_comps(
#'     embedding = pca_emb,
#'     n_neigh_sequence = c(2, 5),
#'     n_repetitions = 3,
#'     # arguments that are passed to the uwot function
#'     umap_arguments = list(
#'         min_dist = 0.3,
#'         metric = "cosine"
#'     )
#' )
#' plot_connected_comps_evolution(nn_conn_comps_obj)
plot_connected_comps_evolution <- function(nn_conn_comps_object) {
    for (n_neighbours in names(nn_conn_comps_object$PCA)) {
        nn_conn_comps_object$PCA[[n_neighbours]] <- as.integer(nn_conn_comps_object$PCA[[n_neighbours]])
    }

    for (n_neighbours in names(nn_conn_comps_object$UMAP)) {
        nn_conn_comps_object$UMAP[[n_neighbours]] <- as.integer(nn_conn_comps_object$UMAP[[n_neighbours]])
    }

    final_comps_df <- reshape2::melt(nn_conn_comps_object)
    colnames(final_comps_df) <- c("n_comps", "n_neigh", "config_name")

    final_comps_df$config_name <- factor(final_comps_df$config_name)
    final_comps_df$n_neigh <- factor(final_comps_df$n_neigh)
    final_comps_df$n_neigh <- factor(final_comps_df$n_neigh, levels(final_comps_df$n_neigh)[stringr::str_order(levels(final_comps_df$n_neigh), numeric = TRUE)])

    max_y <- max(final_comps_df$n_comps)
    max_y <- max_y + ifelse(max_y %% 5, 5 - max_y %% 5, 0)
    min_y <- min(final_comps_df$n_comps)
    min_y <- min_y - min_y %% 5

    # generate additional y ticks that would simulate the exponential tendency
    chosen_breaks <- generate_breaks(min_y, max_y)

    ggplot2::ggplot(
        final_comps_df,
        ggplot2::aes(
            x = .data$n_neigh,
            y = .data$n_comps,
            fill = .data$config_name
        )
    ) +
        ggplot2::geom_hline(
            yintercept = chosen_breaks,
            linetype = "dashed",
            color = "#bcc4be"
        ) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_y_continuous(breaks = chosen_breaks, trans = "log10") +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90,
            vjust = 1,
            hjust = 1
        )) +
        ggplot2::labs(
            x = "# of nearest neighbours",
            y = "# of connected components",
            fill = "configuration"
        ) +
        ggplot2::ggtitle("Distribution of the number of connected components")
}

#### number of neigh, graph type importance ####

get_n_strong_components <- function(nn_matrix,
                                    n_neigh,
                                    prune_snn) {
    g <- igraph::graph_from_adjacency_matrix(
        computeSNN(nn_matrix, n_neigh, prune_snn),
        mode = "undirected",
        weighted = TRUE
    )

    igraph::clusters(g)$no
}

#' Calculate the highest pruning parameter for the SNN graph given Embedding
#'
#' @description Given an embedding, the function calculates the highest pruning
#' parameter for the SNN graph that preserves the connectivity of the graph.
#'
#' @param embedding A matrix associated with a PCA embedding. Embeddings from
#' other dimensionality reduction techniques (such as LSI) can be used.
#' @param n_neigh The number of nearest neighbours.
#'
#' @return The value of the highest pruning parameter.
#' @export
#'
#' @note Given the way the SNN graph is built, the possible values for the pruning
#' parameter are limited and can be determined by the formula `i / (2 * n_neigh - i)`,
#' where `i` is a number of nearest neighbours between 0 and `n_neigh`.
#'
#' @examples
#' set.seed(2024)
#' # create an artificial pca embedding
#' pca_embedding <- matrix(
#'     c(runif(100 * 10), runif(100 * 10, min = 3, max = 4)),
#'     nrow = 200, byrow = TRUE
#' )
#' rownames(pca_embedding) <- as.character(1:200)
#' colnames(pca_embedding) <- paste("PC", 1:10)
#'
#' get_highest_prune_param_embedding(pca_embedding, 5)
get_highest_prune_param_embedding <- function(embedding,
                                              n_neigh) {
    nn_matrix <- getNNmatrix(RANN::nn2(
        embedding,
        k = n_neigh
    )$nn.idx, n_neigh, 0, -1)$nn

    g <- igraph::graph_from_adjacency_matrix(
        nn_matrix,
        mode = "directed"
    )

    target_n_conn_comps <- igraph::clusters(g)$no

    possible_values <- sapply(0:n_neigh, function(i) {
        i / (2 * n_neigh - i)
    })

    start_n <- 1
    stop_n <- length(possible_values)

    while (start_n <= stop_n) {
        middle <- as.integer((stop_n + start_n) / 2)
        current_n_conn_comps <- get_n_strong_components(
            nn_matrix,
            n_neigh,
            possible_values[middle]
        )

        if (current_n_conn_comps > target_n_conn_comps) {
            stop_n <- middle

            if (start_n == stop_n) {
                start_n <- middle - 1
            }
        } else {
            start_n <- middle + 1

            if (start_n == stop_n) {
                break
            }
        }
    }

    return(possible_values[middle])
}

#' Calculate the highest pruning parameter for the SNN graph given NN matrix
#'
#' @description Given a NN adjacency matrix, the function calculates the highest pruning
#' parameter for the SNN graph that preserves the connectivity of the graph.
#'
#' @param nn_matrix The adjacency matrix of the nearest neighbour graph.
#' @param n_neigh The number of nearest neighbours.
#'
#' @return A list with the following fields:
#' - `prune_value`: The value of the highest pruning parameter.
#' - `adj_matrix`: The adjacency matrix of the SNN graph after pruning.
#' @export
#'
#' @note Given the way the SNN graph is built, the possible values for the pruning
#' parameter are limited and can be determined by the formula `i / (2 * n_neigh - i)`,
#' where `i` is a number of nearest neighbours between 0 and `n_neigh`.
#'
#' @examples
#' set.seed(2024)
#' # create an artificial pca embedding
#' pca_embedding <- matrix(
#'     c(runif(100 * 10), runif(100 * 10, min = 3, max = 4)),
#'     nrow = 200, byrow = TRUE
#' )
#' rownames(pca_embedding) <- as.character(1:200)
#' colnames(pca_embedding) <- paste("PC", 1:10)
#'
#' # calculate the nn adjacency matrix
#' nn_matrix <- getNNmatrix(
#'     RANN::nn2(pca_embedding, k = 5)$nn.idx,
#'     5,
#'     0,
#'     -1
#' )$nn
#'
#' get_highest_prune_param(nn_matrix, 5)$prune_value
get_highest_prune_param <- function(nn_matrix,
                                    n_neigh) {
    nn_matrix <- computeSNN(nn_matrix, n_neigh, 0)
    g <- igraph::graph_from_adjacency_matrix(
        nn_matrix,
        mode = "undirected",
        weighted = TRUE
    )

    target_n_conn_comps <- igraph::clusters(g)$no

    possible_values <- sapply(0:n_neigh, function(i) {
        i / (2 * n_neigh - i)
    })

    start_n <- 1
    stop_n <- length(possible_values)
    prev_g <- g

    while (start_n <= stop_n) {
        middle <- as.integer((stop_n + start_n) / 2)
        current_g <- igraph::delete_edges(prev_g, which(igraph::E(prev_g)$weight <= possible_values[middle]))
        gc()

        current_n_conn_comps <- igraph::clusters(current_g)$no

        if (current_n_conn_comps > target_n_conn_comps) {
            stop_n <- middle

            if (start_n == stop_n) {
                start_n <- middle - 1
            }
        } else {
            start_n <- middle + 1
            prev_g <- current_g

            if (start_n == stop_n) {
                break
            }
        }
    }

    return(list(
        prune_value = possible_values[middle],
        adj_matrix = pruneSNN(nn_matrix, possible_values[middle])
    ))
}

assess_nn_stability_pca <- function(embedding,
                                    n_neigh_sequence,
                                    n_repetitions = 100,
                                    seed_sequence = NULL,
                                    ecs_thresh = 1,
                                    graph_type = 2,
                                    prune_value = -1,
                                    clustering_algorithm = 1,
                                    clustering_arguments = list()) {
    ncores <- foreach::getDoParWorkers()
    cell_names <- rownames(embedding)
    partitions_list <- list()

    if (graph_type != 0) {
        partitions_list[[paste("PCA", "snn", sep = "_")]] <- list()
    }

    if (graph_type != 1) {
        partitions_list[[paste("PCA", "nn", sep = "_")]] <- list()
    }

    nn2_res <- RANN::nn2(
        embedding,
        k = max(n_neigh_sequence)
    )$nn.idx

    if (ncores > 1 && is_package_installed("SharedObject")) {
        shared_nn2_res <- SharedObject::share(nn2_res)
    } else {
        shared_nn2_res <- nn2_res
    }

    all_vars <- ls()
    needed_vars <- c(
        "cell_names",
        "shared_nn2_res",
        "graph_type",
        "prune_value"
    )
    neigh_matrices <- foreach::foreach(
        n_neigh = n_neigh_sequence,
        .inorder = TRUE,
        .noexport = all_vars[!(all_vars %in% needed_vars)],
        .export = c("getNNmatrix"),
        .packages = c("ClustAssess")
    ) %dopar% {
        neigh_matrix <- getNNmatrix(shared_nn2_res, n_neigh, 0, -1)
        rownames(neigh_matrix$nn) <- cell_names
        colnames(neigh_matrix$nn) <- cell_names

        if (graph_type > 0) {
            if (prune_value >= 0) {
                neigh_matrix$snn <- computeSNN(neigh_matrix$nn, n_neigh, prune_value)
            } else {
                neigh_matrix$snn <- get_highest_prune_param(
                    nn_matrix = neigh_matrix$nn,
                    n_neigh = n_neigh
                )$adj_matrix
            }
            rownames(neigh_matrix$snn) <- cell_names
            colnames(neigh_matrix$snn) <- cell_names
        }

        neigh_matrix
    }

    names(neigh_matrices) <- as.character(n_neigh_sequence)

    if (ncores > 1 && is_package_installed("SharedObject")) {
        shared_nn2_res <- SharedObject::unshare(nn2_res)
        rm(shared_nn2_res)
    }

    rm(nn2_res)
    gc()

    package_needed <- c()
    if (4 %in% clustering_algorithm) {
        package_needed <- c(package_needed, "leiden")
    }
    for (n_neigh in as.character(n_neigh_sequence)) {
        partitions_list[[paste("PCA", "snn", sep = "_")]][[n_neigh]] <- list()

        if (ncores > 1 && is_package_installed("SharedObject")) {
            shared_neigh_matrix <- SharedObject::share(neigh_matrices[[n_neigh]])
        } else {
            shared_neigh_matrix <- neigh_matrices[[n_neigh]]
        }

        gc()
        # the variables needed in the PSOCK processes
        seed <- NA
        needed_vars <- c(
            "graph_type",
            "algorithm",
            "shared_neigh_matrix",
            "clustering_arguments"
        )
        all_vars <- ls()

        partitions_list_temp <- foreach::foreach(
            seed = seed_sequence,
            .inorder = FALSE,
            .noexport = all_vars[!(all_vars %in% needed_vars)],
            .packages = package_needed
        ) %dopar% {
            # apply the clustering method on the graph specified by the variable `graph_type`
            if (graph_type != 1) {
                cluster_results_nn <- list(
                    mb = do.call(
                        clustering_functions,
                        c(
                            list(
                                object = shared_neigh_matrix$nn,
                                resolution = 0.8,
                                seed = seed
                            ),
                            clustering_arguments
                        )
                    ),
                    freq = 1,
                    seed = seed
                )

                if (graph_type == 0) {
                    return(cluster_results_nn)
                }
            }

            cluster_results_snn <- list(
                mb = do.call(
                    clustering_functions,
                    c(
                        list(
                            object = shared_neigh_matrix$snn,
                            resolution = 0.8,
                            seed = seed
                        ),
                        clustering_arguments
                    )
                ),
                freq = 1,
                seed = seed
            )

            if (graph_type == 1) {
                return(cluster_results_snn)
            }

            return(list(cluster_results_snn, cluster_results_nn))
        }

        neigh_matrices[[n_neigh]] <- NULL

        # merge the partitions that are considered similar by a given ecs threshold
        for (i in seq_along(partitions_list)) {
            partitions_list[[i]][[n_neigh]] <- merge_partitions(
                lapply(partitions_list_temp, function(x) {
                    x[[i]]
                }),
                ecs_thresh = ecs_thresh
            )
        }
    }

    # TODO raise a github issue for SharedObject: is it necesarry to unshare objects??
    # if (ncores > 1) {
    #     shared_neigh_matrix <- SharedObject::unshare(neigh_matrices)
    #     # rm(shared_neigh_matrix)
    #     gc()
    # }

    # create an object showing the number of clusters obtained for each number
    # of neighbours
    nn_object_n_clusters <- list()
    for (config_name in names(partitions_list)) {
        nn_object_n_clusters[[config_name]] <- list()
        for (n_neigh in names(partitions_list[[config_name]])) {
            nn_object_n_clusters[[config_name]][[n_neigh]] <- unlist(lapply(partitions_list[[config_name]][[n_neigh]]$partitions, function(x) {
                rep(length(unique(x$mb)), x$freq)
            }))
        }
    }

    # create an object that contains the ecc of the partitions obtained for each
    # number of neighbours
    nn_ecs_object <- lapply(partitions_list, function(config) {
        lapply(config, function(n_neigh) {
            n_neigh$ecc
            # weighted_element_consistency(
            #     lapply(n_neigh, function(x) {
            #         x$mb
            #     }),
            #     sapply(n_neigh, function(x) {
            #         x$freq
            #     })
            # )
        })
    })

    list(
        n_neigh_k_corresp = nn_object_n_clusters,
        n_neigh_ec_consistency = nn_ecs_object,
        n_different_partitions = lapply(partitions_list, function(config) {
            sapply(config, function(n_neigh) {
                length(n_neigh)
            })
        })
    )
}

assess_nn_stability_umap <- function(embedding,
                                     n_neigh_sequence,
                                     n_repetitions = 100,
                                     seed_sequence = NULL,
                                     ecs_thresh = 1,
                                     graph_type = 2,
                                     prune_value = -1,
                                     clustering_algorithm = 1,
                                     clustering_arguments = list(),
                                     umap_arguments = list()) {
    ncores <- foreach::getDoParWorkers()

    object_names <- c("snn", "nn")
    if (graph_type < 2) {
        object_names <- object_names[graph_type + 1]
    }

    ncells <- nrow(embedding)
    umap_arguments <- process_umap_arguments(umap_arguments, ncells)

    if (ncores > 1 && is_package_installed("SharedObject")) {
        shared_embedding <- SharedObject::share(embedding)
    } else {
        shared_embedding <- embedding
    }

    # the variables needed in the PSOCK processes
    needed_vars <- c(
        "shared_embedding",
        "n_neigh_sequence",
        "graph_reduction_type",
        "graph_type",
        "umap_arguments",
        "clustering_arguments",
        "prune_SNN",
        "prune_value"
    )

    seed <- NA
    all_vars <- ls()

    package_needed <- c("ClustAssess")
    if (4 %in% clustering_algorithm) {
        package_needed <- c(package_needed, "leiden")
    }

    seed_list <- foreach::foreach(
        seed = seed_sequence,
        .inorder = FALSE,
        .noexport = all_vars[!(all_vars %in% needed_vars)],
        .export = c("getNNmatrix"),
        .packages = package_needed
    ) %dopar% {
        # perform the dimensionality reduction
        set.seed(seed)
        row_names <- rownames(shared_embedding)
        embedding <- do.call(uwot::umap, c(list(X = shared_embedding), umap_arguments))
        colnames(embedding) <- paste0("UMAP_", seq_len(ncol(embedding)))
        rownames(embedding) <- row_names

        nn2_res <- RANN::nn2(
            embedding,
            k = max(n_neigh_sequence)
        )$nn.idx

        seed_result <- list()

        for (n_neigh in n_neigh_sequence) {
            neigh_matrix <- getNNmatrix(nn2_res, n_neigh, 0, -1)
            rownames(neigh_matrix$nn) <- row_names
            colnames(neigh_matrix$nn) <- row_names

            if (graph_type > 0) {
                if (prune_value >= 0) {
                    neigh_matrix$snn <- computeSNN(neigh_matrix$nn, n_neigh, prune_value)
                } else {
                    neigh_matrix$snn <- get_highest_prune_param(
                        nn_matrix = neigh_matrix$nn,
                        n_neigh = n_neigh
                    )$adj_matrix
                }
                rownames(neigh_matrix$snn) <- row_names
                colnames(neigh_matrix$snn) <- row_names
            }

            # apply the clustering method on the graph specified by the variable `graph_type`
            if (graph_type > 0) {
                seed_result[["snn"]][[as.character(n_neigh)]] <- list(
                    mb = do.call(
                        clustering_functions,
                        c(
                            list(
                                object = neigh_matrix$snn,
                                resolution = 0.8,
                                seed = seed
                            ),
                            clustering_arguments
                        )
                    ),
                    freq = 1,
                    seed = seed
                )
            }

            if (graph_type %% 2 == 0) {
                seed_result[["nn"]][[as.character(n_neigh)]] <- list(
                    mb = do.call(
                        clustering_functions,
                        c(
                            list(
                                object = neigh_matrix$nn,
                                resolution = 0.8,
                                seed = seed
                            ),
                            clustering_arguments
                        )
                    ),
                    freq = 1,
                    seed = seed
                )
            }
        }

        return(seed_result)
    }

    if (ncores > 1 && is_package_installed("SharedObject")) {
        shared_embedding <- SharedObject::unshare(embedding)
    }

    rm(shared_embedding)
    gc()

    partitions_list <- lapply(object_names, function(gtype) {
        l <- lapply(seq_along(n_neigh_sequence), function(k) {
            lapply(seq_len(n_repetitions), function(seed) {
                seed_list[[seed]][[gtype]][[k]]
            })
        })
        names(l) <- as.character(n_neigh_sequence)
        return(l)
    })

    names(partitions_list) <- paste0("UMAP_", object_names)

    # merge the partitions that are considered similar by a given ecs threshold
    for (i in seq_along(partitions_list)) {
        for (n_neigh in as.character(n_neigh_sequence)) {
            partitions_list[[i]][[n_neigh]] <- merge_partitions(
                partitions_list[[i]][[n_neigh]],
                ecs_thresh = ecs_thresh
            )
        }
    }

    # create an object showing the number of clusters obtained for each number of
    # neighbours
    nn_object_n_clusters <- list()
    for (config_name in names(partitions_list)) {
        nn_object_n_clusters[[config_name]] <- list()
        for (n_neigh in names(partitions_list[[config_name]])) {
            nn_object_n_clusters[[config_name]][[n_neigh]] <- unlist(lapply(partitions_list[[config_name]][[n_neigh]]$partitions, function(x) {
                rep(length(unique(x$mb)), x$freq)
            }))
        }
    }

    # create an object that contains the ecc of the partitions obtained for each
    # number of neighbours
    nn_ecs_object <- lapply(partitions_list, function(config) {
        lapply(config, function(n_neigh) {
            n_neigh$ecc
            # weighted_element_consistency(
            #     lapply(n_neigh, function(x) {
            #         x$mb
            #     }),
            #     sapply(n_neigh, function(x) {
            #         x$freq
            #     })
            # )
        })
    })

    list(
        n_neigh_k_corresp = nn_object_n_clusters,
        n_neigh_ec_consistency = nn_ecs_object,
        n_different_partitions = lapply(partitions_list, function(config) {
            sapply(config, function(n_neigh) {
                length(n_neigh)
            })
        })
    )
}

#' Assess the stability for Graph Building Parameters
#'
#' @description Evaluates clustering stability when changing the values of different
#' parameters involved in the graph building step,
#' namely the base embedding, the graph type and the number of neighbours.
#'

#' @param embedding A matrix associated with a PCA embedding. Embeddings from
#' other dimensionality reduction techniques (such as LSI) can be used.
#' @param n_neigh_sequence A sequence of the number of nearest neighbours.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence A custom seed sequence; if the value is NULL,
#' the sequence will be built starting from 1 with a step of 100.
#' @param graph_reduction_type The graph reduction type, denoting if the graph
#' should be built on either the PCA or the UMAP embedding.
#' @param ecs_thresh The ECS threshold used for merging similar clusterings.
#' @param graph_type Argument indicating whether the graph should be
#' unweighted (0), weighted (1) or both (2).
#' @param prune_value Argument indicating whether to prune the SNN graph. If the
#' value is 0, the graph won't be pruned. If the value is between 0 and 1, the
#' edges with weight under the pruning value will be removed. If the value is
#' -1, the highest pruning value will be calculated automatically and used.
#' @param clustering_algorithm An index indicating which community detection algorithm will
#' be used: Louvain (1), Louvain refined (2), SLM (3) or Leiden (4). More
#' details can be found in the Seurat's `FindClusters` function.
#' @param clustering_arguments A list of arguments that will be passed to the
#' clustering algorithm. See the `FindClusters` function in Seurat for more details.
#' @param umap_arguments Additional arguments passed to the the `uwot::umap` method.
#'
#' @return A list having three fields:
#'
#' * `n_neigh_k_corresp` - list containing the number of the clusters obtained by running
#' the pipeline multiple times with different seed, number of neighbours and graph type (weighted vs unweigted)
#' * `n_neigh_ec_consistency` - list containing the EC consistency of the partitions obtained
#' at multiple runs when changing the number of neighbours or the graph type
#' * `n_different_partitions` - the number of different partitions obtained by each
#' number of neighbours
#'
#' @export
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_emb <- matrix(runif(100 * 30), nrow = 100, byrow = TRUE)
#' rownames(pca_emb) <- as.character(1:100)
#' colnames(pca_emb) <- paste0("PC_", 1:30)
#'
#' nn_stability_obj <- assess_nn_stability(
#'     embedding = pca_emb,
#'     n_neigh_sequence = c(10, 15, 20),
#'     n_repetitions = 10,
#'     graph_reduction_type = "PCA",
#'     clustering_algorithm = 1
#' )
#' plot_n_neigh_ecs(nn_stability_obj)
assess_nn_stability <- function(embedding,
                                n_neigh_sequence,
                                n_repetitions = 100,
                                seed_sequence = NULL,
                                graph_reduction_type = "PCA",
                                ecs_thresh = 1,
                                graph_type = 2,
                                prune_value = -1,
                                clustering_algorithm = 1,
                                clustering_arguments = list(),
                                umap_arguments = list()) {
    # TODO vary by resolution
    # BUG check the performance
    # check parameters
    if (!is.numeric(n_neigh_sequence)) {
        stop("n_neigh_sequence parameter should be numeric")
    }
    # convert number of neighbours to integers
    n_neigh_sequence <- sort(as.integer(n_neigh_sequence))

    if (!is.numeric(n_repetitions) || length(n_repetitions) > 1) {
        stop("n_repetitions parameter should be numeric")
    }
    # convert n_repetitions to integers
    n_repetitions <- as.integer(n_repetitions)

    if (!is.numeric(ecs_thresh) || length(ecs_thresh) > 1) {
        stop("ecs_thresh parameter should be numeric")
    }

    if (!is.numeric(graph_type) || !(graph_type %in% 0:2)) {
        stop("graph_type should be a number between 0 and 2")
    }

    if (!is.matrix(embedding) && !inherits(embedding, "Matrix")) {
        stop("the embedding parameter should be a matrix")
    }

    if (!is.numeric(clustering_algorithm) || length(clustering_algorithm) > 1 || !(clustering_algorithm %in% 1:4)) {
        stop("algorithm should be a number between 1 and 4")
    }

    if (!(graph_reduction_type %in% c("PCA", "UMAP"))) {
        stop("graph_reduction_type parameter should take one of these values: 'PCA' or 'UMAP'")
    }

    ncells <- nrow(embedding)
    n_neigh_sequence <- n_neigh_sequence[which(n_neigh_sequence < ncells)]

    clustering_arguments <- process_clustering_arguments(clustering_arguments, clustering_algorithm)

    if (length(n_neigh_sequence) == 0) {
        warning(glue::glue("The provided values for the `n_neigh_sequence` are greater than the number of cells ({ncells}). For the downstream analysis, we will set `n_neigh` to {ncells - 1}."))
        n_neigh_sequence <- c(ncells - 1)
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

    if (graph_reduction_type == "PCA") {
        return(assess_nn_stability_pca(
            embedding = embedding,
            n_neigh_sequence = n_neigh_sequence,
            n_repetitions = n_repetitions,
            seed_sequence = seed_sequence,
            ecs_thresh = ecs_thresh,
            graph_type = graph_type,
            prune_value = prune_value,
            clustering_algorithm = clustering_algorithm,
            clustering_arguments = clustering_arguments
        ))
    }

    return(assess_nn_stability_umap(
        embedding = embedding,
        n_neigh_sequence = n_neigh_sequence,
        n_repetitions = n_repetitions,
        seed_sequence = seed_sequence,
        ecs_thresh = ecs_thresh,
        graph_type = graph_type,
        prune_value = prune_value,
        clustering_algorithm = clustering_algorithm,
        clustering_arguments = clustering_arguments,
        umap_arguments = umap_arguments
    ))
}

# 2. neighbours <-> number of clusters association

#' Relationship Between Number of Nearest Neighbours and Number of Clusters
#'
#' @description Display the distribution of the
#' number of clusters obtained for each number of neighbours across random seeds.
#'
#' @param nn_object_n_clusters An object or a concatenation of objects returned by the
#' `get_nn_importance` method.
#'
#'
#' @return A ggplot2 object with the distributions displayed as boxplots.
#' @export
#'
#' @note The number of clusters is displayed on a logarithmic scale.
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_emb <- matrix(runif(100 * 30), nrow = 100, byrow = TRUE)
#' rownames(pca_emb) <- as.character(1:100)
#' colnames(pca_emb) <- paste0("PC_", 1:30)
#'
#' nn_stability_obj <- assess_nn_stability(
#'     embedding = pca_emb,
#'     n_neigh_sequence = c(10, 15, 20),
#'     n_repetitions = 10,
#'     graph_reduction_type = "PCA",
#'     clustering_algorithm = 1
#' )
#' plot_n_neigh_k_correspondence(nn_stability_obj)
plot_n_neigh_k_correspondence <- function(nn_object_n_clusters) {
    nn_object_n_clusters <- nn_object_n_clusters$n_neigh_k_corresp
    for (config_name in names(nn_object_n_clusters)) {
        for (n_neighbours in names(nn_object_n_clusters[[config_name]])) {
            nn_object_n_clusters[[config_name]][[n_neighbours]] <- as.integer(nn_object_n_clusters[[config_name]][[n_neighbours]])
        }
    }
    melted_obj <- reshape2::melt(nn_object_n_clusters)
    colnames(melted_obj) <- c("k", "n_neigh", "config_name")

    melted_obj$n_neigh <- factor(melted_obj$n_neigh)
    melted_obj$n_neigh <- factor(melted_obj$n_neigh, levels(melted_obj$n_neigh)[stringr::str_order(levels(melted_obj$n_neigh), numeric = TRUE)])

    max_y <- max(melted_obj$k)
    max_y <- max_y + ifelse(max_y %% 5, 5 - max_y %% 5, 0)
    min_y <- min(melted_obj$k)
    min_y <- min_y - min_y %% 5

    # generate additional y ticks that would simulate the exponential tendency
    chosen_breaks <- generate_breaks(min_y, max_y)

    ggplot2::ggplot(
        melted_obj,
        ggplot2::aes(
            x = .data$n_neigh,
            y = .data$k,
            fill = .data$config_name
        )
    ) +
        ggplot2::geom_hline(
            yintercept = chosen_breaks,
            linetype = "dashed",
            color = "#bcc4be"
        ) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(breaks = chosen_breaks, trans = "log10") +
        ggplot2::labs(
            x = "# of nearest neighbours",
            fill = "configuration"
        ) +
        ggplot2::ggtitle("Distribution of k across different seeds")
}

# 3. ECS distribution across different seeds for different number of neighbours

#' Graph construction parameters - ECC facet
#'
#' @description Display, for all configurations consisting in different number
#' of neighbours, graph types and base embeddings, the EC Consistency of the partitions
#' obtained over multiple runs on an UMAP embedding.
#'
#' @param nn_ecs_object An object or a concatenation of objects returned by the
#' `get_nn_importance` method.
#' @param boxplot_width Used for adjusting the width of the boxplots; the value will
#' be passed to the `width` argument of the `ggplot2::geom_boxplot` method.
#'
#' @return A ggplot2 object.
#' @export
#'
#'
#' @examples
#' set.seed(2024)
#' # create an artificial PCA embedding
#' pca_emb <- matrix(runif(100 * 30), nrow = 100, byrow = TRUE)
#' rownames(pca_emb) <- as.character(1:100)
#' colnames(pca_emb) <- paste0("PC_", 1:30)
#'
#' nn_stability_obj <- assess_nn_stability(
#'     embedding = pca_emb,
#'     n_neigh_sequence = c(10, 15, 20),
#'     n_repetitions = 10,
#'     graph_reduction_type = "PCA",
#'     clustering_algorithm = 1
#' )
#' plot_n_neigh_ecs(nn_stability_obj)
plot_n_neigh_ecs <- function(nn_ecs_object,
                             boxplot_width = 0.5) {
    nn_ecs_object <- nn_ecs_object$n_neigh_ec_consistency
    for (config_name in names(nn_ecs_object)) {
        for (n_neighbours in names(nn_ecs_object[[config_name]])) {
            nn_ecs_object[[config_name]][[n_neighbours]] <- as.numeric(nn_ecs_object[[config_name]][[n_neighbours]])
        }
    }
    melted_obj <- reshape2::melt(nn_ecs_object)
    colnames(melted_obj) <- c("ECC", "n_neigh", "config_name")

    melted_obj$n_neigh <- factor(melted_obj$n_neigh)
    melted_obj$n_neigh <- factor(melted_obj$n_neigh, levels(melted_obj$n_neigh)[stringr::str_order(levels(melted_obj$n_neigh), numeric = TRUE)])

    ggplot2::ggplot(
        melted_obj,
        ggplot2::aes(
            x = .data$n_neigh,
            y = .data$ECC,
            fill = .data$config_name
        )
    ) +
        ggplot2::geom_boxplot(width = boxplot_width) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = "# of nearest neighbours",
            y = "EC consistency",
            fill = "configuration"
        ) +
        ggplot2::ggtitle("Distribution of ECC across different seeds for different # neighbours")
}
