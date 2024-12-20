# TODO check user interrupt
# TODO add logging options to files - maybe verbose lvels
#' @importFrom foreach %dopar%
NULL

ranking_functions <- list(
    "median" = stats::median,
    "max" = max,
    "top_qt" = function(arr) {
        stats::fivenum(arr)[4]
    },
    "top_qt_max" = function(arr) {
        arr_summary <- stats::fivenum(arr)
        arr_summary[4] + arr_summary[5]
    },
    "iqr" = function(arr) {
        arr_summary <- stats::fivenum(arr)
        arr_summary[2] - arr_summary[4]
    },
    "iqr_median" = function(arr) {
        arr_summary <- stats::fivenum(arr)
        arr_summary[3] - (arr_summary[4] - arr_summary[2])
    },
    "iqr_median_coeff" = function(arr) {
        arr_summary <- stats::fivenum(arr)
        arr_summary[3] / (arr_summary[4] - arr_summary[2])
    },
    "mean" = mean
)

# wrapper of the Seurat's `FindClusters` method, that returns
# only the membership vector
seurat_clustering <- function(object, resolution, seed, algorithm = 3, num_start = 10, num_iter = 10, ...) {
    cluster_result <- Seurat::FindClusters(
        object,
        resolution = resolution,
        random.seed = seed,
        algorithm = algorithm,
        n.start = num_start,
        n.iter = num_iter,
        ...
    )
    as.integer(cluster_result[[colnames(cluster_result)[1]]])
}

leiden_clustering <- function(g,
                              resolution,
                              seed,
                              initial_membership = NULL,
                              num_iter = 10,
                              num_starts = 1,
                              ...) {
    best_mb <- NULL
    best_quality <- NULL
    # g <- igraph::graph_from_adjacency_matrix(
    #     adjmatrix = object,
    #     mode = graph_mode,
    #     weighted = TRUE
    # )

    if (num_starts > 1) {
        set.seed(seed)
    }

    for (i in seq_len(num_starts)) {
        cluster_result <- leidenbase::leiden_find_partition(
            igraph = g,
            edge_weights = igraph::E(g)$weight,
            resolution_parameter = resolution,
            num_iter = num_iter,
            seed = ifelse(num_starts == 1, seed, NULL),
            initial_membership = initial_membership,
            ...
        )

        if (is.null(best_quality) || cluster_result$quality > best_quality) {
            best_quality <- cluster_result$quality
            best_mb <- cluster_result$membership
        }
    }

    return(best_mb)
}

clustering_functions <- function(object,
                                 resolution,
                                 seed,
                                 algorithm = 4,
                                 num_iters = 10,
                                 num_starts = 10,
                                 ...) {
    # FIXME change it to 1:3 once you will be able to share an igraph
    if (algorithm %in% 1:4) {
        return(seurat_clustering(
            object = object,
            resolution = resolution,
            seed = seed,
            algorithm = algorithm,
            num_start = num_starts,
            num_iter = num_iters,
            ...
        ))
    }

    return(
        leiden_clustering(
            g = object,
            resolution = resolution,
            seed = seed,
            num_iter = num_iters,
            num_starts = num_starts,
            ...
        )
    )
}

#### Automatic ####
rank_configs <- function(ecc_list, rank_by = "top_qt_max", return_type = "order") {
    if (!(rank_by %in% names(ranking_functions))) {
        rank_by <- "top_qt_max"
    }

    ranking_fun <- ranking_functions[[rank_by]]
    scores <- sapply(ecc_list, ranking_fun)

    if (return_type == "order") {
        return(order(scores, decreasing = TRUE))
    }

    return(length(ecc_list) + 1 - rank(scores, ties.method = "first"))
}

#' Assessment of Stability for Graph Clustering
#' @description Evaluates the stability of different graph clustering methods
#' in the clustering pipeline. The method will iterate through different values of
#' the resolution parameter and compare, using the EC Consistency score, the
#' partitions obtained at different seeds.
#'
#' @param expression_matrix An expression matrix having the features on the rows
#' and the cells on the columns.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user. Defaults to `100`.
#' @param n_neigh_sequence A sequence of the number of nearest neighbours.
#' @param resolution_sequence A sequence of resolution values. The resolution parameter
#' controls the coarseness of the clustering. The higher the resolution, the more
#' clusters will be obtained. The resolution parameter is used in the community
#' detection algorithms.
#' @param features_sets A list of the feature sets. A feature set is a list of genes
#' from the expression matrix that will be used in the dimensionality reduction.
#' @param steps A list with the same names as `feature_sets`. Each name has assigned
#' a ector containing the sizes of the subsets; negative values will
#' be interpreted as using all features.
#' @param seed_sequence A custom seed sequence; if the value is NULL, the
#' sequence will be built starting from 1 with a step of 100.
#' @param graph_reduction_embedding The type of dimensionality reduction used for
#' the graph construction. The options are "PCA" and "UMAP". Defaults to `PCA`.
#' @param include_umap_nn_assessment A boolean value indicating if the UMAP embeddings
#' will be used for the nearest neighbours assessment. Defaults to `FALSE`.
#' @param n_top_configs The number of top configurations that will be used for the
#' downstream analysis in the dimensionality reduction step. Defaults to `3`.
#' @param ranking_criterion The criterion used for ranking the configurations from
#' the dimensionality reduction step. The options are "iqr", "median", "max", "top_qt",
#' "top_qt_max", "iqr_median", "iqr_median_coeff" and "mean". Defaults to `iqr`.
#' @param overall_summary A function used to summarize the stability of the configurations
#' from the dimensionality reduction step across the different resolution values.
#' The options are "median", "max", "top_qt", "top_qt_max", "iqr", "iqr_median",
#' "iqr_median_coeff" and "mean". Defaults to `median`.
#' @param ecs_threshold The ECS threshold used for merging similar clusterings.
#' @param matrix_processing A function that will be used to process the data matrix
#' by using a dimensionality reduction technique. The function should have
#' one parameter, the data matrix, and should return an embedding describing the
#' reduced space. By default, the function will use the precise PCA method with
#' `prcomp`.
#' @param umap_arguments A list containing the arguments that will be passed to the
#' UMAP function. Refer to the `uwot::umap` function for more details.
#' @param prune_value Argument indicating whether to prune the SNN graph. If the
#' value is 0, the graph won't be pruned. If the value is between 0 and 1, the
#' edges with weight under the pruning value will be removed. If the value is
#' -1, the highest pruning value will be calculated automatically and used.
#' @param algorithm_dim_reduction An index indicating the community detection
#' algorithm that will be used in the Dimensionality reduction step.
#' @param algorithm_graph_construct An index indicating the community detection
#' algorithm that will be used in the Graph construction step.
#' @param algorithms_clustering_assessment An index indicating which community
#' detection algorithm will be used for the clustering step: Louvain (1),
#' Louvain refined (2), SLM (3) or Leiden (4). More details can be found in
#' the Seurat's `FindClusters` function.
#' @param clustering_arguments A list containing the arguments that will be passed to the
#' community detection algorithm, such as the number of iterations and the number of starts.
#' Refer to the Seurat's `FindClusters` function for more details.
#' @param verbose Boolean value used for displaying the progress
#' of the assessment.
#' @param temp_file The path to the file where the object will be saved.
#' @param save_temp A boolean value indicating if the object
#' will be saved to a file.
#'
#' @return A list having two fields:
#'
#' * all - a list that contains, for each clustering method and each resolution
#' value, the EC consistency between the partitions obtained by changing the seed
#' * filtered - similar to `all`, but for each configuration, we determine the
#' number of clusters that appears the most and use only the partitions with this
#' size
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(2024)
#' # create an already-transposed artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(20 * 10), runif(30 * 10, min = 3, max = 4)),
#'     nrow = 10, byrow = FALSE
#' )
#' colnames(expr_matrix) <- as.character(seq_len(ncol(expr_matrix)))
#' rownames(expr_matrix) <- paste("feature", seq_len(nrow(expr_matrix)))
#'
#' autom_object <- automatic_stability_assessment(
#'     expression_matrix = expr_matrix,
#'     n_repetitions = 3,
#'     n_neigh_sequence = c(5),
#'     resolution_sequence = c(0.1, 0.5),
#'     features_sets = list(
#'         "set1" = rownames(expr_matrix)
#'     ),
#'     steps = list(
#'         "set1" = c(5, 7)
#'     ),
#'     umap_arguments = list(
#'         # the following parameters have been modified
#'         # from the default values to ensure that
#'         # the function will run under 5 seconds
#'         n_neighbors = 3,
#'         approx_pow = TRUE,
#'         n_epochs = 0,
#'         init = "random",
#'         min_dist = 0.3
#'     ),
#'     n_top_configs = 1,
#'     algorithms_clustering_assessment = 1,
#'     save_temp = FALSE,
#'     verbose = FALSE
#' )
#'
#' # the object can be further used to plot the assessment results
#' plot_feature_overall_stability_boxplot(autom_object$feature_stability)
#' plot_n_neigh_ecs(autom_object$set1$"5"$nn_stability)
#' plot_k_n_partitions(autom_object$set1$"5"$clustering_stability)
#' }
automatic_stability_assessment <- function(expression_matrix,
                                           n_repetitions,
                                           n_neigh_sequence,
                                           resolution_sequence,
                                           features_sets,
                                           steps,
                                           seed_sequence = NULL,
                                           graph_reduction_embedding = "PCA",
                                           include_umap_nn_assessment = FALSE,
                                           n_top_configs = 3,
                                           ranking_criterion = "iqr",
                                           overall_summary = "median",
                                           ecs_threshold = 1, # do we really need it?,
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
                                           algorithm_dim_reduction = 1,
                                           algorithm_graph_construct = 1,
                                           algorithms_clustering_assessment = 1:3,
                                           clustering_arguments = list(),
                                           verbose = TRUE,
                                           temp_file = NULL,
                                           save_temp = TRUE) {
    file_already_exists <- ifelse(is.null(temp_file), TRUE, file.exists(temp_file))
    cell_names <- colnames(expression_matrix)
    if (is.null(cell_names)) {
        cell_names <- paste0("cell_", seq_len(ncol(expression_matrix)))
        colnames(expression_matrix) <- cell_names
    }

    feature_set_names <- names(steps)
    n_configs <- length(feature_set_names)
    total_nconfigs <- 0
    # TODO: check if the name of the feature sets are the same with the steps

    # FEATURE STABILITY
    # sort the steps in an increasing order
    for (set_name in feature_set_names) {
        steps[[set_name]][steps[[set_name]] <= 0 |
            steps[[set_name]] > length(features_sets[[set_name]])] <- length(features_sets[[set_name]])
        steps[[set_name]] <- sort(unique(steps[[set_name]]))
    }

    feature_stability_object <- list(
        by_steps = list(),
        incremental = list(),
        embedding_list = list(),
        pca_list = list()
    )

    if (verbose) {
        print(glue::glue("[{Sys.time()}] Assessing the stability of the dimensionality reduction"))
    }

    for (set_name in feature_set_names) {
        if (verbose) {
            print(set_name)
        }
        temp_object <- assess_feature_stability(
            data_matrix = expression_matrix,
            feature_set = features_sets[[set_name]],
            steps = steps[[set_name]],
            feature_type = set_name,
            resolution = resolution_sequence,
            n_repetitions = n_repetitions,
            seed_sequence = seed_sequence,
            graph_reduction_type = graph_reduction_embedding,
            ecs_thresh = ecs_threshold,
            matrix_processing = matrix_processing,
            umap_arguments = umap_arguments,
            prune_value = prune_value,
            clustering_algorithm = algorithm_dim_reduction,
            clustering_arguments = clustering_arguments,
            verbose = verbose
        )

        feature_stability_object$by_steps[[set_name]] <- temp_object$by_steps[[1]]
        feature_stability_object$incremental[[set_name]] <- temp_object$incremental[[1]]
        feature_stability_object$embedding_list[[set_name]] <- temp_object$embedding_list[[1]]

        # NOTE only for backwards compatability reasons
        if (!is.null(temp_object$pca_list)) {
            feature_stability_object$pca_list[[set_name]] <- temp_object$pca_list[[1]]
        }

        if (!file_already_exists && save_temp) {
            saveRDS(feature_stability_object, temp_file)
        }
    }

    feature_configs <- list(feature_stability = feature_stability_object)
    if (!file_already_exists && save_temp) {
        saveRDS(feature_configs, temp_file)
    }

    if (verbose) {
        print(glue::glue("[{Sys.time()}] Choosing the top {n_top_configs}"))
    }

    for (i in seq_along(feature_set_names)) {
        set_name <- feature_set_names[i]
        # rank the sizes based on the consistency and incremental stability
        ecc_list <- lapply(feature_stability_object$by_steps[[set_name]], function(by_step) {
            sapply(by_step, function(by_res) {
                ranking_functions[[overall_summary]](by_res$ecc)
            })
        })
        ecc_list_medians <- order(sapply(ecc_list, stats::median), decreasing = TRUE)
        ecc_above_median <- ecc_list_medians[seq_len(ceiling((length(ecc_list_medians) + 1) / 2))] # which(ecc_list_medians >= stats::median(ecc_list_medians))

        # ecc_ranking <- rank_configs(ecc_list, rank_by = ranking_criterion, return_type = "rank")
        incremental_list <- c(
            list(sapply(
                feature_stability_object$incremental[[set_name]][[1]],
                function(by_res) {
                    ranking_functions[[overall_summary]](by_res)
                }
            )),
            lapply(
                feature_stability_object$incremental[[set_name]],
                function(by_step) {
                    sapply(by_step, function(by_res) {
                        ranking_functions[[overall_summary]](by_res)
                    })
                }
            )
        )

        incremental_medians <- order(sapply(incremental_list, stats::median), decreasing = TRUE)
        incremental_above_median <- incremental_medians[seq_len(ceiling((length(incremental_medians) + 1) / 2))] # which(incremental_medians >= stats::median(incremental_medians))
        eligible_steps <- intersect(ecc_above_median, incremental_above_median)
        n_top_configs <- min(n_top_configs, length(eligible_steps))
        total_nconfigs <- total_nconfigs + n_top_configs
        feature_stability_ranking <- rank_configs(lapply(eligible_steps, function(i) {
            incremental_list[[i]] + ecc_list[[i]]
        }), rank_by = ranking_criterion, return_type = "order")
        # feature_stability_ranking <- order(

        feature_configs[[set_name]] <- list()

        # generate the pca and umap embeddings for each config
        actual_npcs <- formals(matrix_processing)$actual_npcs
        for (j in seq_len(n_top_configs)) {
            n_steps <- steps[[i]][eligible_steps][feature_stability_ranking[j]]
            current_features <- features_sets[[set_name]][seq_len(n_steps)]
            if (is.null(feature_stability_object$pca_list[[set_name]][[as.character(n_steps)]])) {
                pca_emb <- matrix_processing(expression_matrix[current_features, ])
            } else {
                pca_emb <- feature_stability_object$pca_list[[set_name]][[as.character(n_steps)]]
            }

            umap_emb <- feature_stability_object$embedding_list[[set_name]][[as.character(n_steps)]]

            feature_configs[[set_name]][[as.character(n_steps)]] <-
                list(
                    pca = pca_emb,
                    umap = umap_emb,
                    stable_config = list(
                        feature_set = feature_set_names[i],
                        n_features = n_steps,
                        n_pcs = min(actual_npcs, n_steps %/% 2)
                    )
                )
        }
        feature_configs[[set_name]]$feature_list <- features_sets[[set_name]][seq_len(max(steps[[set_name]]))]
    }

    feature_configs$feature_stability$pca_list <- NULL
    rm(feature_stability_object)
    gc()

    # GRAPH CONSTRUCTION: CONNECTED COMPS
    if (verbose) {
        print(glue::glue("[{Sys.time()}] Assessing the stability of the connected components"))
        pb <- progress::progress_bar$new(
            format = ":featurename - :featuresize [:bar] eta: :eta  total elapsed: :elapsed",
            total = total_nconfigs,
            show_after = 0,
            clear = FALSE,
            width = 80
        )
    }

    for (set_name in feature_set_names) {
        n_names <- length(feature_configs[[set_name]])
        for (n_steps in names(feature_configs[[set_name]])[seq_len(n_names - 1)]) {
            if (verbose) {
                pb$tick(0, tokens = list(featurename = set_name, featuresize = n_steps))
            }

            feature_configs[[set_name]][[n_steps]][["nn_conn_comps"]] <- get_nn_conn_comps(
                embedding = feature_configs[[set_name]][[n_steps]]$pca,
                n_neigh_sequence = n_neigh_sequence,
                n_repetitions = n_repetitions,
                seed_sequence = seed_sequence,
                include_umap = include_umap_nn_assessment,
                umap_arguments = umap_arguments
            )
            if (verbose) {
                pb$tick(tokens = list(featurename = set_name, featuresize = n_steps))
            }

            if (!file_already_exists) {
                saveRDS(feature_configs, temp_file)
            }
        }
    }

    # GRAPH CONSTRUCTION: NN STABILITY
    if (verbose) {
        print(glue::glue("[{Sys.time()}] Assessing the stability of the graph construction parameters"))
        pb <- progress::progress_bar$new(
            format = ":featurename - :featuresize [:bar] eta: :eta  total elapsed: :elapsed",
            total = total_nconfigs,
            show_after = 0,
            clear = FALSE,
            width = 80
        )
    }

    for (set_name in feature_set_names) {
        n_names <- length(feature_configs[[set_name]])
        for (n_steps in names(feature_configs[[set_name]])[seq_len(n_names - 1)]) {
            if (verbose) {
                pb$tick(0, tokens = list(featurename = set_name, featuresize = n_steps))
            }

            feature_configs[[set_name]][[n_steps]][["nn_stability"]] <- assess_nn_stability(
                embedding = feature_configs[[set_name]][[n_steps]]$pca,
                n_neigh_sequence = n_neigh_sequence,
                n_repetitions = n_repetitions,
                seed_sequence = seed_sequence,
                graph_reduction_type = "PCA",
                ecs_thresh = ecs_threshold,
                prune_value = prune_value,
                clustering_algorithm = algorithm_graph_construct,
                clustering_arguments = clustering_arguments
            )

            if (include_umap_nn_assessment) {
                feature_configs[[set_name]][[n_steps]][["nn_stability"]] <- mapply(c,
                    feature_configs[[set_name]][[n_steps]][["nn_stability"]],
                    assess_nn_stability(
                        embedding = feature_configs[[set_name]][[n_steps]]$pca,
                        n_neigh_sequence = n_neigh_sequence,
                        n_repetitions = n_repetitions,
                        seed_sequence = seed_sequence,
                        graph_reduction_type = "UMAP",
                        ecs_thresh = ecs_threshold,
                        prune_value = prune_value,
                        clustering_algorithm = algorithm_graph_construct,
                        umap_arguments = umap_arguments
                    ),
                    SIMPLIFY = FALSE
                )
            }

            if (verbose) {
                pb$tick(tokens = list(featurename = set_name, featuresize = n_steps))
            }

            if (!file_already_exists) {
                saveRDS(feature_configs, temp_file)
            }
        }
    }

    # choose best graph construction parameters for each config
    for (set_name in feature_set_names) {
        n_names <- length(feature_configs[[set_name]])
        for (n_steps in names(feature_configs[[set_name]])[seq_len(n_names - 1)]) {
            nn_ecc_list <- list()
            nn_config_list <- c()
            nn_list <- c()
            index <- 1
            for (config_name in names(feature_configs[[set_name]][[n_steps]]$nn_stability$n_neigh_ec_consistency)) {
                for (n_neigh in names(feature_configs[[set_name]][[n_steps]]$nn_stability$n_neigh_ec_consistency[[config_name]])) {
                    current_ecc <- feature_configs[[set_name]][[n_steps]]$nn_stability$n_neigh_ec_consistency[[config_name]][[n_neigh]]

                    nn_ecc_list[[index]] <- current_ecc
                    index <- index + 1
                    nn_config_list <- c(nn_config_list, config_name)
                    nn_list <- c(nn_list, n_neigh)
                }
            }


            ecc_medians <- sapply(nn_ecc_list, stats::median)
            eligible_configs <- which(ecc_medians >= stats::fivenum(ecc_medians)[4])
            best_config_index <- rank_configs(nn_ecc_list[eligible_configs], rank_by = ranking_criterion, return_type = "order")[1]
            best_config <- nn_config_list[eligible_configs][best_config_index]
            best_nn <- as.integer(nn_list[eligible_configs][best_config_index])

            split_configs <- strsplit(best_config, "_")[[1]]
            base_embedding <- tolower(split_configs[length(split_configs) - 1])
            graph_type <- split_configs[length(split_configs)] # NOTE update if you decide to remove ecs thresh
            # feature_configs[[set_name]][[n_steps]][["adj_matrix"]] <- Seurat::FindNeighbors(
            #     object = feature_configs[[set_name]][[n_steps]][[base_embedding]],
            #     k.param = best_nn,
            #     verbose = FALSE,
            #     nn.method = "rann",
            #     compute.SNN = FALSE
            # )$nn
            feature_configs[[set_name]][[n_steps]][["adj_matrix"]] <- getNNmatrix(
                RANN::nn2(
                    feature_configs[[set_name]][[n_steps]][[base_embedding]],
                    k = best_nn,
                )$nn.idx,
                best_nn,
                0,
                -1
            )$nn
            highest_prune_param <- 0

            if (graph_type == "snn") {
                if (prune_value >= 0) {
                    feature_configs[[set_name]][[n_steps]][["adj_matrix"]] <- computeSNN(
                        feature_configs[[set_name]][[n_steps]][["adj_matrix"]],
                        best_nn,
                        prune_value
                    )
                    highest_prune_param <- prune_value
                } else {
                    highest_prune_param <- get_highest_prune_param(
                        feature_configs[[set_name]][[n_steps]][["adj_matrix"]],
                        best_nn
                    )
                    feature_configs[[set_name]][[n_steps]][["adj_matrix"]] <- highest_prune_param$adj_matrix
                    highest_prune_param <- highest_prune_param$prune_value
                    gc()
                }
            }
            rownames(feature_configs[[set_name]][[n_steps]][["adj_matrix"]]) <- cell_names
            colnames(feature_configs[[set_name]][[n_steps]][["adj_matrix"]]) <- cell_names

            feature_configs[[set_name]][[n_steps]]$stable_config[["base_embedding"]] <- base_embedding
            feature_configs[[set_name]][[n_steps]]$stable_config[["graph_type"]] <- graph_type
            feature_configs[[set_name]][[n_steps]]$stable_config[["n_neighbours"]] <- as.numeric(best_nn)
            feature_configs[[set_name]][[n_steps]]$stable_config[["prune_param"]] <- highest_prune_param

            if (!file_already_exists && save_temp) {
                saveRDS(feature_configs, temp_file)
            }
        }
    }


    # GRAPH CLUSTERING: CLUSTERING METHOD
    if (verbose) {
        print(glue::glue("[{Sys.time()}] Assessing the stability of the graph clustering method"))
    }
    for (set_name in feature_set_names) {
        n_names <- length(feature_configs[[set_name]])
        for (n_steps in names(feature_configs[[set_name]])[seq_len(n_names - 1)]) {
            if (verbose) {
                print(paste(set_name, n_steps))
            }
            feature_configs[[set_name]][[n_steps]][["clustering_stability"]] <- assess_clustering_stability(
                graph_adjacency_matrix = feature_configs[[set_name]][[n_steps]]$adj_matrix,
                resolution = resolution_sequence,
                n_repetitions = n_repetitions,
                seed_sequence = seed_sequence,
                ecs_thresh = ecs_threshold,
                clustering_algorithm = algorithms_clustering_assessment,
                clustering_arguments = clustering_arguments,
                verbose = verbose
            )

            if (!file_already_exists && save_temp) {
                saveRDS(feature_configs, temp_file)
            }
        }
    }


    if (!file_already_exists && save_temp) {
        saveRDS(feature_configs, temp_file)
    }

    # NOTE add this class to uniquely identify the object
    class(feature_configs) <- c(class(feature_configs), "ClustAssess")
    return(feature_configs)
}

