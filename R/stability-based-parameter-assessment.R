# TODO check user interrupt
# TODO add logging options to files - maybe verbose lvels
#' @importFrom foreach %dopar%
NULL

ranking_functions <- list(
  "median" = median,
  "max" = max,
  "top_qt" = function(arr) {
    summary(arr)[5]
  },
  "top_qt_max" = function(arr) {
    summary(arr)[5] + summary(arr)[6]
  },
  "mean" = mean
)

# wrapper of the Seurat's `FindClusters` method, that returns
# only the membership vector
seurat_clustering <- function(object, resolution, seed, algorithm = 4, ...) {
  cluster_result <- Seurat::FindClusters(
    object,
    resolution = resolution,
    random.seed = seed,
    algorithm = algorithm,
    verbose = FALSE,
    ...
  )
  cluster_result[[colnames(cluster_result)[1]]]
}

leiden_clustering <- function(object, resolution, seed, ...) {
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
#' @param graph_adjacency_matrix A square adjacency matrix based on which an igraph
#' object will be built.
#' @param resolution A sequence of resolution values.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence A custom seed sequence; if the value is NULL, the
#' sequence will be built starting from 1 with a step of 100.
#' @param ecs_thresh The ECS threshold used for merging similar clusterings.
#' @param ncores The number of parallel R instances that will run the code.
#' If the value is set to 1, the code will be run sequentially.
#' @param verbose Boolean value used for displaying the progress bar.
#' @param algorithm An index or a list of indexes indicating which community detection
#' algorithm will be used: Louvain (1), Louvain refined (2), SLM (3) or Leiden (4).
#' More details can be found in the Seurat's `FindClusters` function.
#'
#' @return A list having two fields:
#'
#' * all - a list that contains, for each clustering method and each resolution
#' value, the EC consistency between the partitions obtained by changing the seed
#' * filtered - similar to `all`, but for each configuration, we determine the
#' number of clusters that appears the most and use only the partitions with this
#' size
#'
#' @md
#' @export
automatic_stability_assessment <- function(expression_matrix, # expr matrix
                                           n_cores,
                                           n_repetitions,
                                           n_neigh_sequence,
                                           resolution_sequence,
                                           seed_sequence = NULL,
                                           features_sets, # list with the names of the feature sets (HV, MA etc) and the list of top x genes
                                           steps, # list with the names of the feature sets (HV, mA) and the list of steps
                                           graph_reduction_embedding = "PCA",
                                           n_top_configs = 3,
                                           ranking_criterion = "median",
                                           npcs = 30,
                                           verbose = TRUE,
                                           ecs_threshold = 1, # do we really need it?,
                                           ...) {
  # store the additional arguments used by umap in a list
  suppl_args <- list(...)

  feature_set_names <- names(steps)
  n_configs <- length(feature_set_names)
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
    embedding_list = list()
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
      resolution = resolution_sequence,
      steps = steps[[set_name]],
      feature_type = set_name,
      graph_reduction_type = graph_reduction_embedding,
      n_repetitions = n_repetitions,
      seed_sequence = seed_sequence,
      npcs = npcs,
      ecs_thresh = ecs_threshold,
      ncores = n_cores,
      algorithm = 1,
      verbose = verbose,
      ...
      # min_dist = 0.3, n_neighbors = 30, metric = "cosine"
    )

    feature_stability_object$by_steps[[set_name]] <- temp_object$by_steps[[1]]
    feature_stability_object$incremental[[set_name]] <- temp_object$incremental[[1]]
    feature_stability_object$embedding_list[[set_name]] <- temp_object$embedding_list[[1]]
  }

  feature_configs <- list(feature_stability = feature_stability_object)

  for (i in seq_along(feature_set_names)) {
    set_name <- feature_set_names[i]
    # rank the sizes based on the consistency and incremental stability
    ecc_list <- lapply(feature_stability_object$by_steps[[set_name]], function(by_step) {
      sapply(by_step, function(by_res) {
        ranking_functions[[ranking_criterion]](by_res$ecc)
      })
    })

    ecc_ranking <- rank_configs(ecc_list, rank_by = ranking_criterion, return_type = "rank")
    incremental_list <- c(
      list(sapply(
        feature_stability_object$incremental[[set_name]][[1]],
        function(by_res) {
          ranking_functions[[ranking_criterion]](by_res)
        }
      )),
      lapply(
        feature_stability_object$incremental[[set_name]],
        function(by_step) {
          sapply(by_step, function(by_res) {
            ranking_functions[[ranking_criterion]](by_res)
          })
        }
      )
    )
    incremental_ranking <- rank_configs(incremental_list, rank_by = ranking_criterion, return_type = "rank")
    feature_stability_ranking <- order(incremental_ranking + ecc_ranking)

    current_n_top <- min(n_top_configs, length(incremental_ranking))

    feature_configs[[set_name]] <- list()

    # generate the pca and umap embeddings for each config
    for (j in 1:current_n_top) {
      n_steps <- steps[[i]][feature_stability_ranking[j]]
      current_features <- features_sets[[set_name]][1:n_steps]
      suppressWarnings(pca_emb <- Seurat::RunPCA(expression_matrix[current_features, ],
        npcs = min(npcs, n_steps %/% 2),
        approx = FALSE,
        verbose = FALSE
      )@cell.embeddings)

      suppressWarnings(
        umap_emb <- do.call(
          Seurat::RunUMAP,
          c(list(object = pca_emb, verbose = FALSE), suppl_args)
        )@cell.embeddings
      )

      feature_configs[[set_name]][[as.character(n_steps)]] <-
        list(
          pca = pca_emb,
          umap = umap_emb,
          stable_config = list(
            feature_set = feature_set_names[i],
            n_features = n_steps,
            n_pcs = min(npcs, n_steps %/% 2)
          )
        )
    }
    feature_configs[[set_name]]$feature_list <- features_sets[[set_name]][seq_len(max(steps[[set_name]]))]
  }
 
  # GRAPH CONSTRUCTION: CONNECTED COMPS
  if (verbose) {
    print(glue::glue("[{Sys.time()}] Assessing the stability of the connected components"))
    pb <- progress::progress_bar$new(
      format = ":featurename - :featuresize [:bar] eta: :eta  total elapsed: :elapsed",
      total = n_configs * n_top_configs,
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
        ncores = n_cores,
        seed_sequence = seed_sequence,
        ...
      )
      if (verbose) {
        pb$tick(tokens = list(featurename = set_name, featuresize = n_steps))
      }
    }
  }

  # GRAPH CONSTRUCTION: NN STABILITY
  if (verbose) {
    print(glue::glue("[{Sys.time()}] Assessing the stability of the graph construction parameters"))
    pb <- progress::progress_bar$new(
      format = ":featurename - :featuresize [:bar] eta: :eta  total elapsed: :elapsed",
      total = n_configs * n_top_configs,
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

      feature_configs[[set_name]][[n_steps]][["nn_stability"]] <- mapply(c,
        assess_nn_stability(
          embedding = feature_configs[[set_name]][[n_steps]]$pca,
          n_neigh_sequence = n_neigh_sequence,
          n_repetitions = n_repetitions,
          graph_reduction_type = "PCA",
          ecs_thresh = 1,
          ncores = n_cores,
          algorithm = 1,
        ),
        assess_nn_stability(
          embedding = feature_configs[[set_name]][[n_steps]]$pca,
          n_neigh_sequence = n_neigh_sequence,
          n_repetitions = n_repetitions,
          graph_reduction_type = "UMAP",
          ecs_thresh = 1,
          ncores = n_cores,
          algorithm = 1,
          ...
        ),
        SIMPLIFY = FALSE
      )

      if (verbose) {
        pb$tick(tokens = list(featurename = set_name, featuresize = n_steps))
      }
    }
  }

  # choose best graph construction parameters for each config
  for (set_name in feature_set_names) {
    n_names <- length(feature_configs[[set_name]])
    for (n_steps in names(feature_configs[[set_name]])[seq_len(n_names - 1)]) {
      best_ecc <- 0
      best_config <- NULL
      best_nn <- NULL
      for (config_name in names(feature_configs[[set_name]][[n_steps]]$nn_stability$n_neigh_ec_consistency)) {
        for (n_neigh in names(feature_configs[[set_name]][[n_steps]]$nn_stability$n_neigh_ec_consistency[[config_name]])) {
          current_ecc <- ranking_functions[[ranking_criterion]](feature_configs[[set_name]][[n_steps]]$nn_stability$n_neigh_ec_consistency[[config_name]][[n_neigh]])
          if (current_ecc > best_ecc) {
            best_ecc <- current_ecc
            best_config <- config_name
            best_nn <- as.numeric(n_neigh)
          }
        }
      }

      split_configs <- strsplit(best_config, "_")[[1]]
      base_embedding <- tolower(split_configs[length(split_configs) - 2])
      graph_type <- split_configs[length(split_configs) - 1] # update if you decide to remove ecs thresh
      # TODO add automatic prune param
      highest_prune_param <- get_highest_prune_param(
        feature_configs[[set_name]][[n_steps]][[base_embedding]],
        best_nn
      )
      feature_configs[[set_name]][[n_steps]][["adj_matrix"]] <- Seurat::FindNeighbors(
        object = feature_configs[[set_name]][[n_steps]][[base_embedding]],
        k.param = best_nn,
        verbose = FALSE,
        nn.method = "rann",
        prune.SNN = highest_prune_param
      )[[graph_type]]

      feature_configs[[set_name]][[n_steps]]$stable_config[["base_embedding"]] <- base_embedding
      feature_configs[[set_name]][[n_steps]]$stable_config[["graph_type"]] <- graph_type
      feature_configs[[set_name]][[n_steps]]$stable_config[["n_neighbours"]] <- as.numeric(best_nn)
      feature_configs[[set_name]][[n_steps]]$stable_config[["prune_param"]] <- highest_prune_param
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
        ecs_thresh = ecs_threshold,
        ncores = n_cores,
        algorithm = 1:3,
        verbose = verbose
      )
    }
  }

  return(feature_configs)
}
