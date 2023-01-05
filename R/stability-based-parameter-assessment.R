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



#### PCA ####



######################## Graph construction ######################################



############################## Clustering ########################################



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
                                  ecs_threshold = 1, # do we really need it?,
                                  ...) {
  # store the additional arguments used by umap in a list
  suppl_args <- list(...)

  feature_set_names <- names(steps)
  config_names <- paste(names(steps), graph_reduction_embedding, ecs_threshold, sep = "_")
  n_configs <- length(config_names)
  # TODO: check if the name of the feature sets are the same with the steps


  # FEATURE STABILITY
  # sort the steps in an increasing order
  for (set_name in feature_set_names) {
    steps[[set_name]][steps[[set_name]] <= 0 |
      steps[[set_name]] > length(features_sets[[set_name]])] <- length(features_sets[[set_name]])
    steps[[set_name]] <- sort(unique(steps[[set_name]]))
  }

  feature_stability_object <- list(
    steps_stability = list(),
    incremental_stability = list()
  )

  for (set_name in feature_set_names) {
    print(set_name)
    temp_object <- assess_feature_stability(
      data_matrix = expression_matrix,
      feature_set = features_sets[[set_name]],
      steps = steps[[set_name]],
      feature_type = set_name,
      graph_reduction_type = graph_reduction_embedding,
      n_repetitions = n_repetitions,
      seed_sequence = seed_sequence,
      npcs = npcs,
      ecs_thresh = ecs_threshold,
      ncores = n_cores,
      algorithm = 1,
      ...
      # min_dist = 0.3, n_neighbors = 30, metric = "cosine"
    )

    config_name <- names(temp_object$steps_stability)[1]
    feature_stability_object$steps_stability[[config_name]] <- temp_object$steps_stability[[1]]
    feature_stability_object$incremental_stability[[config_name]] <- temp_object$incremental_stability[[1]]
  }

  names(steps) <- config_names
  names(features_sets) <- config_names

  feature_configs <- list(feature_importance = feature_stability_object)

  for (i in 1:length(config_names)) {
    set_name <- config_names[i]
    ecc_list <- lapply(feature_stability_object$steps_stability[[set_name]], function(x) {
      x$ecc
    })
    ecc_ranking <- rank_configs(ecc_list, rank_by = ranking_criterion, return_type = "rank")
    incremetal_list <- c(
      feature_stability_object$incremental_stability[[set_name]][1],
      feature_stability_object$incremental_stability[[set_name]]
    )
    incremental_ranking <- rank_configs(incremetal_list, rank_by = ranking_criterion, return_type = "rank")
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

      suppressWarnings(umap_emb <- do.call(RunUMAP, c(list(object = pca_emb, verbose = FALSE), suppl_args))@cell.embeddings)

      feature_configs[[set_name]][[as.character(n_steps)]] <-
        list(
          pca = pca_emb,
          umap = umap_emb,
          stable_config = list(
            feature_set = feature_set_names[i],
            n_features = n_steps
          )
        )
    }
  }

  # GRAPH CONSTRUCTION: CONNECTED COMPS
  print("Assessing the stability of the connected components")

  for (config_name in config_names) {
    print(config_name)
    for (n_steps in names(feature_configs[[config_name]])) {
      feature_configs[[config_name]][[n_steps]][["nn_conn_comps"]] <- get_nn_conn_comps(
          embedding = feature_configs[[config_name]][[n_steps]]$pca,
          n_neigh_sequence = n_neigh_sequence,
          n_repetitions = n_repetitions,
          ncores = n_cores,
          seed_sequence = seed_sequence,
          # config_name = paste(config_name, n_steps, sep = "_"),
          ...
        )
    }
  }

  # GRAPH CONSTRUCTION: NN IMPORTANCE
  print("Assessing the stability of the graph construction parameters")

  for (config_name in config_names) {
    print(config_name)
    for (n_steps in names(feature_configs[[config_name]])) {
      print(n_steps)
      feature_configs[[config_name]][[n_steps]][["nn_importance"]] <- mapply(c,
        assess_nn_stability(
          embedding = feature_configs[[config_name]][[n_steps]]$pca,
          n_neigh_sequence = n_neigh_sequence,
          n_repetitions = n_repetitions,
          graph_reduction_type = "PCA",
          ecs_thresh = 1,
          ncores = n_cores,
          algorithm = 1,
        ),
        assess_nn_stability(
          embedding = feature_configs[[config_name]][[n_steps]]$pca,
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
    }
  }


  # choose best graph construction parameters for each config
  for (feature_config_name in config_names) {
    for (n_steps in names(feature_configs[[feature_config_name]])) {
      best_ecc <- 0
      best_config <- NULL
      best_nn <- NULL
      for (config_name in names(feature_configs[[feature_config_name]][[n_steps]]$nn_importance$n_neigh_ec_consistency)) {
        for (n_neigh in names(feature_configs[[feature_config_name]][[n_steps]]$nn_importance$n_neigh_ec_consistency[[config_name]])) {
          current_ecc <- ranking_functions[[ranking_criterion]](feature_configs[[feature_config_name]][[n_steps]]$nn_importance$n_neigh_ec_consistency[[config_name]][[n_neigh]])
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
      # TODO add atuomatic prune param
      feature_configs[[feature_config_name]][[n_steps]][["adj_matrix"]] <- FindNeighbors(
        object = feature_configs[[feature_config_name]][[n_steps]][[base_embedding]],
        k.param = best_nn,
        verbose = F,
        nn.method = "rann"
      )[[graph_type]]

      feature_configs[[feature_config_name]][[n_steps]]$stable_config[["base_embedding"]] <- base_embedding
      feature_configs[[feature_config_name]][[n_steps]]$stable_config[["graph_type"]] <- graph_type
      feature_configs[[feature_config_name]][[n_steps]]$stable_config[["n_neighbours"]] <- as.numeric(best_nn)
    }
  }


  # GRAPH CLUSTERING: CLUSTERING METHOD
  print("Assessing the stability of the graph clustering method")

  for (config_name in config_names) {
    print(config_name)
    for (n_steps in names(feature_configs[[config_name]])) {
      print(n_steps)
      feature_configs[[config_name]][[n_steps]][["clustering_importance"]] <- assess_clustering_stability(
        graph_adjacency_matrix = feature_configs[[config_name]][[n_steps]]$adj_matrix,
        resolution = resolution_sequence,
        n_repetitions = n_repetitions,
        ecs_thresh = ecs_threshold,
        ncores = n_cores,
        algorithm = 1:4
      )

      # ecc_clustering_methods <- lapply(
      #   feature_configs[[config_name]][[n_steps]]$clustering_importance$all,
      #   function(clust_method) {
      #     c(sapply(clust_method, c))
      #   }
      # )

      # best_clustering_method <- rank_configs(ecc_clustering_methods, rank_by = ranking_criterion)[1]
      # feature_configs[[config_name]][[n_steps]]$stable_config[["clustering_method"]] <- algorithm_names[best_clustering_method]
    }
  }

  return(feature_configs)
  # GRAPH CLUSTERING: RESOLUTION GRIDSEARCH
  print("Assessing the stability of the number of clusters")
  for (feature_config_name in config_names) {
    print(feature_config_name)
    for (n_steps in names(feature_configs[[feature_config_name]])) {
      print(n_steps)
      base_embedding <- feature_configs[[feature_config_name]][[n_steps]]$stable_config[["base_embedding"]]
      graph_type <- as.numeric(feature_configs[[feature_config_name]][[n_steps]]$stable_config[["graph_type"]] == "snn")
      n_neighbours <- as.numeric(feature_configs[[feature_config_name]][[n_steps]]$stable_config[["n_neighbours"]])
      clustering_alg <- which(algorithm_names == feature_configs[[feature_config_name]][[n_steps]]$stable_config[["clustering_method"]])

      feature_configs[[feature_config_name]][[n_steps]][["resolution_importance"]] <- get_resolution_importance(
        embedding = feature_configs[[feature_config_name]][[n_steps]][[base_embedding]],
        resolution = resolution_sequence,
        n_neigh = n_neighbours,
        n_repetitions = n_repetitions,
        clustering_method = clustering_alg,
        graph_type = graph_type,
        ecs_thresh = ecs_threshold,
        ncores = n_cores
      )
    }
  }

  for (feature_config_name in config_names) {
    for (n_steps in names(feature_configs[[feature_config_name]])) {
      temp_list <- list()

      for (k in names(feature_configs[[feature_config_name]][[n_steps]]$resolution_importance$split_by_k[[1]])) {
        temp_list[[k]] <- feature_configs[[feature_config_name]][[n_steps]]$resolution_importance$split_by_k[[1]][[k]]$partitions[[1]]$mb
      }

      feature_configs[[feature_config_name]][[n_steps]][["stable_partitions"]] <- temp_list
    }
  }

  return(feature_configs)
}
