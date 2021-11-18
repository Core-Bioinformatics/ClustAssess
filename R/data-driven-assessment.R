#' @importFrom foreach %dopar%
NULL

# set these variables as global to solve the foreach warnings
utils::globalVariables(c("i", "seed", "obj"))

# wrapper of the Seurat's `FindClusters` method, that returns
# only the membership vector
seurat_clustering = function(object, resolution, seed, algorithm = 4, ...) {
  cluster.result = Seurat::FindClusters(
    object,
    resolution = resolution,
    random.seed = seed,
    algorithm = algorithm,
    ...
  )
  cluster.result[[colnames(cluster.result)[1]]]
}


# order_list = function(list_obj, order_vec) {
#   new_list = list()
#   index = 1
#
#   for (i in order_vec) {
#     new_list[[index]] = list_obj[[i]]
#     index = index + 1
#   }
#
#   new_list
# }


# generate values from an interval that simulates the loghartimic scale
generate_breaks = function(min.range, max.range) {
  start.point = min.range + 5 - min.range %% 5
  breaks.list = list()
  index = 1
  step = 5
  while (start.point < max.range) {
    breaks.list[[index]] = start.point
    index = index + 1

    start.point = start.point + step

    # the step is doubled after every two iterations
    if (index %% 2) {
      step = step * 2
    }
  }

  c(min.range, unlist(breaks.list), max.range)
}


#### PCA ####

#' Feature Stability Object Generator
#'
#' @description Creates an object that summaries the consistency of the partitions obtained
#' across different runs for varying subsets of a given feature set.
#'
#' @param data_matrix a data matrix having the features on the columns and the observations on the rows.
#' @param feature_set a set of feature names that can be found in the data matrix.
#' @param steps vector containing the sizes of the subsets; negative values will be interpreted as using all features.
#' @param n_repetitions the number of repetitions of applying the pipeline with different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence a custom seed sequence; if the value is NULL, the sequence will be built starting from 1 with a step of 100.
#' @param feature_type a name associated to the feature_set.
#' @param graph_reduction_type the graph reduction type, denoting if the graph should be built on either the PCA or the UMAP embedding.
#' @param npcs the number of principal components.
#' @param ecs_thresh the ECS threshold used for merging similar clusterings.
#' @param ncores the number of parallel R instances that will run the code. If the value is set to 1, the code will be run sequentially.
#' @param ... additional arguments passed to the umap method.
#'
#'
#' @return A list having one field associated with a step value. Each step contains a list with three fields:
#'
#' * ecs - the EC-Consistency of the partitions obtained on all repetitions
#' * embedding - one UMAP embedding generated on the feature subset
#' * most_frequent_partition - the most common partition obtained across repetitions
#' @md
#'
#' @note The algorithm assumes that the feature_set is already when performing the subsetting. For example,
#' if the user wants to analyze highly variable feature set, they should provide them sorted by their variability.
#'
#'
#' @export
#'
#' @examples
#' get_feature_stability_object(data_matrix = as.matrix(mtcars),
#'    feature_set = colnames(mtcars),
#'    steps = -1,
#'    npcs = 5,
#'    n_repetitions = 1)
get_feature_stability_object = function(data_matrix,
                                                 feature_set,
                                                 steps,
                                                 n_repetitions = 30,
                                                 seed_sequence = NULL,
                                                 feature_type = "",
                                                 graph_reduction_type = "PCA",
                                                 npcs = 30,
                                                 ecs_thresh = 1,
                                                 ncores = 1,
                                                 ...) {
  partitions_list = list()

  object_name = paste(feature_type, graph_reduction_type, ecs_thresh, sep = "_")
  return_list = list()
  return_list[[object_name]] = list()

  # create a seed sequence if it's not provided
  if (is.null(seed_sequence)) {
    seed_sequence = seq(from = 1,
                        by = 100,
                        length.out = n_repetitions)
  } else {
    n_repetitions = length(seed_sequence)
  }

  ncores = min(ncores, length(seed_sequence))

  # convert the negative steps to the size of the feature
  steps[steps <= 0 |
          steps > length(feature_set)] = length(feature_set)

  # store the additional arguments used by umap in a list
  suppl_args = list(...)
  i = 1
  while(i <= length(suppl_args)) {
    assign(names(suppl_args)[i], suppl_args[[i]])
    i = i + 1
  }
  umap_arg_names = names(suppl_args)

  for (step in steps) {
    used_features = feature_set[1:step]

    partitions_list[[as.character(step)]] = list()
    return_list[[object_name]][[as.character(step)]] = list()

    # keep only the features that we are using
    trimmed_matrix = data_matrix[, used_features]

    # the variables needed in each PSOCK process
    needed_vars = c("trimmed_matrix", "graph_reduction_type", "npcs", "umap_arg_names")

    if(ncores > 1) {
      # create a parallel backend
      my_cluster <- parallel::makeCluster(
        ncores,
        type = "PSOCK"
      )

      doParallel::registerDoParallel(cl = my_cluster)
    } else {
      # create a sequential backend
      foreach::registerDoSEQ()
    }

    all_vars = ls()

    # send the name of the umap arguments
    partitions_list[[as.character(step)]] = foreach::foreach(seed = seed_sequence,
                                                             .noexport = all_vars[!(all_vars %in% needed_vars)],
                                                             .export = names(suppl_args)) %dopar% { #
                                                               umap_args = list()

                                                               for(umap_arg_name in umap_arg_names) {
                                                                 umap_args[[umap_arg_name]] = get(umap_arg_name)
                                                               }

                                                               # calculate the PCA embedding
                                                               set.seed(seed)
                                                               pca_embedding = irlba::irlba(A = trimmed_matrix, nv = npcs)
                                                               pca_embedding = pca_embedding$u %*% diag(pca_embedding$d)
                                                               colnames(pca_embedding) = paste0("PC_", 1:ncol(pca_embedding))
                                                               rownames(pca_embedding) = rownames(trimmed_matrix)


                                                               if (graph_reduction_type == "UMAP") {
                                                                 # calculate the umap embedding
                                                                 set.seed(seed)
                                                                 umap_embedding = do.call(uwot::umap, c(list(X = pca_embedding), umap_args))
                                                                 colnames(umap_embedding) = paste0("UMAP_", 1:ncol(umap_embedding))
                                                                 rownames(umap_embedding) = rownames(trimmed_matrix)
                                                               }

                                                               # build the nn and snn graphs
                                                               neigh_matrix = Seurat::FindNeighbors(
                                                                 if (graph_reduction_type == "PCA")
                                                                   pca_embedding
                                                                 else
                                                                   umap_embedding,
                                                                 nn.method = "rann",
                                                                 verbose = FALSE
                                                               )$snn

                                                               # apply Leiden Clustering to the generated graph
                                                               cluster_results = Seurat::FindClusters(
                                                                 neigh_matrix,
                                                                 random.seed = seed,
                                                                 algorithm = 4,
                                                                 verbose = FALSE
                                                               )
                                                               cluster_results = cluster_results[, names(cluster_results)[1]]

                                                               # return the clustering
                                                               list(mb = cluster_results,
                                                                    freq = 1,
                                                                    seed = seed)
                                                             }

    # if a parallel backend was created, terminate the processes
    if(ncores > 1)
      parallel::stopCluster(cl = my_cluster)

    # generate an UMAP embedding that will be stored in the returned list
    set.seed(seed_sequence[1])
    pca_embedding = irlba::irlba(A = trimmed_matrix, nv = npcs)
    pca_embedding = pca_embedding$u %*% diag(pca_embedding$d)
    colnames(pca_embedding) = paste0("PC_", 1:ncol(pca_embedding))
    rownames(pca_embedding) = rownames(trimmed_matrix)

    set.seed(seed_sequence[1])
    umap_embedding = uwot::umap(X = pca_embedding,  ...)
    colnames(umap_embedding) = paste0("UMAP_", 1:ncol(umap_embedding))
    rownames(umap_embedding) = rownames(data_matrix)

    # merge the identical partitions into the same object
    partitions_list[[as.character(step)]] = merge_partitions(partitions_list[[as.character(step)]],
                                                             ecs_thresh = ecs_thresh,
                                                             ncores = ncores,
                                                             order = TRUE)

    # compute the EC-consistency of the partition list
    ec_consistency = weighted_element_consistency(lapply(partitions_list[[as.character(step)]], function(x) {
      x$mb
    }),
    sapply(partitions_list[[as.character(step)]], function(x) {
      x$freq
    }),
    ncores = ncores)

    # update the returned object with a list containing the ecc, an UMAP embedding,
    # the most frequent partition and the number of different partitions
    return_list[[object_name]][[as.character(step)]] = list(
      ecc = ec_consistency,
      embedding = umap_embedding,
      most_frequent_partition = partitions_list[[as.character(step)]][[1]],
      n_different_partitions = length(partitions_list[[as.character(step)]])
    )
  }

  return_list
}


#' Feature Stability Boxplot
#'
#' @description Display, for each feature set and for each step, the distribution
#' of the ECC as a boxplot. Above each boxplot there will a number representing
#' the step (or the size of the subset)
#'
#' @param feature_object_list an object or a concatenation of objects returned by the
#' `get_feature_stability_object` method
#'
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' feature_stability_obj = get_feature_stability_object(data_matrix = as.matrix(mtcars),
#'    feature_set = colnames(mtcars),
#'    steps = -1,
#'    npcs = 5,
#'    n_repetitions = 1)
#' plot_feature_stability_boxplot(feature_stability_obj)
plot_feature_stability_boxplot = function(feature_object_list) {
  min_index = -1 # indicates the number of steps that will be displayed on the plot

  # create a dataframe based on the object returned by `get_feature_stability_object`
  for(config_name in names(feature_object_list)) {
    melt_object = reshape2::melt(lapply(feature_object_list[[config_name]], function(x) {
      x$ecc
    }))
    melt_object$L1 = factor(melt_object$L1, levels = stringr::str_sort(unique(melt_object$L1), numeric = T))
    melt_object[["feature_set"]] = rep(config_name, nrow(melt_object))

    temp_steps_df = data.frame(step = levels(melt_object$L1),
                               index = 1:nlevels(melt_object$L1))
    temp_steps_df$feature_set = rep(config_name, nrow(temp_steps_df))

    levels(melt_object$L1) = 1:nlevels(melt_object$L1)


    if (min_index == -1) {
      final_steps_df = temp_steps_df
      final_melt_df = melt_object

      min_index = nrow(temp_steps_df)
    } else {
      if (nrow(temp_steps_df) < min_index) {
        min_index = nrow(temp_steps_df)
      }

      final_steps_df = rbind(final_steps_df, temp_steps_df)
      final_melt_df = rbind(final_melt_df, melt_object)
    }
  }

  # given that the input object can have multiple configurations with different
  # number of steps, we will use only the first `min_index` steps
  final_melt_df = final_melt_df %>% dplyr::filter(as.numeric(.data$L1) <= min_index)
  final_steps_df = final_steps_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)

  final_melt_df$feature_set = factor(final_melt_df$feature_set, levels = names(feature_object_list))
  final_steps_df$feature_set = factor(temp_steps_df$feature_set, levels = names(feature_object_list))

  # generate the coordinates where the sizes of the steps will be displayed
  text_position <-
    stats::aggregate(value ~ L1 + feature_set , final_melt_df, max)

  # return the ggplot object
  ggplot2::ggplot(final_melt_df,
                  ggplot2::aes(
                    x = .data$L1,
                    y = .data$value,
                    fill = .data$feature_set
                  )) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.7),
                          width = 0.5) +
    ggplot2::geom_text(
      data = text_position,
      ggplot2::aes(label = final_steps_df$step),
      position = ggplot2::position_dodge(width = 0.7),
      vjust = -0.5,
      size = 4
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::xlab("# features") +
    ggplot2::ylab("EC consistency")
}


#' Feature Stability - Cluster Membership facet
#'
#' @description Display a facet of plots where each graph is associated with a step of
#' a feature set and illustrates the distribution of the most frequent partition over the UMAP
#' embedding.
#'
#' @param feature_object_list an object or a concatenation of objects returned by the
#' `get_feature_stability_object` method
#' @param text_size the size of the cluster label
#'
#'
#' @return A ggplot facet object
#' @export
#'
#' @examples
#' feature_stability_obj = get_feature_stability_object(data_matrix = as.matrix(mtcars),
#'    feature_set = colnames(mtcars),
#'    steps = c(7, 9),
#'    npcs = 4,
#'    n_repetitions = 1)
#' plot_feature_stability_mb_facet(feature_stability_obj)

plot_feature_stability_mb_facet = function(feature_object_list, text_size = 5) {
  first_temp = T

  for (config_name in names(feature_object_list)) {
    for (steps in names(feature_object_list[[config_name]])) {
      temp_df = data.frame(
        x = feature_object_list[[config_name]][[steps]]$embedding[, 1],
        y = feature_object_list[[config_name]][[steps]]$embedding[, 2],
        mb = feature_object_list[[config_name]][[steps]]$most_frequent_partition$mb
      )

      temp_df[["config_name"]] = rep(config_name, nrow(temp_df))
      temp_df[["steps"]] = rep(steps, nrow(temp_df))

      if (first_temp) {
        first_temp = F
        final_df = temp_df
      } else {
        final_df = rbind(final_df, temp_df)
      }


    }
  }

  final_df$steps = factor(final_df$steps, levels = stringr::str_sort(unique(final_df$steps), numeric = T))
  final_df$mb = factor(final_df$mb)
  text.labs = final_df %>% dplyr::group_by(.data$config_name, .data$steps, .data$mb) %>% dplyr::summarise(mean_x = stats::median(.data$x),
                                                                                                          mean_y = stats::median(.data$y))


  ggplot2::ggplot(final_df,
                  ggplot2::aes(
                    x = .data$x,
                    y = .data$y,
                    color = .data$mb
                  )) +
    ggplot2::geom_point(size = 0.3) +
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
    ggplot2::facet_wrap(~ config_name + steps)
}

#' Feature Stability - ECC facet
#'
#' @description Display a facet of plots where each graph is associated with a step of
#' a feature set and illustrates the distribution of the ECC score over the UMAP
#' embedding.
#'
#' @param feature_object_list an object or a concatenation of objects returned by the
#' `get_feature_stability_object` method
#'
#'
#' @return A ggplot facet object
#' @export
#'
#' @examples
#' feature_stability_obj = get_feature_stability_object(data_matrix = as.matrix(mtcars),
#'    feature_set = colnames(mtcars),
#'    steps = c(7, 9),
#'    npcs = 4,
#'    n_repetitions = 1)
#' plot_feature_stability_ecs_facet(feature_stability_obj)

plot_feature_stability_ecs_facet = function(feature_object_list) {
  first_temp = T

  for (config_name in names(feature_object_list)) {
    for (steps in names(feature_object_list[[config_name]])) {
      temp_df = data.frame(
        x = feature_object_list[[config_name]][[steps]]$embedding[, 1],
        y = feature_object_list[[config_name]][[steps]]$embedding[, 2],
        ecc = feature_object_list[[config_name]][[steps]]$ecc
      )

      temp_df[["config_name"]] = rep(config_name, nrow(temp_df))
      temp_df[["steps"]] = rep(steps, nrow(temp_df))

      if (first_temp) {
        first_temp = F
        final_df = temp_df
      } else {
        final_df = rbind(final_df, temp_df)
      }


    }
  }

  final_df$steps = factor(final_df$steps, levels = stringr::str_sort(unique(final_df$steps), numeric = T))


  ggplot2::ggplot(final_df,
                  ggplot2::aes(
                    x = .data$x,
                    y = .data$y,
                    color = .data$ecc
                  )) +
    ggplot2::geom_point(size = 0.3) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(color = "ECC",
                  x = "UMAP_1",
                  y = "UMAP_2") +
    ggplot2::facet_wrap(~ config_name + steps)
}


#' Feature Stability Incremental Boxplot
#'
#' @description Perform an incremental ECS between two consecutive steps. The
#' ECS distribution will be displayed as a boxplot. Above each boxplot there will be
#' a pair of numbers representing the two steps that are compared.
#'
#' @param feature_object_list an object or a concatenation of objects returned by the
#' `get_feature_stability_object` method
#' @param dodge_width used for adjusting the horizontal position of the boxplot; the value
#' will be passed in the `width` argument of the `ggplot2::position_dodge` method
#'
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' feature_stability_obj = get_feature_stability_object(data_matrix = as.matrix(mtcars),
#'    feature_set = colnames(mtcars),
#'    steps = c(7, 9),
#'    npcs = 4,
#'    n_repetitions = 1)
#' plot_feature_stability_ecs_incremental(feature_stability_obj, 0.7)

plot_feature_stability_ecs_incremental = function(feature_object_list, dodge_width = 0.7) {
  min_index = -1

  first_df = T

  for (config_name in names(feature_object_list)) {
    nsteps = length(feature_object_list[[config_name]]) - 1

    # treat the case with only one step

    for (i in 1:nsteps) {
      temp_df = data.frame(
        ecs = element_sim_elscore(
          feature_object_list[[config_name]][[i]]$most_frequent_partition$mb,
          feature_object_list[[config_name]][[i +
                                                1]]$most_frequent_partition$mb
        ),
        index = i,
        feature_set = config_name
      )

      if (first_df) {
        first_df = F

        steps_df = data.frame(
          step = paste(
            names(feature_object_list[[config_name]])[i],
            names(feature_object_list[[config_name]])[i +
                                                        1],
            sep = "-"
          ),
          index = i,
          feature_set = config_name
        )

        ecs_df = temp_df
      } else {
        ecs_df = rbind(ecs_df, temp_df)

        steps_df = rbind(steps_df, c(paste(
          names(feature_object_list[[config_name]])[i],
          names(feature_object_list[[config_name]])[i +
                                                      1],
          sep = "-"
        ),
        i,
        config_name))
      }

    }

    if (min_index == -1 || nsteps < min_index) {
      min_index = nsteps
    }
  }

  ecs_df = ecs_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)
  steps_df = steps_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)

  ecs_df$feature_set = factor(ecs_df$feature_set, levels = names(feature_object_list))
  steps_df$feature_set = factor(steps_df$feature_set, levels = names(feature_object_list))

  ecs_df$index = factor(ecs_df$index)

  text_position <-
    stats::aggregate(ecs ~ index + feature_set , ecs_df, max)


  ggplot2::ggplot(ecs_df,
                  ggplot2::aes(
                    x = .data$index,
                    y = .data$ecs,
                    fill = .data$feature_set
                  )) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = dodge_width)) +
    ggplot2::geom_text(
      data = text_position,
      ggplot2::aes(label = steps_df$step),
      position = ggplot2::position_dodge(width = dodge_width),
      vjust = -0.5,
      size = 4
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::xlab("# features") +
    ggplot2::ylab("EC similiarity")
}


######################## Graph construction ######################################

#### connected components ####

#' Nearest Neighbors - Connected Components Link
#'
#' @description One of the steps in the PhenoGraph pipeline is building a graph
#' by using the kNN method on a reduced-space embedding. This method returns an
#' object that describes the relationship between different number of nearest
#' neighbors and the connectivity of the graph. In the context of graph clustering,
#' the number of connected components can be used as a
#' lower bound for the number of clusters. The calculations are performed multiple
#' times by changing the seed at each repetition.
#'
#' @param object if the graph reduction type is PCA, the object should be a data matrix; in
#' the case of UMAP, the user could also provide a PCA embedding.
#' @param n_neigh_sequence a sequence of the number of nearest neighbors.
#' @param config_name user specified string that uniquely describes the embedding characteristics.
#' @param n_repetitions the number of repetitions of applying the pipeline with different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence a custom seed sequence; if the value is NULL, the sequence will be built starting from 1 with a step of 100.
#' @param graph_reduction_type the graph reduction type, denoting if the graph should be built on either the PCA or the UMAP embedding.
#' @param ncores the number of parallel R instances that will run the code. If the value is set to 1, the code will be run sequentially.
#' @param ... additional arguments passed to the `irlba::irlba` or the `uwot::umap` method, depending on the value of graph_reduction_type.
#'
#'
#'
#' @return A list having one field associated with a number of nearest neighbors.
#' Each value contains an array of the number of connected components
#' obtained on the specified number of repetitions.
#'
#' @export
#'
#' @examples
#' get_nn_conn_comps(object = as.matrix(mtcars),
#'     n_neigh_sequence = c(2,5,15),
#'     config_name = "mt_cars",
#'     n_repetitions = 1,
#'     graph_reduction_type = "UMAP",
#'     min_dist = 0.3,
#'     n_neighbors = 30,
#'     metric = "cosine")

get_nn_conn_comps = function(object,
                             n_neigh_sequence,
                             config_name = "",
                             n_repetitions = 30,
                             seed_sequence = NULL,
                             graph_reduction_type = "UMAP",
                             ncores = 1,
                             ...) {
  # create a seed sequence if it's not provided
  if (is.null(seed_sequence)) {
    seed_sequence = seq(from = 1,
                        by = 100,
                        length.out = n_repetitions)
  }

  ncores = min(ncores, length(seed_sequence))

  # store the additional arguments used by umap or irlba in a list
  suppl_args = list(...)
  for(i in 1:length(suppl_args)) {
    assign(names(suppl_args)[i], suppl_args[[i]])
  }
  suppl_arg_names = names(suppl_args)

  nn_conn_comps_list = list()

  if(ncores > 1) {
    # create a parallel backend
    my_cluster <- parallel::makeCluster(
      ncores,
      type = "PSOCK"
    )

    doParallel::registerDoParallel(cl = my_cluster)
  } else {
    # create a sequential backend
    foreach::registerDoSEQ()
  }

  all_vars = ls()
  # the variables needed in each PSOCK process
  needed_vars = c("object", "n_neigh_sequence", "graph_reduction_type", "suppl_arg_names")

  # send the name of the dim reduction arguments
  nn_conn_comps_list_temp = foreach::foreach(seed = seed_sequence,
                                             .noexport = all_vars[!(all_vars %in% needed_vars)],
                                             .export = names(suppl_args)) %dopar% { #
                                               dim_red_args = list()

                                               for(suppl_arg_name in suppl_arg_names) {
                                                 dim_red_args[[suppl_arg_name]] = get(suppl_arg_name)
                                               }

                                               # perform the UMAP / PCA dimensionality reduction
                                               set.seed(seed)
                                               if (graph_reduction_type == "UMAP") {
                                                 embedding = do.call(uwot::umap, c(list(X=object), dim_red_args))
                                                 colnames(embedding) = paste0("UMAP_", 1:ncol(embedding))
                                               } else {
                                                 irlba_object = do.call(irlba::irlba, c(list(A=object), dim_red_args))
                                                 embedding = irlba_object$u %*% diag(irlba_object$d)
                                                 colnames(embedding) = paste0("PC_", 1:ncol(embedding))
                                               }
                                               rownames(embedding) = rownames(object)

                                               # for each neighbour, return the number of connected components
                                               # of the generated graph
                                               sapply(n_neigh_sequence, function(n_neigh) {
                                                 g = igraph::graph_from_adjacency_matrix(
                                                   Seurat::FindNeighbors(
                                                     embedding,
                                                     k.param = n_neigh,
                                                     nn.method = "rann",
                                                     compute.SNN = F,
                                                     verbose = F
                                                   )$nn
                                                 )

                                                 length(unique(igraph::clusters(g)$membership))
                                               })
                                             }

  # if a parallel backend was created, terminate the processes
  if(ncores > 1)
    parallel::stopCluster(cl = my_cluster)

  # store the results obtained for each number of neighbours in different lists
  for (i in 1:length(n_neigh_sequence)) {
    n_neigh = n_neigh_sequence[i]
    nn_conn_comps_list[[as.character(n_neigh)]] = sapply(nn_conn_comps_list_temp, function(x) { x[i]})
  }

  nn_conn_comps_list = list(nn_conn_comps_list)
  names(nn_conn_comps_list) = c(paste(config_name, graph_reduction_type, sep = "_"))

  nn_conn_comps_list
}

#' Graphical representation of the link between NN and graph connectivity
#'
#' @description Display, for each number of neighbors, the distribution of the
#' number connected components obtained at different seeds. The distribution is
#' displayed as a boxplot.
#'
#' @param nn_conn_comps_object an object or a concatenation of objects returned by the
#' `get_nn_conn_comps` method
#'
#'
#' @return A ggplot object
#' @export
#'
#' @note The number of connected components is displayed on a logarithmic scale.
#'
#' @examples
#' nn_conn_comps_obj = get_nn_conn_comps(object = as.matrix(mtcars),
#'     n_neigh_sequence = c(2,5,15),
#'     config_name = "mt_cars",
#'     n_repetitions = 1,
#'     graph_reduction_type = "UMAP",
#'     min_dist = 0.3,
#'     n_neighbors = 30,
#'     metric = "cosine")
#' plot_connected_comps_evolution(nn_conn_comps_obj)

plot_connected_comps_evolution = function(nn_conn_comps_object) {
  final_comps_df = reshape2::melt(nn_conn_comps_object)
  colnames(final_comps_df) = c("n_comps", "n_neigh", "config_name")

  final_comps_df$config_name = factor(final_comps_df$config_name)
  final_comps_df$n_neigh = factor(final_comps_df$n_neigh)
  final_comps_df$n_neigh = factor(final_comps_df$n_neigh, levels(final_comps_df$n_neigh)[stringr::str_order(levels(final_comps_df$n_neigh), numeric = T)])

  max_y = max(final_comps_df$n_comps)
  max_y = max_y + ifelse(max_y %% 5, 5 - max_y %% 5, 0)
  min_y = min(final_comps_df$n_comps)
  min_y = min_y - min_y %% 5


  chosen_breaks = generate_breaks(min_y, max_y)



  ggplot2::ggplot(
    final_comps_df,
    ggplot2::aes(
      x = .data$n_neigh,
      y = .data$n_comps,
      fill = .data$config_name
    )
  ) +
    ggplot2::geom_hline(yintercept = chosen_breaks,
                        linetype = "dashed",
                        color = "#bcc4be") +
    ggplot2::geom_boxplot() +
    ggplot2::scale_y_continuous(breaks = chosen_breaks, trans = "log10")  +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 1,
      hjust = 1
    )) +
    ggplot2::labs(x = "# of nearest neighbors",
                  y = "# of connected components",
                  fill = "configuration") +
    ggplot2::ggtitle("Distribution of the number of connected components")
}

#### number of neigh, graph type importance ####

#' Graph building Parameters Importance
#'
#' @description Creates an object that contain stability measurements and assessments
#' when changing the values of different parameters involved in the Graph Building step,
#' namely the base embedding, the graph type and the number of neighbours.
#'
#' @param object if the graph reduction type is PCA, the object should be a normalized data matrix; in
#' the case of UMAP, the user could also provide a matrix representing the PCA embedding.
#' @param n_neigh_sequence a sequence of the number of nearest neighbours.
#' @param n_repetitions the number of repetitions of applying the pipeline with different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence a custom seed sequence; if the value is NULL, the sequence will be built starting from 1 with a step of 100.
#' @param graph_reduction_type the graph reduction type, denoting if the graph should be built on either the PCA or the UMAP embedding.
#' @param ecs_thresh the ECS threshold used for merging similar clusterings.
#' @param ncores the number of parallel R instances that will run the code. If the value is set to 1, the code will be run sequentially.
#' @param transpose decides whether the input object will be transposed or not
#' @param graph_type argument indicating whether the graph should be
#' unweighted (0), weighted (1) or both (2).
#' @param ... additional arguments passed to the `irlba::irlba` or the `uwot::umap` method, depending on the value of graph_reduction_type.
#'
#'
#'
#' @return A list having four fields:
#'
#' * n_neigh_k_corresp - list containing the number of the clusters obtained by running
#' the pipeline multiple times with different seed, number of neighbors and graph type (weighted vs unweigted)
#' * n_neigh_ec_consistency - list containing the EC consistency of the partitions obtained
#' at multiple runs when changing the number of neighbors or the graph type
#' * partitions_list - the merged list of partitions obtained on the specified number of repetitions
#' * graph_reduction_type - a string representing the graph reduction type
#'
#' @md
#' @export
#'
#' @examples
#' get_nn_importance(object = as.matrix(mtcars),
#'     n_neigh_sequence = c(2,5,15),
#'     n_repetitions = 1,
#'     graph_reduction_type = "UMAP",
#'     min_dist = 0.3,
#'     n_neighbors = 30,
#'     metric = "cosine")

get_nn_importance = function(object,
                             n_neigh_sequence,
                             n_repetitions = 30,
                             seed_sequence = NULL,
                             graph_reduction_type = "PCA",
                             ecs_thresh = 0.99,
                             ncores = 1,
                             transpose = (graph_reduction_type == "PCA"),
                             graph_type = 2,
                             ...) {
  partitions_list = list()

  if(graph_type != 1)
    partitions_list[[paste(graph_reduction_type, "nn", ecs_thresh, sep = "_")]] = list()

  if(graph_type != 0)
    partitions_list[[paste(graph_reduction_type, "snn", ecs_thresh, sep = "_")]] = list()

  # transpose the matrix if the parameter is set to TRUE
  if(transpose)
    object = t(object)

  # create a seed sequence if it's not provided
  if (is.null(seed_sequence)) {
    seed_sequence = seq(from = 1,
                        by = 100,
                        length.out = n_repetitions)
  }

  ncores = min(ncores, length(seed_sequence))

  # additional arguments that will be used by umap or irlba
  suppl_args = list(...)
  for(i in 1:length(suppl_args)) {
    assign(names(suppl_args)[i], suppl_args[[i]])
  }
  suppl_arg_names = names(suppl_args)

  for (n_neigh in n_neigh_sequence) {
    partitions_list[[paste(graph_reduction_type, "snn", ecs_thresh, sep = "_")]][[as.character(n_neigh)]] = list()

    if(ncores > 1) {
      # create a parallel backend
      my_cluster <- parallel::makeCluster(
        ncores,
        type = "PSOCK"
      )

      doParallel::registerDoParallel(cl = my_cluster)
    } else {
      # create a sequential backend
      foreach::registerDoSEQ()
    }

    # the variables needed in the PSOCK processes
    needed_vars = c("object", "n_neigh", "graph_reduction_type", "suppl_arg_names", "graph_type")
    all_vars = ls()

    partitions_list_temp = foreach::foreach(seed = seed_sequence,
                                            .noexport = all_vars[!(all_vars %in% needed_vars)],
                                            .export = names(suppl_args)) %dopar% { #
                                              dim_red_args = list()

                                              for(suppl_arg_name in suppl_arg_names) {
                                                dim_red_args[[suppl_arg_name]] = get(suppl_arg_name)
                                              }

                                              # perform the dimensionality reduction
                                              set.seed(seed)
                                              if (graph_reduction_type == "UMAP") {
                                                embedding = do.call(uwot::umap, c(list(X=object), dim_red_args))
                                                colnames(embedding) = paste0("UMAP_", 1:ncol(embedding))
                                              } else {
                                                irlba_object = do.call(irlba::irlba, c(list(A=object), dim_red_args))
                                                embedding = irlba_object$u %*% diag(irlba_object$d)
                                                colnames(embedding) = paste0("PC_", 1:ncol(embedding))
                                              }
                                              rownames(embedding) = rownames(object)

                                              # build the nn and snn graphs
                                              neigh_matrix = Seurat::FindNeighbors(
                                                embedding,
                                                k.param = n_neigh,
                                                nn.method = "rann",
                                                verbose = F
                                              )

                                              # apply the Leiden clustering on the graph specified by the variable `graph_type`
                                              if(graph_type != 1) {
                                                cluster_results_nn = Seurat::FindClusters(
                                                  neigh_matrix$nn,
                                                  random.seed = seed,
                                                  algorithm = 4,
                                                  verbose = F
                                                )
                                                cluster_results_nn = list(mb = cluster_results_nn[, names(cluster_results_nn)[1]],
                                                                          freq = 1,
                                                                          seed = seed)

                                                if(graph_type == 0)
                                                  return(cluster_results_nn)
                                              }

                                              cluster_results_snn = Seurat::FindClusters(
                                                neigh_matrix$snn,
                                                random.seed = seed,
                                                algorithm = 4,
                                                verbose = F
                                              )
                                              cluster_results_snn = list(mb = cluster_results_snn[, names(cluster_results_snn)[1]],
                                                                         freq = 1,
                                                                         seed = seed)

                                              if(graph_type == 1)
                                                return(cluster_results_snn)

                                              return(list(cluster_results_nn, cluster_results_snn))
                                            }

    # if a parallel backend was created, terminate the processes
    if(ncores > 1)
      parallel::stopCluster(cl = my_cluster)

    # merge the partitions that are considered similar by a given ecs threshold
    for(i in 1:length(partitions_list)) {
      partitions_list[[i]][[as.character((n_neigh))]] = merge_partitions(lapply(partitions_list_temp, function(x) { x[[i]]}),
                                                                         ecs_thresh = ecs_thresh,
                                                                         ncores = ncores)
    }
  }

  # create an object showing the number of clusters obtained for each number of
  # neighbours
  nn_object_n_clusters = list()
  for (config_name in names(partitions_list)) {
    nn_object_n_clusters[[config_name]] = list()
    for (n_neigh in names(partitions_list[[config_name]])) {
      nn_object_n_clusters[[config_name]][[n_neigh]] = unlist(lapply(partitions_list[[config_name]][[n_neigh]], function(x) {
        rep(length(unique(x$mb)), x$freq)
      }))
    }
  }

  # create an object that contains the ecc of the partitions obtained for each
  # number of neighbours
  nn_ecs_object = lapply(partitions_list, function(config) {
    lapply(config, function(n_neigh) {
      weighted_element_consistency(lapply(n_neigh, function(x) {
        x$mb
      }),
      sapply(n_neigh, function(x) {
        x$freq
      }),
      ncores = ncores)
    })
  })

  list(
    n_neigh_k_corresp = nn_object_n_clusters,
    n_neigh_ec_consistency = nn_ecs_object,
    n_different_partitions = lapply(partitions_list, function(config) {
      sapply(config, function(n_neigh) { length(n_neigh)})
    }),
    graph_reduction_type = graph_reduction_type
  )
}


# 2. neighbors <-> number of clusters association


#' Graphical representation of the link between NN and number of clusters
#'
#' @description Display, for each number of neighbors, the distribution of the
#' number of clusters obtained at different seeds. The distribution is
#' displayed as a boxplot.
#'
#' @param nn_object_n_clusters an object or a concatenation of objects returned by the
#' `get_nn_importance` method
#'
#'
#' @return A ggplot object
#' @export
#'
#' @note The number of clusters is displayed on a logarithmic scale.
#'
#' @examples
#' nn_importance_obj = get_nn_importance(object = as.matrix(mtcars),
#'     n_neigh_sequence = c(2,5,15),
#'     n_repetitions = 1,
#'     graph_reduction_type = "UMAP",
#'     min_dist = 0.3,
#'     n_neighbors = 30,
#'     metric = "cosine")
#' plot_n_neigh_k_correspondence(nn_importance_obj)

plot_n_neigh_k_correspondence = function(nn_object_n_clusters) {
  melted_obj = reshape2::melt(nn_object_n_clusters$n_neigh_k_corresp)
  colnames(melted_obj) = c("k", "n_neigh", "config_name")


  melted_obj$n_neigh = factor(melted_obj$n_neigh)

  melted_obj$n_neigh = factor(melted_obj$n_neigh, levels(melted_obj$n_neigh)[stringr::str_order(levels(melted_obj$n_neigh), numeric = T)])

  max_y = max(melted_obj$k)
  max_y = max_y + ifelse(max_y %% 5, 5 - max_y %% 5, 0)
  min_y = min(melted_obj$k)
  min_y = min_y - min_y %% 5


  chosen_breaks = generate_breaks(min_y, max_y)

  ggplot2::ggplot(melted_obj,
                  ggplot2::aes(
                    x = .data$n_neigh,
                    y = .data$k,
                    fill = .data$config_name
                  )) +
    ggplot2::geom_hline(yintercept = chosen_breaks,
                        linetype = "dashed",
                        color = "#bcc4be") +
    ggplot2::geom_boxplot() +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(breaks = chosen_breaks, trans = "log10") +
    ggplot2::labs(x = "# of nearest neighbors",
                  fill = "configuration") +
    ggplot2::ggtitle("Distribution of k across different seeds")
}

# 3. ECS distribution across different seeds for different number of neighbors

#' Graph construction parameteres - ECC facet
#'
#' @description Display, for all configurations consisting in different number
#' of neighbors, graph types and base embeddings, the EC Consistency of the partitions
#' obtained over multiple runs on an UMAP embedding.
#'
#' @param nn_ecs_object an object or a concatenation of objects returned by the
#' `get_nn_importance` method
#'
#'
#' @return A ggplot facet object
#' @export
#'
#'
#' @examples
#' nn_importance_obj = get_nn_importance(object = as.matrix(mtcars),
#'     n_neigh_sequence = c(2,5,15),
#'     n_repetitions = 1,
#'     graph_reduction_type = "UMAP",
#'     min_dist = 0.3,
#'     n_neighbors = 30,
#'     metric = "cosine")
#' plot_n_neigh_ecs(nn_importance_obj)

plot_n_neigh_ecs = function(nn_ecs_object) {
  melted_obj = reshape2::melt(nn_ecs_object$n_neigh_ec_consistency)
  colnames(melted_obj) = c("ECS", "n_neigh", "config_name")


  melted_obj$n_neigh = factor(melted_obj$n_neigh)

  melted_obj$n_neigh = factor(melted_obj$n_neigh, levels(melted_obj$n_neigh)[stringr::str_order(levels(melted_obj$n_neigh), numeric = T)])

  ggplot2::ggplot(melted_obj,
                  ggplot2::aes(
                    x = .data$n_neigh,
                    y = .data$ECS,
                    fill = .data$config_name
                  )) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "# of nearest neighbors",
                  y = "EC consistency",
                  fill = "configuration") +
    ggplot2::ggtitle("Distribution of ECS across different seeds for different # neighbors")
}


############################## Clustering ########################################

#### clustering methods ####

#' Clustering Method Importance
#'
#' @description Analyzes the importance of choosing a specific graph clustering
#' method in the pipeline. The method will iterate through different values of
#' the resolution parameter and compare, using the EC Consistency score, the
#' partitions obtained at different seeds.
#'
#' @param graph_adjacency_matrix a square adjacency matrix based on which an igraph
#' object will be built.
#' @param resolution a sequence of resolution values.
#' @param n_repetitions the number of repetitions of applying the pipeline with different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence a custom seed sequence; if the value is NULL, the sequence will be built starting from 1 with a step of 100.
#' @param ecs_thresh the ECS threshold used for merging similar clusterings.
#' @param ncores the number of parallel R instances that will run the code. If the value is set to 1, the code will be run sequentially.
#' @param verbose boolean value used for displaying the progress bar
#'
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
#'
#'
#' @examples
#' adj_matrix = Seurat::FindNeighbors(as.matrix(mtcars),
#'     k.param = 25,
#'     nn.method = "rann",
#'     verbose = FALSE,
#'     compute.SNN = FALSE)$nn
#' get_clustering_difference_object(graph_adjacency_matrix = adj_matrix,
#'     resolution = 0.5,
#'     n_repetitions = 1,
#'     verbose = FALSE)


get_clustering_difference_object = function(graph_adjacency_matrix,
                                            resolution,
                                            n_repetitions = 30,
                                            seed_sequence = NULL,
                                            ecs_thresh = 1,
                                            ncores = 1,
                                            verbose = TRUE) {
  # create a seed sequence if it's not provided
  if (is.null(seed_sequence)) {
    seed_sequence = seq(from = 1,
                        by = 100,
                        length.out = n_repetitions)

  }

  algorithm_names = c("Louvain", "Louvain.refined", "SLM", "Leiden")
  result_object = list()

  for (i in 1:4) {
    # if verbose is set to TRUE, create a progress bar
    if(verbose) {
      message(paste("generating partitions for", algorithm_names[i], "method"))
      pb <- progress::progress_bar$new(
        format = "resolution = :what [:bar] :current/:total eta: :eta  total elapsed:  :elapsed",
        total = length(resolution),
        clear = FALSE,
        width = 80
      )

      pb$tick(0)

    }

    # calculate the partitions for each resolution value
    result_object[[algorithm_names[i]]] = sapply(resolution, function(res) {
      if(verbose) { pb$tick(tokens = list(what = res)) }
      res_list = list()
      res_list[[as.character(res)]] = get_resolution_partitions(
        graph_adjacency_matrix,
        resolution = res,
        clustering_function = seurat_clustering,
        algorithm = i,
        seed_sequence = seed_sequence,
        ncores = ncores,
        ecs_thresh = ecs_thresh
      )

      res_list
    })

  }

  # for each resolution value, determine the number of clusters that appears most frequently
  # filter the list of partition based on that number of clusters and calculate ecc
  filtered_result = lapply(result_object, function(config_object) {
    lapply(config_object, function(res) {
      corresp_k = names(which.max(unlist(lapply(res, function(res_cluster) {
        sum(sapply(res_cluster, function(partition) {
          partition$freq
        }))
      }))))[1]

      list(ec_consistency = weighted_element_consistency(lapply(res[[corresp_k]], function(x) {
        x$mb
      }),
      sapply(res[[corresp_k]], function(x) {
        x$freq
      }),
      ncores = ncores),
      k = corresp_k)
    })
  })

  # calculate ecc for all the partitions created with a resolution value
  all_result = lapply(result_object, function(config_object) {
    lapply(config_object, function(res) {
      k_names = names(res)
      partition_list = list()
      for (k in k_names) {
        partition_list = c(partition_list, res[[k]])
      }

      weighted_element_consistency(lapply(partition_list, function(x) {
        x$mb
      }),
      sapply(partition_list, function(x) {
        x$freq
      }),
      ncores = ncores)
    })
  })

  list(filtered = filtered_result,
       all = all_result)
}


#' Graphical representation of the clustering method importance
#'
#' @description Display, for each clustering method and each resolution value,
#' the distribution of the EC consistency as a boxplot. We will use the `filtered`
#' field of the object returned by the `get_clustering_difference_object` method.
#' Above each boxplot, the number of clusters will be displayed.
#'
#' @param clustering_difference_object an object returned by the
#' `get_clustering_difference_object` method
#'
#'
#' @return A ggplot object
#' @export
#'
#'
#' @examples
#' adj_matrix = Seurat::FindNeighbors(as.matrix(mtcars),
#'     k.param = 25,
#'     nn.method = "rann",
#'     verbose = FALSE,
#'     compute.SNN = FALSE)$nn
#' clust_diff_obj = get_clustering_difference_object(graph_adjacency_matrix = adj_matrix,
#'     resolution = c(0.5, 0.8, 1),
#'     n_repetitions = 1,
#'     verbose = FALSE)
#' plot_clustering_difference_boxplot(clust_diff_obj)

plot_clustering_difference_boxplot = function(clustering_difference_object) {
  k_labels = reshape2::melt(lapply(clustering_difference_object$filtered, function(config_object) {
    lapply(config_object, function(res) {
      res$k
    })
  })) %>% dplyr::arrange(.data$L1)

  melted_obj = reshape2::melt(lapply(clustering_difference_object$filtered, function(config_object) {
    lapply(config_object, function(res) {
      res$ec_consistency
    })
  }))

  colnames(melted_obj) = c("EC_consistency", "resolution", "clustering_method")


  text_position <-
    stats::aggregate(EC_consistency ~ resolution + clustering_method ,
                     melted_obj,
                     max)

  ggplot2::ggplot(
    melted_obj,
    ggplot2::aes(
      x = .data$resolution,
      y = .data$EC_consistency,
      fill = .data$clustering_method
    )
  ) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::theme_classic() +
    ggplot2::labs(y = "EC Consistency",
                  fill = "clustering method") +
    ggplot2::geom_text(
      data = text_position,
      ggplot2::aes(label = k_labels$value),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = -0.5,
      size = 3
    )
}


#' Graphical representation of the clustering method importance
#'
#' @description Display, for each clustering method and each resolution value,
#' the distribution of the EC consistency on a given embedding. We will use the `all`
#' field of the object returned by the `get_clustering_difference_object` method.
#'
#' @param clustering_difference_object an object returned by the
#' `get_clustering_difference_object` method
#' @param embedding an embedding (only the first two dimensions will be used for
#' visualisation)
#' @param low_limit the lowest value of ECC that will be displayed on the embedding
#' @param high_limit the highest value of ECC that will be displayed on the embedding
#' @param grid boolean value indicating whether the facet should be a grid (where each
#' row is associated with a resolution value and each column with a clustering method) or
#' a wrap
#'
#' @return A ggplot facet object
#' @export
#'
#'
#' @examples
#' adj_matrix = Seurat::FindNeighbors(as.matrix(mtcars),
#'     k.param = 25,
#'     nn.method = "rann",
#'     verbose = FALSE,
#'     compute.SNN = FALSE)$nn
#' clust_diff_obj = get_clustering_difference_object(graph_adjacency_matrix = adj_matrix,
#'     resolution = c(0.5, 0.8, 1),
#'     n_repetitions = 1,
#'     verbose = FALSE)
#' plot_clustering_difference_facet(clust_diff_obj, as.matrix(mtcars[,1:2]))

plot_clustering_difference_facet = function(clustering_difference_object,
                                            embedding,
                                            low_limit = 0,
                                            high_limit = 1,
                                            grid = T) {
  npoints = nrow(embedding)
  if (length(clustering_difference_object$all[[1]][[1]]) != npoints) {
    stop("The provided embedding and the consistency arrays must have the same number of elements!")
  }

  if (ncol(embedding) < 2) {
    stop("The embedding should have at least two dimensions!")
  }



  melt_obj = reshape2::melt(clustering_difference_object$all)

  n_embeeding_repetitions = nrow(melt_obj) / npoints

  melt_obj["x"] = rep(embedding[, 1], n_embeeding_repetitions)
  melt_obj["y"] = rep(embedding[, 2], n_embeeding_repetitions)


  melt_obj["value"][melt_obj["value"] < low_limit] = low_limit
  melt_obj["value"][melt_obj["value"] > high_limit] = high_limit



  return_plot = ggplot2::ggplot(melt_obj,
                                ggplot2::aes(
                                  x = .data$x,
                                  y = .data$y,
                                  color = .data$value
                                )) +
    ggplot2::geom_point(size = 0.2) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(color = "ECC")

  if (grid) {
    return(return_plot + ggplot2::facet_grid(L2 ~ L1))
  }

  return_plot + ggplot2::facet_wrap(~ L2 + L1)
}

#### resolution ####

# using a graph-based clustering algorithm, determine the partitions
# obtained on a graph using different resolution values and different seeds
get_resolution_partitions = function(clustered_object,
                                     resolution,
                                     clustering_function,
                                     seed_sequence,
                                     ncores = 1,
                                     ecs_thresh = 0.99,
                                     ...) {
  different_partitions = list()
  ncores = min(ncores, length(seed_sequence))

  # additional arguments used by the clustering method
  suppl_args = list(...)
  i = 1
  while(i <= length(suppl_args)) {
    assign(names(suppl_args)[i], suppl_args[[i]])
    i = i + 1
  }

  suppl_arg_names = names(suppl_args)

  if(ncores > 1) {
    # create a parallel backend
    my_cluster <- parallel::makeCluster(
      ncores,
      type = "PSOCK"
    )

    doParallel::registerDoParallel(cl = my_cluster)
  } else {
    # create a sequential backend
    foreach::registerDoSEQ()
  }

  # the variables needed in the PSOCK processes
  needed_vars = c("resolution", "clustering_function", "clustered_object", "suppl_arg_names")
  all_vars = ls()

  different_partitions_temp = foreach::foreach(seed = seed_sequence,
                                               .noexport = all_vars[!(all_vars %in% needed_vars)],
                                               .export = names(suppl_args)) %dopar% { #
                                                 clust_args = list()

                                                 seed = get("seed")

                                                 for(suppl_arg_name in suppl_arg_names) {
                                                   clust_args[[suppl_arg_name]] = get(suppl_arg_name)
                                                 }

                                                 # apply the clustering, which should return a membership vector
                                                 do.call(clustering_function, c(list(object = clustered_object,
                                                                                     resolution = resolution,
                                                                                     seed = seed),
                                                                                clust_args))
                                               }

  # if a parallel backend was created, terminate the processes
  if(ncores > 1)
    parallel::stopCluster(cl = my_cluster)

  # group the partitions by the number of clusters
  for(i in 1:length(seed_sequence)) {
    k = as.character(length(unique(different_partitions_temp[[i]])))

    if (!(k %in% names(different_partitions))) {
      different_partitions[[k]] = list()
      different_partitions[[k]][[1]] = list(mb = different_partitions_temp[[i]],
                                            freq = 1,
                                            seed = seed_sequence[i])
    } else {
      index = length(different_partitions[[k]])
      different_partitions[[k]][[index+1]] = list(mb = different_partitions_temp[[i]],
                                                  freq = 1,
                                                  seed = seed_sequence[i])
    }
  }

  # merge the partitions using the ecs threshold
  for(k in names(different_partitions)) {
    different_partitions[[k]] = merge_partitions(different_partitions[[k]],
                                                 ecs_thresh = ecs_thresh,
                                                 ncores = ncores,
                                                 order = TRUE)
  }

  different_partitions
}

#' Resolution gridsearch
#'
#' @description Perform a gridsearch over the resolution, number of neighbors
#' and graph type.
#'
#' @param embedding the base embedding for the graph construction.
#' @param resolution a sequence of resolution values.
#' @param n_neigh a value or a sequence of number of neighbors used for graph construction.
#' @param n_repetitions the number of repetitions of applying the pipeline with different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence a custom seed sequence; if the value is NULL, the sequence will be built starting from 1 with a step of 100.
#' @param clustering_method the graph clustering method that should be applied;
#' the user should provide an integer between 1 and 4, where encodings are the
#' same as specified in the Seurat's `FindClusters` method. The arguments allows
#' selecting multiple algorithms.
#' @param graph_type argument indicating whether the graph should be
#' unweighted (0), weighted (1) or both (2).
#' @param object_name user specified string that uniquely describes the embedding characteristics.
#' @param ecs_thresh the ECS threshold used for merging similar clusterings.
#' @param ncores the number of parallel R instances that will run the code. If the value is set to 1, the code will be run sequentially.
#'
#'
#'
#' @return A five-level list. The hierarchy is as follows:
#'
#' * the configuration name: concatenation between the object name provided by
#' the user, the number of neighbors, the graph type and the clustering method
#' * the resolution value \eqn{\gamma}
#' * the number of clusters *k* that can be obtained using the specified resolution
#' * the partitions obtained with resolution \eqn{\gamma} and have *k* clusters
#' * the structure of a partitions, which consists in having a `mb` field with
#' the flat membership vector, `freq` denoting its frequency and `seed`, that is
#' the seed used to obtain this partition in this configuration.
#'
#' @md
#' @export
#'
#' @examples
#' wrapper_gridsearch(embedding = as.matrix(mtcars),
#'    resolution = 0.8,
#'    n_neigh = 25,
#'    n_repetitions = 1,
#'    clustering_method = 4,
#'    graph_type = 0,
#'    object_name = "mt_cars")

wrapper_gridsearch = function(embedding,
                              resolution,
                              n_neigh,
                              n_repetitions = 100,
                              seed_sequence = NULL,
                              clustering_method = 4,
                              graph_type = 0,
                              object_name = NULL,
                              ecs_thresh = 1,
                              ncores = 1) {
  # another parameters
  # prune ? for FindNeighbors, SNN case
  # num_iter for FindClusters


  different_partitions = list()

  if (class(resolution) != "numeric") {
    stop("resolution parameter should take a numeric value")
  }

  if (!(class(n_neigh) %in% c("numeric", "integer"))) {
    stop("n_neigh parameter should take a numeric value")
  }

  if (!(class(clustering_method) %in% c("numeric", "integer")) ||
      sum(clustering_method %in% 1:4) != length(clustering_method)) {
    stop("clustering_method should be a number between 1 and 4")
  }


  if (class(graph_type) != "numeric" || !(graph_type %in% 0:2)) {
    stop("graph_type should be a number between 0 and 2")
  }

  if (length(clustering_method) > 1) {
    for (i in 1:length(clustering_method)) {
      different_partitions = c(
        different_partitions,
        wrapper_gridsearch(
          embedding = embedding,
          resolution = resolution,
          n_neigh = n_neigh,
          n_repetitions = n_repetitions,
          seed_sequence = seed_sequence,
          clustering_method = clustering_method[i],
          graph_type = graph_type,
          object_name = object_name,
          ncores = ncores,
          ecs_thresh = ecs_thresh
        )
      )

    }

    return(different_partitions)
  }

  if (length(n_neigh) > 1) {
    for (i in 1:length(n_neigh)) {
      different_partitions = c(
        different_partitions,
        wrapper_gridsearch(
          embedding = embedding,
          resolution = resolution,
          n_neigh = n_neigh[i],
          n_repetitions = n_repetitions,
          seed_sequence = seed_sequence,
          clustering_method = clustering_method,
          graph_type = graph_type,
          object_name = object_name,
          ncores = ncores,
          ecs_thresh = ecs_thresh
        )
      )

    }

    return(different_partitions)
  }

  if (graph_type == 2) {
    different_partitions = c(
      wrapper_gridsearch(
        embedding = embedding,
        resolution = resolution,
        n_neigh = n_neigh,
        n_repetitions = n_repetitions,
        seed_sequence = seed_sequence,
        clustering_method = clustering_method,
        graph_type = 0,
        object_name = object_name,
        ecs_thresh = ecs_thresh,
        ncores = ncores
      ),
      wrapper_gridsearch(
        embedding = embedding,
        resolution = resolution,
        n_neigh = n_neigh,
        n_repetitions = n_repetitions,
        seed_sequence = seed_sequence,
        clustering_method = clustering_method,
        graph_type = 1,
        object_name = object_name,
        ecs_thresh = ecs_thresh,
        ncores = ncores
      )
    )

    return(different_partitions)
  }

  graph_type_name = ifelse(graph_type == 0, "NN", "SNN")
  algorithm_name = switch(clustering_method,
                          "Louvain",
                          "Louvain.refined",
                          "SLM",
                          "Leiden")

  details = paste(n_neigh, graph_type_name, algorithm_name, sep = "_")

  if (length(resolution) > 1) {
    for (i in 1:length(resolution)) {
      res = resolution[i]

      if (i == 1) {
        different_partitions = wrapper_gridsearch(
          embedding = embedding,
          resolution = res,
          n_neigh = n_neigh,
          n_repetitions = n_repetitions,
          seed_sequence = seed_sequence,
          clustering_method = clustering_method,
          graph_type = graph_type,
          object_name = object_name,
          ncores = ncores,
          ecs_thresh = ecs_thresh
        )

        used_name = names(different_partitions)[1]
      } else {
        different_partitions[[used_name]] = c(
          different_partitions[[used_name]],
          wrapper_gridsearch(
            embedding = embedding,
            resolution = res,
            n_neigh = n_neigh,
            n_repetitions = n_repetitions,
            seed_sequence = seed_sequence,
            clustering_method = clustering_method,
            graph_type = graph_type,
            object_name = object_name,
            ncores = ncores,
            ecs_thresh = ecs_thresh
          )[[used_name]]
        )
      }
    }

    return(different_partitions)
  }

  if (is.null(object_name)) {
    object_name = details
  } else {
    object_name = paste(object_name, details, sep = "_")
  }


  # create a seed sequence if it's not provided
  if (is.null(seed_sequence)) {
    seed_sequence = seq(from = 1,
                        by = 100,
                        length.out = n_repetitions)

  }



  adj_matrix = Seurat::FindNeighbors(
    embedding,
    k.param = n_neigh,
    nn.method = "rann",
    compute.SNN = (graph_type != 0),
    verbose = F
  )[[ifelse(graph_type == 0, "nn", "snn")]]

  different_partitions = get_resolution_partitions(
    adj_matrix,
    resolution = resolution,
    clustering_function = seurat_clustering,
    algorithm = clustering_method,
    seed_sequence = seed_sequence,
    ncores = ncores,
    ecs_thresh = ecs_thresh
  )


  ret_object = list()
  ret_object[[object_name]] = list()
  ret_object[[object_name]][[as.character(round(resolution, digits = 3))]] = different_partitions


  ret_object
}


#' Merge resolutions
#'
#' @description Merge partitions obtained with different resolution values
#' that have the same number of clusters.
#'
#' @param res_obj A list associated to a configuration field from the object
#' returned by the `wrapper_gridsearch` method.
#' @param ncores the number of parallel R instances that will run the code. If the value is set to 1, the code will be run sequentially.
#'
#'
#' @return A list having one field assigned to each number of clusters. A number
#' of cluster will contain a list of all merged partitions. To avoid duplicates,
#' `merged_partitions` with threshold 1 is applied.
#'
#' @md
#' @export
#'
#'
#' @examples
#' merge_resolutions(list(
#'  "0.8" = list(
#'    "2" = list(
#'      list(mb = c(1,1,1,2),
#'           freq = 2),
#'      list(mb = c(1,1,2,2),
#'           freq = 5)
#'    )
#'  ),
#'  "0.9" = list(
#'    "3" = list(
#'      list(mb = c(1,2,3,3),
#'           freq = 1)
#'    ),
#'    "2" = list(
#'      list(mb = c(2,2,3,3),
#'           freq = 2)
#'    )
#'  )
#'))


merge_resolutions = function(res_obj,
                             ncores = 1) {
  clusters_obj = list()
  indices = list()

  for (res.val in names(res_obj)) {
    for (k in names(res_obj[[res.val]])) {
      if (!(k %in% names(clusters_obj))) {
        clusters_obj[[k]] = res_obj[[res.val]][[k]]
        #indices[[k]] = length(clusters_obj[[k]])

      } else {
        clusters_obj[[k]] = c(clusters_obj[[k]], res_obj[[res.val]][[k]])
      }
    }
  }

  ncores = min(ncores, length(clusters_obj))
  k_vals = names(clusters_obj)

  needed_vars = c("clusters_obj")

  if(ncores > 1) {
    my_cluster <- parallel::makeCluster(
      ncores,
      type = "PSOCK"
    )

    doParallel::registerDoParallel(cl = my_cluster)
  } else {
    foreach::registerDoSEQ()
  }


  all_vars = ls()

  clusters_obj = foreach::foreach(i = 1:length(clusters_obj),
                                  .noexport = all_vars[!(all_vars %in% needed_vars)],
                                  .export = c("merge_identical_partitions", "are_identical_memberships")) %dopar% {
                                    merge_identical_partitions(clusters_obj[[i]])
                                  }

  names(clusters_obj) = k_vals

  if(ncores > 1)
    parallel::stopCluster(cl = my_cluster)

  clusters_obj
}


#' Correspondence between the resolution parameter and the number of clusters
#'
#' @description For each configuration provided in the object, display what
#' number of clusters appear for different values of the resolution parameters.
#' Different shapes of points indicate different configuration, while the color
#' illustrates the frequency of the most common partition. The frequency is calculated
#' as the fraction between the number of total appearances and the number of
#' runs.
#'
#' @param res_object_list an object or a concatenation of objects returned by the
#' `wrapper_gridsearch` method
#' @param res_object_names custom names that the user could assing to each
#' configuration; if not specified, the plot will use the generated configuration
#' names
#' @param given_height used for adjusting the vertical position of the boxplot; the value
#' will be passed in the `height` argument of the `ggstance::position_dodgev` method
#'
#'
#' @return A ggplot object
#' @export
#'
#'
#' @examples
#' gridsearch_result = wrapper_gridsearch(embedding = as.matrix(mtcars),
#'    resolution = 0.8,
#'    n_neigh = 25,
#'    n_repetitions = 1,
#'    clustering_method = 4,
#'    graph_type = 0,
#'    object_name = "mt_cars")

#' plot_k_resolution_corresp(gridsearch_result)

plot_k_resolution_corresp = function(res_object_list,
                                     res_object_names = NULL,
                                     given_height = 0.7) {
  if (is.null(res_object_names)) {
    res_object_names = names(res_object_list)
  }

  object_number = length(res_object_list)

  for (i in 1:object_number) {
    res_object = res_object_list[[i]]

    n.runs = sum(sapply(res_object[[names(res_object)[1]]], function(x) {
      sum(sapply(x, function(y) {
        y$freq
      }))
    }))

    list.appereances = lapply(res_object, function(x) {
      lapply(x, function(y) {
        y[[1]]$freq
      })
    })
    temp.df.appereances = reshape2::melt(list.appereances)
    colnames(temp.df.appereances) = c("freq", "number_clusters", "resolution_value")


    temp.df.appereances[["configuration"]] = rep(res_object_names[i], nrow(temp.df.appereances))
    # temp.df.appereances$res.value = format(as.double(temp.df.appereances$res.value), digits = digits)
    temp.df.appereances$freq = temp.df.appereances$freq / n.runs

    if (i == 1) {
      final.df = temp.df.appereances
    } else {
      final.df = rbind(final.df, temp.df.appereances)
    }
  }

  final.df[["configuration"]] = factor(final.df[["configuration"]])
  final.df[["number_clusters"]] = factor(as.numeric(final.df[["number_clusters"]]))


  ggplot2::ggplot(
    final.df,
    ggplot2::aes(
      y = .data$number_clusters,
      x = .data$resolution_value,
      color = .data$freq,
      shape = .data$configuration,
      group = interaction(.data$number_clusters, .data$configuration)
    )
  ) +
    ggplot2::geom_hline(
      yintercept = unique(final.df$number_clusters),
      linetype = "dashed",
      color = "#e3e3e3"
    ) +
    ggplot2::geom_vline(
      xintercept = unique(final.df$resolution_value),
      linetype = "dashed",
      color = "#e3e3e3"
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = given_height),  #ggstance::position_dodgev(height = given.height),
                        size = 2.5) + #position_jitterdodge(jitter.width=0.85))
    ggplot2::theme_classic() +
    #ggplot2::theme(axis.title.y = ggtext::element_markdown()) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = "resolution",
                  y = "k",
                  color = "frequency") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))

}


#' Reproducibility display of the number of clusters
#'
#' @description For each configuration provided in the object, display how
#' many different partitions with the same number of clusters can be obtained
#' by changing the seed. The color gradient suggests the frequency of the most
#' common partition relative to the total number of appearances of that specific
#' number of clusters.
#'
#' @param partition_obj_list an object or a concatenation of objects returned by the
#' `merge_resolutions` method
#' @param object_names custom names that the user could assing to each
#' configuration; if not specified, the plot will use the generated configuration
#' names
#'
#' @return A ggplot object
#' @export
#'
#'
#' @examples
#' gridsearch_result = wrapper_gridsearch(embedding = as.matrix(mtcars),
#'    resolution = 0.8,
#'    n_neigh = 25,
#'    n_repetitions = 1,
#'    clustering_method = 4,
#'    graph_type = 0,
#'    object_name = "mt_cars")

#' merged_resolutions = merge_resolutions(gridsearch_result[[1]])
#' plot_k_n_partitions(list(merged_resolutions), "mt_cars_merged")

plot_k_n_partitions = function(partition_obj_list, object_names = NULL) {
  if (is.null(object_names)) {
    object_names = names(partition_obj_list)
  }

  n.objects = length(partition_obj_list)

  max_n_part = 0

  for (i in 1:n.objects) {
    partition_object = partition_obj_list[[i]]

    unique.partitions.temp.df = reshape2::melt(lapply(partition_object, function(x) {
      length(x)
    }))
    colnames(unique.partitions.temp.df) = c("n.partitions", "n.clusters")

    unique.partitions.temp.df[["configuration"]] = rep(object_names[i], nrow(unique.partitions.temp.df))

    unique.partitions.temp.df[["first.occ"]] = as.numeric(lapply(partition_object, function(x) {
      max(sapply(x, function(y) {
        y$freq
      }))
    }))
    unique.partitions.temp.df[["total.occ"]] = as.numeric(lapply(partition_object, function(x) {
      sum(sapply(x, function(y) {
        y$freq
      }))
    }))
    unique.partitions.temp.df[["frequency"]] = unique.partitions.temp.df$first.occ / unique.partitions.temp.df$total.occ

    max_n_part = max(c(max(
      unique.partitions.temp.df$n.partitions
    ), max_n_part))

    if (i == 1) {
      unique.partitions.final.df = unique.partitions.temp.df
    } else {
      unique.partitions.final.df = rbind(unique.partitions.final.df, unique.partitions.temp.df)
    }
  }

  unique.partitions.final.df$configuration = factor(unique.partitions.final.df$configuration)
  unique.partitions.final.df$n.clusters = factor(unique.partitions.final.df$n.clusters,
                                                 levels = stringr::str_sort(unique(
                                                   unique.partitions.final.df$n.clusters
                                                 ), numeric = T))


  ggplot2::ggplot(
    unique.partitions.final.df,
    ggplot2::aes(
      x = .data$n.clusters,
      y = .data$n.partitions,
      shape = .data$configuration,
      color = .data$frequency
    ),
    group = interaction(.data$n.clusters, .data$precision)
  ) +
    ggplot2::scale_y_continuous(breaks = seq(from = 0, to = max_n_part, by = 2)) +
    ggplot2::geom_hline(
      yintercept = seq(from = 0, to = max_n_part, by = 2),
      linetype = "dashed",
      color = "#e3e3e3"
    ) +
    ggplot2::geom_vline(
      xintercept = unique(unique.partitions.final.df$n.clusters),
      linetype = "dashed",
      color = "#e3e3e3"
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(0.9), size = 2.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_viridis_c() +
    ggplot2::xlab("k") +
    ggplot2::ylab("# partitions")
}
