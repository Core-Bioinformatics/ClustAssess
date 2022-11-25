#' @importFrom foreach %dopar%
NULL

algorithm_names <- c("Louvain", "Louvain.refined", "SLM", "Leiden")
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

# generate values from an interval that simulates the logarithmic scale
generate_breaks <- function(min.range, max.range) {
  start.point <- min.range + 5 - min.range %% 5
  breaks.list <- list()
  index <- 1
  step <- 5
  while (start.point < max.range) {
    breaks.list[[index]] <- start.point
    index <- index + 1

    start.point <- start.point + step

    # the step is doubled after every two iterations
    if (index %% 2) {
      step <- step * 2
    }
  }

  c(min.range, unlist(breaks.list), max.range)
}

#### PCA ####

#' Evaluate Feature Set Stability
#'
#' @description Evaluate the stability of clusterings obtained
#' based on incremental subsets of a given feature set.
#'
#' @param data_matrix A data matrix having the features on the rows
#' and the observations on the columns.
#' @param feature_set A set of feature names that can be found on the rownames of
#' the data matrix.
#' @param steps Vector containing the sizes of the subsets; negative values will
#' be interpreted as using all features.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence A custom seed sequence; if the value is NULL, the
#' sequence will be built starting from 1 with a step of 100.
#' @param feature_type A name associated to the feature_set.
#' @param graph_reduction_type The graph reduction type, denoting if the graph
#' should be built on either the PCA or the UMAP embedding.
#' @param npcs The number of principal components.
#' @param ecs_thresh The ECS threshold used for merging similar clusterings.
#' @param ncores The number of parallel R instances that will run the code. If
#' the value is set to 1, the code will be run sequentially.
#' @param algorithm An index indicating which community detection algorithm will
#' be used: Louvain (1), Louvain refined (2), SLM (3) or Leiden (4). More details
#' can be found in the Seurat's `FindClusters` function.
#' @param ... additional arguments passed to the umap method.
#'
#'
#' @return A list having one field associated with a step value. Each step
#' contains a list with three fields:
#'
#' * ecc - the EC-Consistency of the partitions obtained on all repetitions
#' * embedding - one UMAP embedding generated on the feature subset
#' * most_frequent_partition - the most common partition obtained across repetitions
#' @md
#'
#' @note The algorithm assumes that the feature_set is already sorted when performing
#' the subsetting. For example, if the user wants to analyze highly variable feature set,
#' they should provide them sorted by their variability.
#'
#'
#' @export
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(100 * 10), runif(100 * 10, min = 3, max = 4)), nrow = 200, byrow = TRUE)
#' rownames(expr_matrix) <- as.character(1:200)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- get_feature_stability(
#'   data_matrix = t(expr_matrix),
#'   feature_set = colnames(expr_matrix),
#'   feature_type = "feature_name",
#'   steps = 5,
#'   npcs = 2,
#'   n_repetitions = 10,
#'   algorithm = 1,
#'   # the following parameters are used by the umap function and are not mandatory
#'   n_neighbors = 3,
#'   approx_pow = TRUE,
#'   n_epochs = 0,
#'   init = "random",
#'   min_dist = 0.3
#' )
#' plot_feature_stability_boxplot(feature_stability_result)
get_feature_stability <- function(data_matrix,
                                  feature_set,
                                  steps,
                                  feature_type,
                                  n_repetitions = 100,
                                  seed_sequence = NULL,
                                  graph_reduction_type = "PCA",
                                  npcs = 30,
                                  ecs_thresh = 1,
                                  ncores = 1,
                                  algorithm = 4,
                                  ...) {
  # check parameters
  if (!is.matrix(data_matrix) && !methods::is(data_matrix, "Matrix")) {
    stop("object parameter should be a matrix")
  }

  # transpose the matrix to have observations on rows and features on columns
  data_matrix <- t(data_matrix)

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

  if (!is.numeric(npcs) || length(npcs) > 1) {
    stop("npcs parameter should be numeric")
  }
  # convert npcs to integers
  npcs <- as.integer(npcs)

  if (!is.numeric(ecs_thresh)) {
    stop("ecs_thresh parameter should be numeric")
  }

  if (!is.numeric(algorithm) || length(algorithm) > 1 || !(algorithm %in% 1:4)) {
    stop("algorithm should be a number between 1 and 4")
  }

  if (!is.numeric(ncores) || length(ncores) > 1) {
    stop("ncores parameter should be numeric")
  }
  # convert number of cores to integers
  ncores <- as.integer(ncores)

  partitions_list <- list()

  object_name <- paste(feature_type, graph_reduction_type, ecs_thresh, sep = "_")
  steps_ecc_list <- list()
  steps_ecc_list[[object_name]] <- list()

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

  ncores <- min(ncores, length(seed_sequence), parallel::detectCores())

  # convert the negative steps to the size of the feature
  steps[steps <= 0 |
    steps > length(feature_set)] <- length(feature_set)
  steps <- sort(unique(steps))

  # store the additional arguments used by umap in a list
  suppl_args <- list(...)
  i <- 1
  while (i <= length(suppl_args)) {
    assign(names(suppl_args)[i], suppl_args[[i]])
    i <- i + 1
  }

  if (graph_reduction_type == "UMAP") {
    suppl_args[["n_threads"]] <- 1
    suppl_args[["n_sgd_threads"]] <- 1
  }
  umap_arg_names <- names(suppl_args)

  for (step in steps) {
    used_features <- feature_set[1:step]
    actual_npcs <- min(npcs, step %/% 2)

    partitions_list[[as.character(step)]] <- list()
    steps_ecc_list[[object_name]][[as.character(step)]] <- list()

    # keep only the features that we are using
    trimmed_matrix <- data_matrix[, used_features]

    # the variables needed in each PSOCK process
    needed_vars <- c("trimmed_matrix", "graph_reduction_type", "actual_npcs", "suppl_args", "algorithm")

    if (ncores > 1) {
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

    all_vars <- ls()
    seed <- 0

    # send the name of the umap arguments
    partitions_list[[as.character(step)]] <- foreach::foreach(
      seed = seed_sequence,
      .noexport = all_vars[!(all_vars %in% needed_vars)]
    ) %dopar% {
      # calculate the PCA embedding
      set.seed(seed)
      pca_embedding <- irlba::irlba(A = trimmed_matrix, nv = actual_npcs)
      pca_embedding <- pca_embedding$u %*% diag(pca_embedding$d)
      colnames(pca_embedding) <- paste0("PC_", 1:ncol(pca_embedding))
      rownames(pca_embedding) <- rownames(trimmed_matrix)

      if (graph_reduction_type == "UMAP") {
        # calculate the umap embedding
        set.seed(seed)
        umap_embedding <- do.call(uwot::umap, c(list(X = pca_embedding), suppl_args))
        colnames(umap_embedding) <- paste0("UMAP_", 1:ncol(umap_embedding))
        rownames(umap_embedding) <- rownames(trimmed_matrix)
      }

      # build the nn and snn graphs
      neigh_matrix <- Seurat::FindNeighbors(
        if (graph_reduction_type == "PCA") {
          pca_embedding
        } else {
          umap_embedding
        },
        nn.method = "rann",
        verbose = FALSE
      )$snn

      # apply graph clustering to the snn graph
      cluster_results <- Seurat::FindClusters(
        neigh_matrix,
        random.seed = seed,
        algorithm = algorithm,
        verbose = FALSE
      )
      cluster_results <- cluster_results[, names(cluster_results)[1]]

      # return the clustering
      list(
        mb = cluster_results,
        freq = 1,
        seed = seed
      )
    }

    # if a parallel backend was created, terminate the processes
    if (ncores > 1) {
      parallel::stopCluster(cl = my_cluster)
    }

    # generate an UMAP embedding that will be stored in the returned list
    set.seed(seed_sequence[1])
    pca_embedding <- irlba::irlba(A = trimmed_matrix, nv = actual_npcs)
    pca_embedding <- pca_embedding$u %*% diag(pca_embedding$d)
    colnames(pca_embedding) <- paste0("PC_", 1:ncol(pca_embedding))
    rownames(pca_embedding) <- rownames(trimmed_matrix)

    set.seed(seed_sequence[1])
    umap_embedding <- uwot::umap(X = pca_embedding, n_threads = 1, n_sgd_threads = 1, ...)
    colnames(umap_embedding) <- paste0("UMAP_", 1:ncol(umap_embedding))
    rownames(umap_embedding) <- rownames(data_matrix)

    # copy_p = partitions_list[[as.character(step)]]
    # merge the identical partitions into the same object
    partitions_list[[as.character(step)]] <- merge_partitions(partitions_list[[as.character(step)]],
      ecs_thresh = ecs_thresh,
      ncores = ncores,
      order = TRUE
    )

    # compute the EC-consistency of the partition list
    ec_consistency <- weighted_element_consistency(lapply(partitions_list[[as.character(step)]], function(x) {
      x$mb
    }),
    sapply(partitions_list[[as.character(step)]], function(x) {
      x$freq
    }),
    ncores = ncores
    )

    # update the returned object with a list containing the ecc, an UMAP embedding,
    # the most frequent partition and the number of different partitions
    steps_ecc_list[[object_name]][[as.character(step)]] <- list(
      ecc = ec_consistency,
      embedding = umap_embedding,
      most_frequent_partition = partitions_list[[as.character(step)]][[1]],
      # part_list = partitions_list[[as.character(step)]],
      # orig_pat_list = copy_p,
      n_different_partitions = length(partitions_list[[as.character(step)]])
    )
  }

  incremental_ecs_list <- list()
  incremental_ecs_list[[object_name]] <- list()

  if (length(steps) > 1) {
    for (i in 2:length(steps)) {
      incremental_ecs_list[[object_name]][[paste(steps[i - 1], steps[i], sep = "-")]] <- corrected_l1_mb(
        steps_ecc_list[[object_name]][[as.character(steps[i - 1])]]$most_frequent_partition$mb,
        steps_ecc_list[[object_name]][[as.character(steps[i])]]$most_frequent_partition$mb
      )
    }
  }

  list(
    steps_stability = steps_ecc_list,
    incremental_stability = incremental_ecs_list
  )
}

#' Feature Stability Boxplot
#'
#' @description Display EC consistency for each feature set and for each step.
#' Above each boxplot there is a number representing
#' the step (or the size of the subset)
#'
#' @param feature_object_list An object or a concatenation of objects returned by the
#' `get_feature_stability` method
#' @param text_size The size of the labels above boxplots
#' @param boxplot_width Used for adjusting the width of the boxplots; the value will
#' be passed to the `width` argument of the `ggplot2::geom_boxplot` method.
#' @param dodge_width Used for adjusting the horizontal position of the boxplot; the value
#' will be passed to the `width` argument of the `ggplot2::position_dodge` method.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(100 * 10), runif(100 * 10, min = 3, max = 4)), nrow = 200, byrow = TRUE)
#' rownames(expr_matrix) <- as.character(1:200)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- get_feature_stability(
#'   data_matrix = t(expr_matrix),
#'   feature_set = colnames(expr_matrix),
#'   feature_type = "feature_name",
#'   steps = 5,
#'   npcs = 2,
#'   n_repetitions = 10,
#'   algorithm = 1,
#'   # the following parameters are used by the umap function and are not mandatory
#'   n_neighbors = 3,
#'   approx_pow = TRUE,
#'   n_epochs = 0,
#'   init = "random",
#'   min_dist = 0.3
#' )
#' plot_feature_stability_boxplot(feature_stability_result)
plot_feature_stability_boxplot <- function(feature_object_list,
                                           text_size = 4,
                                           boxplot_width = 0.4,
                                           dodge_width = 0.7) {
  min_index <- -1 # indicates the number of steps that will be displayed on the plot
  feature_object_list <- feature_object_list$steps_stability

  # create a dataframe based on the object returned by `get_feature_stability`
  for (config_name in names(feature_object_list)) {
    melt_object <- reshape2::melt(lapply(feature_object_list[[config_name]], function(x) {
      x$ecc
    }))
    melt_object$L1 <- factor(melt_object$L1, levels = stringr::str_sort(unique(melt_object$L1), numeric = T))
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

  # generate the coordinates where the sizes of the steps will be displayed
  text_position <-
    stats::aggregate(value ~ L1 + feature_set, final_melt_df, max)
  text_position$value <- max(text_position$value) + 0.02

  # return the ggplot object
  ggplot2::ggplot(
    final_melt_df,
    ggplot2::aes(
      x = .data$L1,
      y = .data$value,
      fill = .data$feature_set
    )
  ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = dodge_width),
      width = boxplot_width
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
    ggplot2::ylab("EC consistency")
}

#' Feature Stability - Cluster Membership Facet Plot
#'
#' @description Display a facet of plots where each subpanel is associated with
#' a feature set and illustrates the distribution of the most frequent partition
#' over the UMAP embedding.
#'
#' @param feature_object_list An object or a concatenation of objects returned by the
#' `get_feature_stability` method
#' @param text_size The size of the cluster label
#' @param n_facet_cols The number of facet's columns.
#' @param point_size The size of the points displayed on the plot.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(100 * 10), runif(50 * 10, min = 5, max = 7)), nrow = 150, byrow = TRUE)
#' rownames(expr_matrix) <- as.character(1:150)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- get_feature_stability(
#'   data_matrix = t(expr_matrix),
#'   feature_set = colnames(expr_matrix),
#'   feature_type = "feature_name",
#'   steps = c(5, 10),
#'   npcs = 2,
#'   n_repetitions = 3,
#'   algorithm = 1,
#'   # the following parameters are used by the umap function and are not mandatory
#'   n_neighbors = 3,
#'   approx_pow = TRUE,
#'   n_epochs = 0,
#'   init = "random",
#'   min_dist = 0.3
#' )
#' plot_feature_stability_mb_facet(feature_stability_result)
plot_feature_stability_mb_facet <- function(feature_object_list,
                                            text_size = 5,
                                            n_facet_cols = 3,
                                            point_size = 0.3) {
  if (!is.numeric(text_size) || length(text_size) > 1) {
    stop("text_size parameter should be numeric")
  }

  first_temp <- T
  feature_object_list <- feature_object_list$steps_stability

  for (config_name in names(feature_object_list)) {
    for (steps in names(feature_object_list[[config_name]])) {
      temp_df <- data.frame(
        x = feature_object_list[[config_name]][[steps]]$embedding[, 1],
        y = feature_object_list[[config_name]][[steps]]$embedding[, 2],
        mb = feature_object_list[[config_name]][[steps]]$most_frequent_partition$mb
      )

      temp_df[["config_name"]] <- rep(config_name, nrow(temp_df))
      temp_df[["steps"]] <- rep(steps, nrow(temp_df))

      if (first_temp) {
        first_temp <- F
        final_df <- temp_df
      } else {
        final_df <- rbind(final_df, temp_df)
      }
    }
  }

  final_df$steps <- factor(final_df$steps, levels = stringr::str_sort(unique(final_df$steps), numeric = T))
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
#' @description Display a facet of plots where each subpanel is associated with a
#' feature set and illustrates the distribution of the EC consistency score
#' over the UMAP embedding.
#'
#' @param feature_object_list An object or a concatenation of objects returned by the
#' `get_feature_stability` method
#' @param n_facet_cols The number of facet's columns.
#' @param point_size The size of the points displayed on the plot.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(100 * 10), runif(50 * 10, min = 5, max = 7)), nrow = 150, byrow = TRUE)
#' rownames(expr_matrix) <- as.character(1:150)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- get_feature_stability(
#'   data_matrix = t(expr_matrix),
#'   feature_set = colnames(expr_matrix),
#'   feature_type = "feature_name",
#'   steps = c(5, 10),
#'   npcs = 2,
#'   n_repetitions = 3,
#'   algorithm = 1,
#'   # the following parameters are used by the umap function and are not mandatory
#'   n_neighbors = 3,
#'   approx_pow = TRUE,
#'   n_epochs = 0,
#'   init = "random",
#'   min_dist = 0.3
#' )
#' plot_feature_stability_ecs_facet(feature_stability_result)
plot_feature_stability_ecs_facet <- function(feature_object_list,
                                             n_facet_cols = 3,
                                             point_size = 0.3) {
  first_temp <- T
  feature_object_list <- feature_object_list$steps_stability

  for (config_name in names(feature_object_list)) {
    for (steps in names(feature_object_list[[config_name]])) {
      temp_df <- data.frame(
        x = feature_object_list[[config_name]][[steps]]$embedding[, 1],
        y = feature_object_list[[config_name]][[steps]]$embedding[, 2],
        ecc = feature_object_list[[config_name]][[steps]]$ecc
      )

      temp_df[["config_name"]] <- rep(config_name, nrow(temp_df))
      temp_df[["steps"]] <- rep(steps, nrow(temp_df))

      if (first_temp) {
        first_temp <- F
        final_df <- temp_df
      } else {
        final_df <- rbind(final_df, temp_df)
      }
    }
  }

  final_df$steps <- factor(final_df$steps, levels = stringr::str_sort(unique(final_df$steps), numeric = T))

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

#' Feature Stability Incremental Boxplot
#'
#' @description Perform an incremental ECS between two consecutive feature steps.
#'
#'
#' @param feature_object_list An object or a concatenation of objects returned by the
#' `get_feature_stability` method.
#' @param dodge_width Used for adjusting the horizontal position of the boxplot; the value
#' will be passed to the `width` argument of the `ggplot2::position_dodge` method.
#' @param text_size The size of the labels above boxplots.
#' @param boxplot_width Used for adjusting the width of the boxplots; the value will
#' be passed to the `width` argument of the `ggplot2::geom_boxplot` method.
#'
#' @return A ggplot2 object with ECS distribution will be displayed as a
#' boxplot. Above each boxplot there will be a pair of numbers representing the
#' two steps that are compared.
#' @export
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(25 * 10), runif(75 * 10, min = 5, max = 7)), nrow = 100, byrow = TRUE)
#' rownames(expr_matrix) <- as.character(1:100)
#' colnames(expr_matrix) <- paste("feature", 1:10)
#'
#' feature_stability_result <- get_feature_stability(
#'   data_matrix = t(expr_matrix),
#'   feature_set = colnames(expr_matrix),
#'   feature_type = "feature_name",
#'   steps = c(5, 10),
#'   npcs = 2,
#'   n_repetitions = 3,
#'   algorithm = 1,
#'   # the following parameters are used by the umap function and are not mandatory
#'   n_neighbors = 3,
#'   approx_pow = TRUE,
#'   n_epochs = 0,
#'   init = "random",
#'   min_dist = 0.3
#' )
#' plot_feature_stability_ecs_incremental(feature_stability_result)
plot_feature_stability_ecs_incremental <- function(feature_object_list,
                                                   dodge_width = 0.7,
                                                   text_size = 4,
                                                   boxplot_width = 0.4) {
  if (!is.numeric(dodge_width) || length(dodge_width) > 1) {
    stop("dodge_width parameter should be numeric")
  }

  min_index <- -1
  feature_object_list <- feature_object_list$incremental_stability
  first_df <- T

  for (config_name in names(feature_object_list)) {
    n.pairs <- length(feature_object_list[[config_name]])
    pairs.names <- names(feature_object_list[[config_name]])

    # treat the case with only one step
    for (i in 1:n.pairs) {
      first.n.steps <- strsplit(pairs.names[i], "-")[[1]][1]
      second.n.steps <- strsplit(pairs.names[i], "-")[[1]][2]

      temp_df <- data.frame(
        ecs = feature_object_list[[config_name]][[i]],
        index = i,
        feature_set = config_name
      )

      if (first_df) {
        first_df <- F

        steps_df <- data.frame(
          step = paste(
            first.n.steps,
            second.n.steps,
            sep = "-\n"
          ),
          index = i,
          feature_set = config_name,
          stringsAsFactors = FALSE
        )

        ecs_df <- temp_df
      } else {
        ecs_df <- rbind(ecs_df, temp_df)

        steps_df <- rbind(steps_df, c(
          paste(
            first.n.steps,
            second.n.steps,
            sep = "-\n"
          ),
          i,
          config_name
        ))
      }
    }

    if (min_index == -1 || n.pairs < min_index) {
      min_index <- n.pairs
    }
  }

  ecs_df <- ecs_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)
  steps_df <- steps_df %>% dplyr::filter(as.numeric(.data$index) <= min_index)

  ecs_df$feature_set <- factor(ecs_df$feature_set, levels = names(feature_object_list))
  steps_df$feature_set <- factor(steps_df$feature_set, levels = names(feature_object_list))

  ecs_df$index <- factor(ecs_df$index)

  text_position <-
    stats::aggregate(ecs ~ index + feature_set, ecs_df, max)
  text_position$ecs <- max(text_position$ecs) + 0.05

  ggplot2::ggplot(
    ecs_df,
    ggplot2::aes(
      x = .data$index,
      y = .data$ecs,
      fill = .data$feature_set
    )
  ) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = dodge_width), width = boxplot_width) +
    ggplot2::geom_text(
      data = text_position,
      ggplot2::aes(label = steps_df$step),
      position = ggplot2::position_dodge(width = dodge_width),
      size = text_size
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::xlab("# features") +
    ggplot2::ylab("EC similiarity")
}

######################## Graph construction ######################################

#### connected components ####

#' Relationship Between Nearest Neighbors and Connected Components
#'
#' @description One of the steps in the clustering pipeline is building a
#' k-nearest neighbor graph on a reduced-space embedding. This method assesses
#' the relationship between different number of nearest
#' neighbors and the connectivity of the graph. In the context of graph clustering,
#' the number of connected components can be used as a
#' lower bound for the number of clusters. The calculations are performed multiple
#' times by changing the seed at each repetition.
#'
#' @param object A data matrix. If the graph reduction type is PCA, the object
#' should be an expression matrix, with features on rows and observations on columns;
#' in the case of UMAP, the user could also provide a matrix associated to a PCA embedding.
#' See also the transpose argument.
#' @param n_neigh_sequence A sequence of the number of nearest neighbors.
#' @param config_name User specified string that uniquely describes the embedding characteristics.
#' @param n_repetitions The number of repetitions of applying the pipeline with different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence A custom seed sequence; if the value is NULL, the sequence will be built starting from 1 with a step of 100.
#' @param graph_reduction_type The graph reduction type, denoting if the graph should be built on either the PCA or the UMAP embedding.
#' @param transpose Logical: whether the input object will be transposed or not.
#' Set to FALSE if the input is an observations X features matrix, and set to TRUE
#' if the input is a features X observations matrix.
#' @param ncores The number of parallel R instances that will run the code. If the value is set to 1, the code will be run sequentially.
#' @param ... Additional arguments passed to the `irlba::irlba` or the `uwot::umap` method, depending on the value of graph_reduction_type.
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
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(50 * 10), runif(50 * 10, min = 1, max = 2)), nrow = 100, byrow = TRUE)
#' rownames(expr_matrix) <- as.character(1:100)
#'
#' # the graph reduction type is PCA, so we can provide the expression matrix as argument
#' nn_conn_comps_obj <- get_nn_conn_comps(
#'   object = expr_matrix,
#'   n_neigh_sequence = c(2, 3, 5),
#'   config_name = "example_config",
#'   n_repetitions = 10,
#'   graph_reduction_type = "PCA",
#'   transpose = FALSE,
#'   # the following parameter is used by the irlba function and is not mandatory
#'   nv = 3
#' )
#' plot_connected_comps_evolution(nn_conn_comps_obj)
get_nn_conn_comps <- function(object,
                              n_neigh_sequence,
                              config_name = "",
                              n_repetitions = 100,
                              seed_sequence = NULL,
                              graph_reduction_type = "UMAP",
                              transpose = (graph_reduction_type == "PCA"),
                              ncores = 1,
                              ...) {
  # check parameters
  if (!is.numeric(n_neigh_sequence)) {
    stop("n_neigh_sequence parameter should be numeric")
  }
  # convert number of neighbors to integers
  n_neigh_sequence <- as.integer(n_neigh_sequence)

  if (!is.numeric(ncores) || length(ncores) > 1) {
    stop("ncores parameter should be numeric")
  }
  # convert number of cores to integers
  ncores <- as.integer(ncores)

  if (!is.numeric(n_repetitions) || length(n_repetitions) > 1) {
    stop("n_repetitions parameter should be numeric")
  }
  # convert n_repetitions to integers
  n_repetitions <- as.integer(n_repetitions)

  if (!is.matrix(object) && !methods::is(object, "Matrix")) {
    stop("object parameter should be a matrix")
  }

  if (!(graph_reduction_type %in% c("PCA", "UMAP"))) {
    stop("graph_reduction_type parameter should take one of these values: 'PCA' or 'UMAP'")
  }

  if (!is.logical(transpose)) {
    stop("tranpose parameter should be logical")
  }

  if (!is.character(config_name)) {
    stop("config_name parameter should be a string")
  }

  # transpose the matrix if the parameter is set to TRUE
  if (transpose) {
    object <- t(object)
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

  ncores <- min(ncores, length(seed_sequence), parallel::detectCores())

  # store the additional arguments used by umap or irlba in a list
  suppl_args <- list(...)
  for (i in 1:length(suppl_args)) {
    assign(names(suppl_args)[i], suppl_args[[i]])
  }

  if (graph_reduction_type == "UMAP") {
    suppl_args[["n_threads"]] <- 1
    suppl_args[["n_sgd_threads"]] <- 1
  }
  suppl_arg_names <- names(suppl_args)

  nn_conn_comps_list <- list()

  if (ncores > 1) {
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

  all_vars <- ls()
  # the variables needed in each PSOCK process
  needed_vars <- c("object", "n_neigh_sequence", "graph_reduction_type", "suppl_args")

  seed <- NA

  # send the name of the dim reduction arguments
  nn_conn_comps_list_temp <- foreach::foreach(
    seed = seed_sequence,
    .noexport = all_vars[!(all_vars %in% needed_vars)]
  ) %dopar% { #
    # dim_red_args = list()
    #
    # for(suppl_arg_name in suppl_arg_names) {
    #   dim_red_args[[suppl_arg_name]] = get(suppl_arg_name)
    # }

    # perform the UMAP / PCA dimensionality reduction
    set.seed(seed)
    if (graph_reduction_type == "UMAP") {
      embedding <- do.call(uwot::umap, c(list(X = object), suppl_args))
      colnames(embedding) <- paste0("UMAP_", 1:ncol(embedding))
    } else {
      irlba_object <- do.call(irlba::irlba, c(list(A = object), suppl_args))
      embedding <- irlba_object$u %*% diag(irlba_object$d)
      colnames(embedding) <- paste0("PC_", 1:ncol(embedding))
    }
    rownames(embedding) <- rownames(object)

    # for each neighbour, return the number of connected components
    # of the generated graph
    sapply(n_neigh_sequence, function(n_neigh) {
      g <- igraph::graph_from_adjacency_matrix(
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
  if (ncores > 1) {
    parallel::stopCluster(cl = my_cluster)
  }

  # store the results obtained for each number of neighbours in different lists
  for (i in 1:length(n_neigh_sequence)) {
    n_neigh <- n_neigh_sequence[i]
    nn_conn_comps_list[[as.character(n_neigh)]] <- sapply(nn_conn_comps_list_temp, function(x) {
      x[i]
    })
  }

  nn_conn_comps_list <- list(nn_conn_comps_list)
  names(nn_conn_comps_list) <- c(paste(config_name, graph_reduction_type, sep = "_"))

  nn_conn_comps_list
}

#' Relationship Between Number of Nearest Neighbors and Graph Connectivity
#'
#' @description Display the distribution of the number connected components
#' obtained for each number of neighbors across random seeds.
#'
#' @param nn_conn_comps_object An object or a concatenation of objects returned by the
#' `get_nn_conn_comps` method.
#'
#'
#' @return A ggplot2 object with boxplots for the connected component distributions.
#' @export
#'
#' @note The number of connected components is displayed on a logarithmic scale.
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(50 * 10), runif(50 * 10, min = 1, max = 2)), nrow = 100, byrow = TRUE)
#' rownames(expr_matrix) <- as.character(1:100)
#'
#' # the graph reduction type is PCA, so we can provide the expression matrix as argument
#' nn_conn_comps_obj <- get_nn_conn_comps(
#'   object = expr_matrix,
#'   n_neigh_sequence = c(2, 3, 5),
#'   config_name = "example_config",
#'   n_repetitions = 10,
#'   graph_reduction_type = "PCA",
#'   transpose = FALSE,
#'   # the following parameter is used by the irlba function and is not mandatory
#'   nv = 3
#' )
#' plot_connected_comps_evolution(nn_conn_comps_obj)
plot_connected_comps_evolution <- function(nn_conn_comps_object) {
  final_comps_df <- reshape2::melt(nn_conn_comps_object)
  colnames(final_comps_df) <- c("n_comps", "n_neigh", "config_name")

  final_comps_df$config_name <- factor(final_comps_df$config_name)
  final_comps_df$n_neigh <- factor(final_comps_df$n_neigh)
  final_comps_df$n_neigh <- factor(final_comps_df$n_neigh, levels(final_comps_df$n_neigh)[stringr::str_order(levels(final_comps_df$n_neigh), numeric = T)])

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
      x = "# of nearest neighbors",
      y = "# of connected components",
      fill = "configuration"
    ) +
    ggplot2::ggtitle("Distribution of the number of connected components")
}

#### number of neigh, graph type importance ####

#' Assess Graph Building Parameters
#'
#' @description Evaluates clustering stability when changing the values of different
#' parameters involved in the graph building step,
#' namely the base embedding, the graph type and the number of neighbours.
#'
#' @param object The data matrix. If the graph reduction type is PCA, the object
#' should be an expression matrix, with features on rows and observations on columns;
#' in the case of UMAP, the user could also provide a matrix associated to a PCA embedding.
#' See also the transpose argument.
#' @param n_neigh_sequence A sequence of the number of nearest neighbours.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence A custom seed sequence; if the value is NULL,
#' the sequence will be built starting from 1 with a step of 100.
#' @param graph_reduction_type The graph reduction type, denoting if the graph
#' should be built on either the PCA or the UMAP embedding.
#' @param ecs_thresh The ECS threshold used for merging similar clusterings.
#' @param ncores The number of parallel R instances that will run the code.
#' If the value is set to 1, the code will be run sequentially.
#' @param transpose Logical: whether the input object will be transposed or not.
#' Set to FALSE if the input is an observations X features matrix, and set to TRUE
#' if the input is a features X observations matrix.
#' @param graph_type Argument indicating whether the graph should be
#' unweighted (0), weighted (1) or both (2).
#' @param algorithm An index indicating which community detection algorithm will
#' be used: Louvain (1), Louvain refined (2), SLM (3) or Leiden (4). More details
#' can be found in the Seurat's `FindClusters` function.
#' @param ... Additional arguments passed to the `irlba::irlba` or the `uwot::umap`
#' method, depending on the value of graph_reduction_type.
#'
#'
#'
#' @return A list having three fields:
#'
#'
#' * n_neigh_k_corresp - list containing the number of the clusters obtained by running
#' the pipeline multiple times with different seed, number of neighbors and graph type (weighted vs unweigted)
#' * n_neigh_ec_consistency - list containing the EC consistency of the partitions obtained
#' at multiple runs when changing the number of neighbors or the graph type
#' * n_different_partitions - the number of different partitions obtained by each
#' number of neighbors
#'
#' @md
#' @export
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(100 * 10), runif(100 * 10, min = 5, max = 6)), nrow = 200)
#' rownames(expr_matrix) <- as.character(1:200)
#'
#' nn_importance_obj <- get_nn_importance(
#'   object = expr_matrix,
#'   n_neigh_sequence = c(10, 15, 20),
#'   n_repetitions = 10,
#'   graph_reduction_type = "PCA",
#'   algorithm = 1,
#'   transpose = FALSE, # the matrix is already observations x features, so we won't transpose it
#'   # the following parameter is used by the irlba function and is not mandatory
#'   nv = 2
#' )
#' plot_n_neigh_ecs(nn_importance_obj)
get_nn_importance <- function(object,
                              n_neigh_sequence,
                              n_repetitions = 100,
                              seed_sequence = NULL,
                              graph_reduction_type = "PCA",
                              ecs_thresh = 1,
                              ncores = 1,
                              transpose = (graph_reduction_type == "PCA"),
                              graph_type = 2,
                              algorithm = 4,
                              config_name = "",
                              ...) {
  # check parameters
  if (!is.logical(transpose)) {
    stop("tranpose parameter should be logical")
  }

  if (!is.numeric(n_neigh_sequence)) {
    stop("n_neigh_sequence parameter should be numeric")
  }
  # convert number of neighbors to integers
  n_neigh_sequence <- as.integer(n_neigh_sequence)

  if (!is.numeric(ncores) || length(ncores) > 1) {
    stop("ncores parameter should be numeric")
  }
  # convert number of cores to integers
  ncores <- as.integer(ncores)

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

  if (!is.matrix(object) && !methods::is(object, "Matrix")) {
    stop("object parameter should be a matrix")
  }

  if (!is.numeric(algorithm) || length(algorithm) > 1 || !(algorithm %in% 1:4)) {
    stop("algorithm should be a number between 1 and 4")
  }

  if (!(graph_reduction_type %in% c("PCA", "UMAP"))) {
    stop("graph_reduction_type parameter should take one of these values: 'PCA' or 'UMAP'")
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

  partitions_list <- list()

  if (graph_type != 1) {
    partitions_list[[paste(graph_reduction_type, "nn", ecs_thresh, sep = "_")]] <- list()
  }

  if (graph_type != 0) {
    partitions_list[[paste(graph_reduction_type, "snn", ecs_thresh, sep = "_")]] <- list()
  }

  # transpose the matrix if the parameter is set to TRUE
  if (transpose) {
    object <- t(object)
  }

  ncores <- min(ncores, length(seed_sequence), parallel::detectCores())
  print(nrow(object))

  # additional arguments that will be used by umap or irlba
  suppl_args <- list(...)
  i <- 1
  while (i < length(suppl_args)) {
    assign(names(suppl_args)[i], suppl_args[[i]])
    i <- i + 1
  }

  if (graph_reduction_type == "UMAP") {
    suppl_args[["n_threads"]] <- 1
    suppl_args[["n_sgd_threads"]] <- 1
  }

  suppl_arg_names <- names(suppl_args)

  for (n_neigh in n_neigh_sequence) {
    partitions_list[[paste(graph_reduction_type, "snn", ecs_thresh, sep = "_")]][[as.character(n_neigh)]] <- list()

    if (ncores > 1) {
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
    needed_vars <- c("object", "n_neigh", "graph_reduction_type", "graph_type", "algorithm", "suppl_args")
    all_vars <- ls()
    seed <- NA

    partitions_list_temp <- foreach::foreach(
      seed = seed_sequence,
      .noexport = all_vars[!(all_vars %in% needed_vars)]
    ) %dopar% {
      # perform the dimensionality reduction
      set.seed(seed)
      if (graph_reduction_type == "UMAP") {
        embedding <- do.call(uwot::umap, c(list(X = object), suppl_args))
        colnames(embedding) <- paste0("UMAP_", 1:ncol(embedding))
      } else {
        irlba_object <- do.call(irlba::irlba, c(list(A = object), suppl_args))
        embedding <- irlba_object$u %*% diag(irlba_object$d)
        colnames(embedding) <- paste0("PC_", 1:ncol(embedding))
      }
      rownames(embedding) <- rownames(object)

      # build the nn and snn graphs
      neigh_matrix <- Seurat::FindNeighbors(
        embedding,
        k.param = n_neigh,
        nn.method = "rann",
        verbose = F
      )

      # apply the clustering method on the graph specified by the variable `graph_type`
      if (graph_type != 1) {
        cluster_results_nn <- Seurat::FindClusters(
          neigh_matrix$nn,
          random.seed = seed,
          algorithm = algorithm,
          verbose = F
        )
        cluster_results_nn <- list(
          mb = cluster_results_nn[, names(cluster_results_nn)[1]],
          freq = 1,
          seed = seed
        )

        if (graph_type == 0) {
          return(cluster_results_nn)
        }
      }

      cluster_results_snn <- Seurat::FindClusters(
        neigh_matrix$snn,
        random.seed = seed,
        algorithm = algorithm,
        verbose = F
      )
      cluster_results_snn <- list(
        mb = cluster_results_snn[, names(cluster_results_snn)[1]],
        freq = 1,
        seed = seed
      )

      if (graph_type == 1) {
        return(cluster_results_snn)
      }

      return(list(cluster_results_nn, cluster_results_snn))
    }

    # if a parallel backend was created, terminate the processes
    if (ncores > 1) {
      parallel::stopCluster(cl = my_cluster)
    }

    print("merge partitions")
    # merge the partitions that are considered similar by a given ecs threshold
    for (i in 1:length(partitions_list)) {
      partitions_list[[i]][[as.character((n_neigh))]] <- merge_partitions(lapply(partitions_list_temp, function(x) {
        x[[i]]
      }),
      ecs_thresh = ecs_thresh,
      ncores = ncores
      )
    }
  }

  names(partitions_list) <- paste(config_name, names(partitions_list), sep = "_")

  # create an object showing the number of clusters obtained for each number of
  # neighbours
  nn_object_n_clusters <- list()
  for (config_name in names(partitions_list)) {
    nn_object_n_clusters[[config_name]] <- list()
    for (n_neigh in names(partitions_list[[config_name]])) {
      nn_object_n_clusters[[config_name]][[n_neigh]] <- unlist(lapply(partitions_list[[config_name]][[n_neigh]], function(x) {
        rep(length(unique(x$mb)), x$freq)
      }))
    }
  }

  print("calculating weighted ecc")

  # create an object that contains the ecc of the partitions obtained for each
  # number of neighbours
  nn_ecs_object <- lapply(partitions_list, function(config) {
    lapply(config, function(n_neigh) {
      weighted_element_consistency(lapply(n_neigh, function(x) {
        x$mb
      }),
      sapply(n_neigh, function(x) {
        x$freq
      }),
      ncores = ncores
      )
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

# 2. neighbors <-> number of clusters association

#' Relationship Between Number of Nearest Neighbors and Number of Clusters
#'
#' @description Display the distribution of the
#' number of clusters obtained for each number of neighbors across random seeds.
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
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(runif(100 * 10), nrow = 100)
#' rownames(expr_matrix) <- as.character(1:100)
#'
#' nn_importance_obj <- get_nn_importance(
#'   object = expr_matrix,
#'   n_neigh_sequence = c(2, 5),
#'   n_repetitions = 5,
#'   graph_reduction_type = "PCA",
#'   algorithm = 1,
#'   transpose = FALSE, # the matrix is already observations x features, so we won't transpose it
#'   # the following parameter is used by the irlba function and is not mandatory
#'   nv = 2
#' )
#' plot_n_neigh_k_correspondence(nn_importance_obj)
plot_n_neigh_k_correspondence <- function(nn_object_n_clusters) {
  melted_obj <- reshape2::melt(nn_object_n_clusters$n_neigh_k_corresp)
  colnames(melted_obj) <- c("k", "n_neigh", "config_name")

  melted_obj$n_neigh <- factor(melted_obj$n_neigh)
  melted_obj$n_neigh <- factor(melted_obj$n_neigh, levels(melted_obj$n_neigh)[stringr::str_order(levels(melted_obj$n_neigh), numeric = T)])

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
      x = "# of nearest neighbors",
      fill = "configuration"
    ) +
    ggplot2::ggtitle("Distribution of k across different seeds")
}

# 3. ECS distribution across different seeds for different number of neighbors

#' Graph construction parameters - ECC facet
#'
#' @description Display, for all configurations consisting in different number
#' of neighbors, graph types and base embeddings, the EC Consistency of the partitions
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
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(100 * 10), runif(100 * 10, min = 5, max = 6)), nrow = 200)
#' rownames(expr_matrix) <- as.character(1:200)
#'
#' nn_importance_obj <- get_nn_importance(
#'   object = expr_matrix,
#'   n_neigh_sequence = c(10, 15, 20),
#'   n_repetitions = 10,
#'   graph_reduction_type = "PCA",
#'   algorithm = 1,
#'   transpose = FALSE, # the matrix is already observations x features, so we won't transpose it
#'   # the following parameter is used by the irlba function and is not mandatory
#'   nv = 2
#' )
#' plot_n_neigh_ecs(nn_importance_obj)
plot_n_neigh_ecs <- function(nn_ecs_object,
                             boxplot_width = 0.5) {
  melted_obj <- reshape2::melt(nn_ecs_object$n_neigh_ec_consistency)
  colnames(melted_obj) <- c("ECC", "n_neigh", "config_name")

  melted_obj$n_neigh <- factor(melted_obj$n_neigh)
  melted_obj$n_neigh <- factor(melted_obj$n_neigh, levels(melted_obj$n_neigh)[stringr::str_order(levels(melted_obj$n_neigh), numeric = T)])

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
      x = "# of nearest neighbors",
      y = "EC consistency",
      fill = "configuration"
    ) +
    ggplot2::ggtitle("Distribution of ECC across different seeds for different # neighbors")
}

############################## Clustering ########################################

#### clustering methods ####

#' Graph Clustering Method Stability
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
#'
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(runif(100 * 10), nrow = 100)
#' rownames(expr_matrix) <- as.character(1:100)
#'
#' adj_matrix <- Seurat::FindNeighbors(expr_matrix,
#'   k.param = 10,
#'   nn.method = "rann",
#'   verbose = FALSE,
#'   compute.SNN = FALSE
#' )$nn
#' clust_diff_obj <- get_clustering_difference(
#'   graph_adjacency_matrix = adj_matrix,
#'   resolution = c(0.5, 1),
#'   n_repetitions = 10,
#'   algorithm = 1:2,
#'   verbose = FALSE
#' )
#' plot_clustering_difference_boxplot(clust_diff_obj)
get_clustering_difference <- function(graph_adjacency_matrix,
                                      resolution,
                                      n_repetitions = 100,
                                      seed_sequence = NULL,
                                      ecs_thresh = 1,
                                      ncores = 1,
                                      algorithm = 1:4,
                                      verbose = TRUE) {
  # check the parameters
  if (!is.numeric(resolution)) {
    stop("resolution parameter should be numeric")
  }

  if (!is.numeric(ecs_thresh) || length(ecs_thresh) > 1) {
    stop("ecs_thresh parameter should be numeric")
  }

  if (!is.numeric(ncores) || length(ncores) > 1) {
    stop("ncores parameter should be numeric")
  }
  # convert ncores to an integer
  ncores <- as.integer(ncores)

  if (!is.numeric(n_repetitions) || length(n_repetitions) > 1) {
    stop("n_repetitions parameter should be numeric")
  }
  # convert n_repetitions to an integer
  n_repetitions <- as.integer(n_repetitions)

  if (!is.numeric(algorithm) || !(all(algorithm %in% 1:4))) {
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

  for (i in algorithm) {
    # if verbose is set to TRUE, create a progress bar
    if (verbose) {
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
    result_object[[algorithm_names[i]]] <- sapply(resolution, function(res) {
      if (verbose) {
        pb$tick(tokens = list(what = res))
      }
      res_list <- list()
      res_list[[as.character(res)]] <- get_resolution_partitions(
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
  filtered_result <- lapply(result_object, function(config_object) {
    lapply(config_object, function(res) {
      corresp_k <- names(which.max(unlist(lapply(res, function(res_cluster) {
        sum(sapply(res_cluster$partitions, function(partition) {
          partition$freq
        }))
      }))))[1]

      # list(ec_consistency = weighted_element_consistency(lapply(res[[corresp_k]]$partitions, function(x) {
      #   x$mb
      # }),
      # sapply(res[[corresp_k]]$partitions, function(x) {
      #   x$freq
      # }),
      # ncores = ncores),
      list(ec_consistency = res[[corresp_k]]$ecc, k = corresp_k)
    })
  })

  # calculate ecc for all the partitions created with a resolution value
  all_result <- lapply(result_object, function(config_object) {
    lapply(config_object, function(res) {
      k_names <- names(res)
      partition_list <- list()
      for (k in k_names) {
        partition_list <- c(partition_list, res[[k]]$partitions)
      }

      weighted_element_consistency(lapply(partition_list, function(x) {
        x$mb
      }),
      sapply(partition_list, function(x) {
        x$freq
      }),
      ncores = ncores
      )
    })
  })

  list(
    filtered = filtered_result,
    all = all_result
  )
}

#' Clustering Method Stability Boxplot
#'
#' @description Display EC consistency across clustering method and resolution
#' values. The `filtered` field of the object returned by the
#' `get_clustering_difference_object` method is used.
#' Above each boxplot, the number of clusters is displayed.
#'
#' @param clustering_difference_object An object returned by the
#' `get_clustering_difference_object` method.
#' @param text_size The size of the labels above boxplots.
#' @param boxplot_width Used for adjusting the width of the boxplots; the value will
#' be passed to the `width` argument of the `ggplot2::geom_boxplot` method.
#' @param dodge_width Used for adjusting the horizontal position of the boxplot; the value
#' will be passed to the `width` argument of the `ggplot2::position_dodge` method.
#'
#' @return A ggplot2 object with the EC consistency distributions. Higher
#' consistency indicates a more stable clustering.
#' @export
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(runif(500 * 10), nrow = 500)
#' rownames(expr_matrix) <- as.character(1:500)
#'
#' adj_matrix <- Seurat::FindNeighbors(expr_matrix,
#'   k.param = 10,
#'   nn.method = "rann",
#'   verbose = FALSE,
#'   compute.SNN = FALSE
#' )$nn
#' clust_diff_obj <- get_clustering_difference(
#'   graph_adjacency_matrix = adj_matrix,
#'   resolution = c(0.5, 1),
#'   n_repetitions = 10,
#'   algorithm = 1:2,
#'   verbose = FALSE
#' )
#' plot_clustering_difference_boxplot(clust_diff_obj)
plot_clustering_difference_per_resolution_boxplot <- function(clustering_difference_object,
                                                              text_size = 3,
                                                              dodge_width = 0.9,
                                                              boxplot_width = 0.4) {
  # get the labels representing the number of clusters
  k_labels <- reshape2::melt(lapply(clustering_difference_object$filtered, function(config_object) {
    lapply(config_object, function(res) {
      res$k
    })
  })) %>% dplyr::arrange(.data$L1)

  melted_obj <- reshape2::melt(lapply(clustering_difference_object$filtered, function(config_object) {
    lapply(config_object, function(res) {
      res$ec_consistency
    })
  }))

  colnames(melted_obj) <- c("EC_consistency", "resolution", "clustering_method")

  # calculate the coordinates where the k labels will be displayed
  text_position <-
    stats::aggregate(
      EC_consistency ~ resolution + clustering_method,
      melted_obj,
      max
    )
  text_position$EC_consistency <- max(text_position$EC_consistency) + 0.02

  ggplot2::ggplot(
    melted_obj,
    ggplot2::aes(
      x = .data$resolution,
      y = .data$EC_consistency,
      fill = .data$clustering_method
    )
  ) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = dodge_width), width = boxplot_width) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      y = "EC Consistency",
      fill = "clustering method"
    ) +
    ggplot2::geom_text(
      data = text_position,
      ggplot2::aes(label = k_labels$value),
      position = ggplot2::position_dodge(width = dodge_width),
      size = text_size
    )
}

plot_clustering_difference_boxplot <- function(clustering_object) {
  per_method_ecc <- reshape2::melt(lapply(clustering_object$all, function(clust_method) {
    c(sapply(clust_method, c))
  }))

  ggplot2::ggplot(
    per_method_ecc,
    ggplot2::aes(
      x = .data$L1,
      y = .data$value
    )
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::labs( # position = ggplot2::position_dodge(width = dodge_width), width = boxplot_width) +
      y = "EC Consistency",
      x = ""
    ) +
    ggplot2::theme_classic()
}

#' Clustering Method Stability Facet Plot
#'
#' @description Display the distribution of the EC consistency for each
#' clustering method and each resolution value on a given embedding The `all`
#' field of the object returned by the `get_clustering_difference_object` method is used.
#'
#' @param clustering_difference_object An object returned by the
#' `get_clustering_difference_object` method.
#' @param embedding An embedding (only the first two dimensions will be used for
#' visualization).
#' @param low_limit The lowest value of ECC that will be displayed on the embedding.
#' @param high_limit The highest value of ECC that will be displayed on the embedding.
#' @param grid Boolean value indicating whether the facet should be a grid (where each
#' row is associated with a resolution value and each column with a clustering method) or
#' a wrap.
#'
#' @return A ggplot2 object.
#' @export
#'
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(c(runif(250 * 10), runif(250 * 10, min = 5, max = 7)), nrow = 500)
#' rownames(expr_matrix) <- as.character(1:500)
#'
#' pca_embedding <- irlba::irlba(expr_matrix, nv = 2)
#' pca_embedding <- pca_embedding$u %*% diag(pca_embedding$d)
#' rownames(pca_embedding) <- as.character(1:500)
#'
#' adj_matrix <- Seurat::FindNeighbors(pca_embedding,
#'   k.param = 10,
#'   nn.method = "rann",
#'   verbose = FALSE,
#'   compute.SNN = FALSE
#' )$nn
#' clust_diff_obj <- get_clustering_difference(
#'   graph_adjacency_matrix = adj_matrix,
#'   resolution = c(0.5, 1),
#'   n_repetitions = 10,
#'   algorithm = 1:2,
#'   verbose = FALSE
#' )
#' plot_clustering_difference_facet(clust_diff_obj, pca_embedding)
plot_clustering_difference_facet <- function(clustering_difference_object,
                                             embedding,
                                             low_limit = 0,
                                             high_limit = 1,
                                             grid = TRUE) {
  # check parameters
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
  if (length(clustering_difference_object$all[[1]][[1]]) != npoints) {
    stop("The provided embedding and the consistency arrays must have the same number of elements!")
  }

  if (ncol(embedding) < 2) {
    stop("The embedding should have at least two dimensions!")
  }

  melt_obj <- reshape2::melt(clustering_difference_object$all)

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

#### resolution ####

# using a graph-based clustering algorithm, determine the partitions
# obtained on a graph using different resolution values and different seeds
get_resolution_partitions <- function(clustered_object,
                                      resolution,
                                      clustering_function,
                                      seed_sequence,
                                      ncores = 1,
                                      ecs_thresh = 1,
                                      ...) {
  different_partitions <- list()
  ncores <- min(ncores, length(seed_sequence), parallel::detectCores())

  # additional arguments used by the clustering method
  suppl_args <- list(...)
  i <- 1
  while (i <= length(suppl_args)) {
    assign(names(suppl_args)[i], suppl_args[[i]])
    i <- i + 1
  }

  suppl_arg_names <- names(suppl_args)

  if (ncores > 1) {
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
  needed_vars <- c("resolution", "clustering_function", "clustered_object", "suppl_args")
  all_vars <- ls()

  seed <- 0

  different_partitions_temp <- foreach::foreach(
    seed = seed_sequence,
    .noexport = all_vars[!(all_vars %in% needed_vars)]
  ) %dopar% { #
    seed <- get("seed")

    # apply the clustering, which should return a membership vector
    do.call(clustering_function, c(
      list(
        object = clustered_object,
        resolution = resolution,
        seed = seed
      ),
      suppl_args
    ))
  }

  # if a parallel backend was created, terminate the processes
  if (ncores > 1) {
    parallel::stopCluster(cl = my_cluster)
  }

  # group the partitions by the number of clusters
  for (i in 1:length(seed_sequence)) {
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

  # merge the partitions using the ecs threshold
  for (k in names(different_partitions)) {
    different_partitions[[k]] <- list(partitions = merge_partitions(different_partitions[[k]],
      ecs_thresh = ecs_thresh,
      ncores = ncores,
      order = TRUE
    ))

    # compute the EC-consistency of the partition list
    ec_consistency <- weighted_element_consistency(lapply(different_partitions[[k]]$partitions, function(x) {
      x$mb
    }),
    sapply(different_partitions[[k]]$partitions, function(x) {
      x$freq
    }),
    ncores = ncores
    )

    different_partitions[[k]][["ecc"]] <- ec_consistency
  }

  different_partitions
}

#' Evaluate Stability Across Resolution, Number of Neighbors, and Graph Type
#'
#' @description Perform a grid search over the resolution, number of neighbors
#' and graph type.
#'
#' @param embedding The base embedding for the graph construction.
#' @param resolution A sequence of resolution values.
#' @param n_neigh A value or a sequence of number of neighbors used for graph construction.
#' @param n_repetitions The number of repetitions of applying the pipeline with
#' different seeds; ignored if seed_sequence is provided by the user.
#' @param seed_sequence A custom seed sequence; if the value is NULL, the sequence
#' will be built starting from 1 with a step of 100.
#' @param clustering_method An index or a list of indexes indicating which community detection
#' algorithm will be used: Louvain (1), Louvain refined (2), SLM (3) or Leiden (4).
#' More details can be found in the Seurat's `FindClusters` function.
#' @param graph_type Argument indicating whether the graph should be
#' unweighted (0), weighted (1) or both (2).
#' @param object_name User specified string that uniquely describes the
#' embedding characteristics.
#' @param ecs_thresh The ECS threshold used for merging similar clusterings.
#' @param ncores The number of parallel R instances that will run the code.
#' If the value is set to 1, the code will be run sequentially.
#'
#'
#'
#' @return A list having two fields:
#'
#' * split_by_resolution: A five-level list. The hierarchy is as follows:
#'
#'     * the configuration name: concatenation between the object name provided by
#' the user, the number of neighbors, the graph type and the clustering method
#'     * the resolution value \eqn{\gamma}
#'     * the number of clusters *k* that can be obtained using the specified resolution
#'     * a `partitions` field containing the partitions obtained with resolution \eqn{\gamma} and have *k* clusters
#'     along with a `ecc` field that contains the Element-Centric Consistency of the partition list
#'     * the structure of a partition, which consists in having a `mb` field with
#' the flat membership vector, `freq` denoting its frequency and `seed`, that is
#' the seed used to obtain this partition in this configuration.
#'
#' * split_by_k: has a similar structure, but the resolution level is removed.
#' The partitions obtained in a configuration with the same number of clusters
#' will be merged into the same list.
#'
#' @md
#' @export
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(runif(500 * 10), nrow = 500)
#'
#' # get the PCA embedding of the data
#' pca_embedding <- irlba::irlba(expr_matrix, nv = 2)
#' pca_embedding <- pca_embedding$u %*% diag(pca_embedding$d)
#' rownames(pca_embedding) <- as.character(1:500)
#'
#' # run the function on the pca embedding
#' resolution_result <- get_resolution_importance(
#'   embedding = pca_embedding,
#'   resolution = c(0.8, 1),
#'   n_neigh = c(5, 7),
#'   n_repetitions = 5,
#'   clustering_method = 1,
#'   graph_type = 2,
#'   object_name = "name_example"
#' )
#'
#' plot_k_resolution_corresp(resolution_result)
get_resolution_importance <- function(embedding,
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

  different_partitions <- list()

  # check the parameters
  if (!is.numeric(resolution)) {
    stop("resolution parameter should be numeric")
  }

  if (!is.numeric(n_neigh)) {
    stop("n_neigh parameter should be numeric")
  }

  # convert number of neighbors to integers
  n_neigh <- sapply(n_neigh, as.integer)

  if (!is.numeric(clustering_method) ||
    !all(clustering_method %in% 1:4)) {
    stop("clustering_method should be a number between 1 and 4")
  }

  if (!is.numeric(graph_type) || length(graph_type) > 1 || !(graph_type %in% 0:2)) {
    stop("graph_type should be a number between 0 and 2")
  }

  if (!is.numeric(ecs_thresh) || length(ecs_thresh) > 1) {
    stop("ecs_thresh parameter should be numeric")
  }

  if (!is.numeric(ncores) || length(ncores) > 1) {
    stop("ncores parameter should be numeric")
  }
  # convert ncores to an integer
  ncores <- as.integer(ncores)

  if (!is.numeric(n_repetitions) || length(n_repetitions) > 1) {
    stop("n_repetitions parameter should be numeric")
  }
  # convert n_repetitions to an integer
  n_repetitions <- as.integer(n_repetitions)

  if (!is.null(object_name) && !is.character(object_name)) {
    stop("object_name parameter should be a string")
  }

  if (length(clustering_method) > 1) {
    # generate a list of partitions for each clustering method chosen by the user
    different_partitions <- get_resolution_importance(
      embedding = embedding,
      resolution = resolution,
      n_neigh = n_neigh,
      n_repetitions = n_repetitions,
      seed_sequence = seed_sequence,
      clustering_method = clustering_method[1],
      graph_type = graph_type,
      object_name = object_name,
      ncores = ncores,
      ecs_thresh = ecs_thresh
    )

    i <- 2
    while (i <= length(clustering_method)) {
      temp_res_imp <- get_resolution_importance(
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

      different_partitions$split_by_resolution <- c(
        different_partitions$split_by_resolution,
        temp_res_imp$split_by_resolution
      )
      different_partitions$split_by_k <- c(
        different_partitions$split_by_k,
        temp_res_imp$split_by_k
      )
      i <- i + 1
    }

    return(different_partitions)
  }

  if (length(n_neigh) > 1) {
    # generate a list of partitions for each number of neighbours chosen by the user
    different_partitions <- get_resolution_importance(
      embedding = embedding,
      resolution = resolution,
      n_neigh = n_neigh[1],
      n_repetitions = n_repetitions,
      seed_sequence = seed_sequence,
      clustering_method = clustering_method,
      graph_type = graph_type,
      object_name = object_name,
      ncores = ncores,
      ecs_thresh = ecs_thresh
    )

    i <- 2
    while (i <= length(n_neigh)) {
      temp_res_imp <- get_resolution_importance(
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
      different_partitions$split_by_resolution <- c(
        different_partitions$split_by_resolution,
        temp_res_imp$split_by_resolution
      )
      different_partitions$split_by_k <- c(
        different_partitions$split_by_k,
        temp_res_imp$split_by_k
      )
      i <- i + 1
    }

    return(different_partitions)
  }

  if (graph_type == 2) {
    # generate a list of partitions for both types of graphs: nn and snn
    different_partitions <- mapply(c,
      get_resolution_importance(
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
      get_resolution_importance(
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
      ),
      SIMPLIFY = FALSE
    )

    return(different_partitions)
  }

  graph_type_name <- ifelse(graph_type == 0, "NN", "SNN")
  algorithm_name <- switch(clustering_method,
    "Louvain",
    "Louvain.refined",
    "SLM",
    "Leiden"
  )

  details <- paste(n_neigh, graph_type_name, algorithm_name, sep = "_")

  # generate a list of partitions for each resolution value
  for (i in 1:length(resolution)) {
    res <- resolution[i]

    if (i == 1) {
      different_partitions <- one_resolution_importance(
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

      used_name <- names(different_partitions)[1]
    } else {
      different_partitions[[used_name]] <- c(
        different_partitions[[used_name]],
        one_resolution_importance(
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

  return_object <- list(
    split_by_resolution = different_partitions,
    split_by_k = merge_resolutions(different_partitions[[1]],
      ncores = ncores,
      object_name = used_name
    )
  )

  return(return_object)
}

# call the get_resolution_partitions using a single resolution value
# and the clustering methods provided by Seurat
one_resolution_importance <- function(embedding,
                                      resolution,
                                      n_neigh,
                                      n_repetitions = 100,
                                      seed_sequence = NULL,
                                      clustering_method = 4,
                                      graph_type = 0,
                                      object_name = NULL,
                                      ecs_thresh = 1,
                                      ncores = 1) {
  graph_type_name <- ifelse(graph_type == 0, "NN", "SNN")
  algorithm_name <- switch(clustering_method,
    "Louvain",
    "Louvain.refined",
    "SLM",
    "Leiden"
  )

  details <- paste(n_neigh, graph_type_name, algorithm_name, sep = "_")

  if (is.null(object_name)) {
    object_name <- details
  } else {
    object_name <- paste(object_name, details, sep = "_")
  }

  # create a seed sequence if it's not provided
  if (is.null(seed_sequence)) {
    seed_sequence <- seq(
      from = 1,
      by = 100,
      length.out = n_repetitions
    )
  }

  # calculate the adjacency matrix of the nn or snn graph, depending on the
  # graph type
  adj_matrix <- Seurat::FindNeighbors(
    embedding,
    k.param = n_neigh,
    nn.method = "rann",
    compute.SNN = (graph_type != 0),
    verbose = F
  )[[ifelse(graph_type == 0, "nn", "snn")]]

  # generate the list of partitions
  different_partitions <- get_resolution_partitions(
    adj_matrix,
    resolution = resolution,
    clustering_function = seurat_clustering,
    algorithm = clustering_method,
    seed_sequence = seed_sequence,
    ncores = ncores,
    ecs_thresh = ecs_thresh
  )

  ret_object <- list()
  ret_object[[object_name]] <- list()
  ret_object[[object_name]][[as.character(round(resolution, digits = 3))]] <- different_partitions

  ret_object
}

#' Merge Resolutions
#'
#' @description Merge partitions obtained with different resolution values
#' that have the same number of clusters.
#'
#' @param res_obj A list associated to a configuration field from the object
#' returned by the `get_resolution_importance` method.
#' @param ncores The number of parallel R instances that will run the code.
#' If the value is set to 1, the code will be run sequentially.
#'
#'
#' @return A list having one field assigned to each number of clusters. A number
#' of cluster will contain a list of all merged partitions. To avoid duplicates,
#' `merged_partitions` with threshold 1 is applied.
#'
#' @md
#' @keywords internal
merge_resolutions <- function(res_obj,
                              object_name,
                              ncores = 1) {
  if (!is.numeric(ncores) || length(ncores) > 1) {
    stop("ncores parameter should be numeric")
  }

  # convert ncores to an integer
  ncores <- as.integer(ncores)

  clusters_obj <- list()
  indices <- list()

  # concatenate all partitions having the same number of clusters into a list
  for (res.val in names(res_obj)) {
    for (k in names(res_obj[[res.val]])) {
      if (!(k %in% names(clusters_obj))) {
        clusters_obj[[k]] <- res_obj[[res.val]][[k]]$partitions
      } else {
        clusters_obj[[k]] <- c(clusters_obj[[k]], res_obj[[res.val]][[k]]$partitions)
      }
    }
  }

  ncores <- min(ncores, length(clusters_obj), parallel::detectCores())
  k_vals <- names(clusters_obj)

  if (ncores > 1) {
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

  all_vars <- ls()
  needed_vars <- c("clusters_obj")

  i <- 1

  # merge identical partitions from each list associated with a number of clusters
  clusters_obj <- foreach::foreach(
    i = 1:length(clusters_obj),
    .noexport = all_vars[!(all_vars %in% needed_vars)],
    .export = c("merge_identical_partitions", "are_identical_memberships")
  ) %dopar% {
    merge_identical_partitions(clusters_obj[[i]])
  }

  names(clusters_obj) <- k_vals

  # if a parallel backend was created, terminate the processes
  if (ncores > 1) {
    parallel::stopCluster(cl = my_cluster)
  }

  for (k in names(clusters_obj)) {
    clusters_obj[[k]] <- list(partitions = clusters_obj[[k]])

    # compute the EC-consistency of the partition list
    ec_consistency <- weighted_element_consistency(lapply(clusters_obj[[k]]$partitions, function(x) {
      x$mb
    }),
    sapply(clusters_obj[[k]]$partitions, function(x) {
      x$freq
    }),
    ncores = ncores
    )

    clusters_obj[[k]][["ecc"]] <- ec_consistency
  }

  return_object <- list()
  return_object[[object_name]] <- clusters_obj
  return_object
}

#' Correspondence Between Resolution and the Number of Clusters
#'
#' @description For each configuration provided in the res_object_list, display what
#' number of clusters appear for different values of the resolution parameters.
#'
#'
#' @param res_object_list An object returned by the
#' `get_resolution_importance` method.
#' @param res_object_names Custom names that the user could assing to each
#' configuration; if not specified, the plot will use the generated configuration
#' names.
#' @param given_height Used for adjusting the vertical position of the boxplot; the value
#' will be passed in the `width` argument of the `ggplot::position_dodge` method.
#' @param colour_information String that specifies the information type that will be
#' illustrated using gradient colour: either `freq_k` for the frequency of the number
#' of clusters or `ecc` for the Element-Centric Consistency
#' of the partitions obtained when the resolution is fixed.
#' @param pt_size_range Indicates the minimum and the maximum size a point on the plot can have.
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
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(runif(500 * 10), nrow = 500)
#'
#' # get the PCA embedding of the data
#' pca_embedding <- irlba::irlba(expr_matrix, nv = 2)
#' pca_embedding <- pca_embedding$u %*% diag(pca_embedding$d)
#' rownames(pca_embedding) <- as.character(1:500)
#'
#' # run the function on the pca embedding
#' resolution_result <- get_resolution_importance(
#'   embedding = pca_embedding,
#'   resolution = c(0.8, 1),
#'   n_neigh = c(5, 7),
#'   n_repetitions = 5,
#'   clustering_method = 1,
#'   graph_type = 2,
#'   object_name = "name_example"
#' )
#'
#' plot_k_resolution_corresp(resolution_result)
plot_k_resolution_corresp <- function(res_object_list,
                                      res_object_names = NULL,
                                      colour_information = c("ecc", "freq_k"),
                                      given_height = 0.7,
                                      pt_size_range = c(1.5, 4)) {
  if (length(colour_information) > 1) {
    colour_information <- colour_information[1]
  }

  if (!(colour_information %in% c("ecc", "freq_k"))) {
    stop("colour_information can be either `ecc` or `freq_k`")
  }

  res_object_list <- res_object_list$split_by_resolution

  # use the names of the fields from the list
  if (is.null(res_object_names)) {
    res_object_names <- names(res_object_list)
  }

  object_number <- length(res_object_list)

  # create a dataframe that contains the number of cases when, for a given resolution,
  # a number of clusters was obtained
  for (i in 1:object_number) {
    res_object <- res_object_list[[i]]

    n.runs <- sum(sapply(res_object[[names(res_object)[1]]], function(x) {
      sum(sapply(x$partitions, function(y) {
        y$freq
      }))
    }))

    list.appereances <- lapply(res_object, function(x) {
      lapply(x, function(y) {
        y$partitions[[1]]$freq
      })
    })
    temp.df.appereances <- reshape2::melt(list.appereances)
    colnames(temp.df.appereances) <- c("freq_partition", "number_clusters", "resolution_value")

    temp.df.appereances[["freq_k"]] <- unlist(lapply(res_object, function(x) {
      lapply(x, function(y) {
        sum(sapply(y$partitions, function(z) {
          z$freq
        }))
      })
    }))

    temp.df.appereances[["configuration"]] <- rep(res_object_names[i], nrow(temp.df.appereances))
    temp.df.appereances$freq_partition <- temp.df.appereances$freq_partition / temp.df.appereances$freq_k
    temp.df.appereances$freq_k <- temp.df.appereances$freq_k / n.runs
    temp.df.appereances$ecc <- unlist(lapply(res_object, function(x) {
      sapply(x, function(k) {
        mean(k$ecc)
      })
    }))

    if (i == 1) {
      final.df <- temp.df.appereances
    } else {
      final.df <- rbind(final.df, temp.df.appereances)
    }
  }

  final.df[["configuration"]] <- factor(final.df[["configuration"]])
  final.df[["number_clusters"]] <- factor(as.numeric(final.df[["number_clusters"]]))

  p1 <- ggplot2::ggplot(
    final.df,
    ggplot2::aes(
      y = .data$number_clusters,
      x = .data$resolution_value,
      shape = .data$configuration,
      group = interaction(.data$number_clusters, .data$configuration)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_bw()

  p2 <- ggplot2::ggplot(
    final.df,
    ggplot2::aes(
      y = .data$number_clusters,
      x = .data$resolution_value,
      size = .data$freq_partition,
      group = interaction(.data$number_clusters, .data$configuration)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(size = "partition\nfrequency")

  p3 <- ggplot2::ggplot(
    final.df,
    ggplot2::aes(
      y = .data$number_clusters,
      x = .data$resolution_value,
      color = .data$freq_k,
      group = interaction(.data$number_clusters, .data$configuration)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(name = ifelse(colour_information == "ecc", "ECC", "frequency_k")) +
    ggplot2::theme_bw()

  l1 <- gtable::gtable_filter(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p1)), "guide-box")
  l2 <- gtable::gtable_filter(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p2)), "guide-box")
  l3 <- gtable::gtable_filter(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p3)), "guide-box")
  legends_width <- max(as.numeric(l1$widths), as.numeric(l2$widths) + as.numeric(l3$widths))
  design <- "AA\nBC"
  legends_plot <- patchwork::wrap_plots(l1, l2, l3, design = design)

  main_plot <- ggplot2::ggplot(
    final.df,
    ggplot2::aes(
      y = .data$number_clusters,
      x = .data$resolution_value,
      size = .data$freq_partition,
      color = .data[[colour_information]],
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
    ggplot2::geom_point(position = ggplot2::position_dodge(width = given_height)) + # ggstance::position_dodgev(height = given.height),
    ggplot2::theme_classic() +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(
      x = "resolution",
      y = "k"
    ) +
    ggplot2::scale_size_continuous(range = pt_size_range) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ), legend.position = "none")

  return(patchwork::wrap_plots(main_plot, legends_plot, widths = ggplot2::unit(c(1, legends_width), c("null", "cm"))))
}

#' Relationship Between the Number of Clusters and the Number of Unique Partitions
#'
#' @description For each configuration provided in partition_obj_list, display how
#' many different partitions with the same number of clusters can be obtained
#' by changing the seed.
#'
#' @param partition_obj_list An object or a concatenation of objects returned by the
#' `merge_resolutions` method.
#' @param object_names Custom names that the user could assing to each
#' configuration; if not specified, the plot will use the generated configuration
#' names.
#' @param colour_information String that specifies the information type that will be
#' illustrated using gradient colour: either `frequency_partition` for the frequency of the
#' most common partition or `ecc` for the Element-Centric Consistency
#' of the partitions obtained when the the number of clusters is fixed.
#' @param pt_size_range Indicates the minimum and the maximum size a point on the plot can have.
#'
#' @return A ggplot2 object. The color gradient suggests the frequency of the most
#' common partition relative to the total number of appearances of that specific
#' number of clusters or the Element-Centric Consistency of the partitions. The size
#' illustrates the frequency of the partitions with *k* clusters relative to the
#' total number of partitions.
#' @export
#'
#'
#' @examples
#' set.seed(2021)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(runif(500 * 10), nrow = 500)
#'
#' # get the PCA embedding of the data
#' pca_embedding <- irlba::irlba(expr_matrix, nv = 2)
#' pca_embedding <- pca_embedding$u %*% diag(pca_embedding$d)
#' rownames(pca_embedding) <- as.character(1:500)
#'
#' # run the function on the pca embedding
#' resolution_result <- get_resolution_importance(
#'   embedding = pca_embedding,
#'   resolution = c(0.8, 1),
#'   n_neigh = c(5, 7),
#'   n_repetitions = 5,
#'   clustering_method = 1,
#'   graph_type = 2,
#'   object_name = "name_example"
#' )
#'
#' plot_k_n_partitions(resolution_result)
plot_k_n_partitions <- function(partition_obj_list,
                                object_names = NULL,
                                pt_size_range = c(1.5, 4),
                                colour_information = c("ecc", "frequency_partition")) {
  if (length(colour_information) > 1) {
    colour_information <- colour_information[1]
  }

  if (!(colour_information %in% c("ecc", "frequency_partition"))) {
    stop("colour_information can be either `ecc` or `frequency_partition`")
  }

  partition_obj_list <- partition_obj_list$split_by_k

  # use the names of the fields from the list
  if (is.null(object_names)) {
    object_names <- names(partition_obj_list)
  }

  n.objects <- length(partition_obj_list)

  max_n_part <- 0

  # creates a dataframe that contains, for each configuration
  # the number of different partitions with a given number of clusters that are obtained
  for (i in 1:n.objects) {
    partition_object <- partition_obj_list[[i]]

    unique.partitions.temp.df <- reshape2::melt(lapply(partition_object, function(x) {
      length(x$partitions)
    }))
    colnames(unique.partitions.temp.df) <- c("n.partitions", "n.clusters")

    unique.partitions.temp.df[["configuration"]] <- rep(object_names[i], nrow(unique.partitions.temp.df))

    unique.partitions.temp.df[["first.occ"]] <- as.numeric(lapply(partition_object, function(x) {
      max(sapply(x$partitions, function(y) {
        y$freq
      }))
    }))
    unique.partitions.temp.df[["total.occ"]] <- as.numeric(lapply(partition_object, function(x) {
      sum(sapply(x$partitions, function(y) {
        y$freq
      }))
    }))
    unique.partitions.temp.df[["frequency_partition"]] <- unique.partitions.temp.df$first.occ / unique.partitions.temp.df$total.occ
    unique.partitions.temp.df[["ecc"]] <- sapply(partition_object, function(x) {
      mean(x$ecc)
    })
    overall_total_occ <- sum(unique.partitions.temp.df$total.occ)
    unique.partitions.temp.df[["frequency_k"]] <- unique.partitions.temp.df$total.occ / overall_total_occ

    max_n_part <- max(c(max(
      unique.partitions.temp.df$n.partitions
    ), max_n_part))

    if (i == 1) {
      unique.partitions.final.df <- unique.partitions.temp.df
    } else {
      unique.partitions.final.df <- rbind(unique.partitions.final.df, unique.partitions.temp.df)
    }
  }

  unique.partitions.final.df$configuration <- factor(unique.partitions.final.df$configuration)
  unique.partitions.final.df$n.clusters <- factor(unique.partitions.final.df$n.clusters,
    levels = stringr::str_sort(unique(
      unique.partitions.final.df$n.clusters
    ), numeric = T)
  )
  p1 <- ggplot2::ggplot(
    unique.partitions.final.df,
    ggplot2::aes(
      x = .data$n.clusters,
      y = .data$n.partitions,
      shape = .data$configuration,
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_bw()

  p2 <- ggplot2::ggplot(
    unique.partitions.final.df,
    ggplot2::aes(
      x = .data$n.clusters,
      y = .data$n.partitions,
      size = .data$frequency_k
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_bw()

  p3 <- ggplot2::ggplot(
    unique.partitions.final.df,
    ggplot2::aes(
      x = .data$n.clusters,
      y = .data$n.partitions,
      color = .data[[colour_information]]
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(name = ifelse(colour_information == "ecc", "ECC", "partition\nfrequency")) +
    ggplot2::theme_bw()

  l1 <- gtable::gtable_filter(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p1)), "guide-box")
  l2 <- gtable::gtable_filter(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p2)), "guide-box")
  l3 <- gtable::gtable_filter(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p3)), "guide-box")
  legends_width <- max(as.numeric(l1$widths), as.numeric(l2$widths) + as.numeric(l3$widths))
  design <- "AA\nBC"
  legends_plot <- patchwork::wrap_plots(l1, l2, l3, design = design)

  main_plot <- ggplot2::ggplot(
    unique.partitions.final.df,
    ggplot2::aes(
      x = .data$n.clusters,
      y = .data$n.partitions,
      shape = .data$configuration,
      size = .data$frequency_k,
      color = .data[[colour_information]]
    ),
    group = interaction(.data$n.clusters, .data$precision)
  ) +
    ggplot2::scale_y_continuous(breaks = seq(from = 0, to = max_n_part, by = 2)) +
    ggplot2::scale_size_continuous(range = pt_size_range) +
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
    ggplot2::geom_point(position = ggplot2::position_dodge(0.9)) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_viridis_c(name = ifelse(colour_information == "ecc", "ECC", "partition\nfrequency")) +
    # ggplot2::theme()
    ggplot2::xlab("k") +
    ggplot2::ylab("# partitions") +
    ggplot2::theme(legend.position = "none")

  return(patchwork::wrap_plots(main_plot, legends_plot, widths = ggplot2::unit(c(1, legends_width), c("null", "cm"))))
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


automatic_clustassess <- function(expression_matrix, # expr matrix
                                  n_cores,
                                  n_repetitions,
                                  n_neigh_sequence,
                                  resolution_sequence,
                                  seed_sequence = NULL,
                                  features_sets, # list with the names of the feature sets (HV, MA etc) and the list of top x genes
                                  steps, # list with the names of the feature sets (HV, mA) and the list of steps
                                  graph_reduction_embedding = "PCA",
                                  n_top_configs = 3,
                                  approx = FALSE,
                                  ranking_criterion = "median",
                                  npcs = 30,
                                  ecs_threshold = 1, # do we really need it?,
                                  ...) {
  # store the additional arguments used by umap in a list
  suppl_args <- list(...)

  feature_set_names <- names(steps)
  config_names <- paste(names(steps), graph_reduction_embedding, ecs_threshold, sep = "_")
  n_configs <- length(config_names)
  # to-do: check if the name of the feature sets are the same with the steps


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
    temp_object <- get_feature_stability(
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
      suppressWarnings(pca_emb <- RunPCA(expression_matrix[current_features, ],
        npcs = min(npcs, n_steps %/% 2),
        approx = approx,
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
      feature_configs[[config_name]][[n_steps]][["nn_conn_comps"]] <- c(
        get_nn_conn_comps(
          object = feature_configs[[config_name]][[n_steps]]$pca,
          n_neigh_sequence = n_neigh_sequence,
          n_repetitions = n_repetitions,
          graph_reduction_type = "UMAP",
          ncores = n_cores,
          seed_sequence = seed_sequence,
          config_name = paste(config_name, n_steps, sep = "_"),
          ...
        ),
        get_nn_conn_comps(
          object = expression_matrix,
          n_neigh_sequence = n_neigh_sequence,
          n_repetitions = n_repetitions,
          graph_reduction_type = "PCA",
          seed_sequence = seed_sequence,
          config_name = paste(config_name, n_steps, sep = "_"),
          ncores = n_cores,
          nv = min(as.numeric(n_steps) %/% 2, npcs)
        )
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
        get_nn_importance(
          object = expression_matrix[current_features, ],
          n_neigh_sequence = n_neigh_sequence,
          n_repetitions = n_repetitions,
          graph_reduction_type = "PCA",
          ecs_thresh = 1,
          ncores = n_cores,
          algorithm = 1,
          config_name = paste(config_name, n_steps, sep = "_"),
          nv = min(as.numeric(n_steps) %/% 2, npcs)
        ),
        get_nn_importance(
          object = cuomo@reductions$pca@cell.embeddings,
          n_neigh_sequence = n_neigh_sequence,
          n_repetitions = n_repetitions,
          graph_reduction_type = "UMAP",
          ecs_thresh = 1,
          ncores = n_cores,
          config_name = paste(config_name, n_steps, sep = "_"),
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
          if (current_ecc >= best_ecc) {
            best_ecc <- current_ecc
            best_config <- config_name
            best_nn <- as.numeric(n_neigh)
          }
        }
      }

      split_configs <- strsplit(best_config, "_")[[1]]
      base_embedding <- tolower(split_configs[length(split_configs) - 2])
      graph_type <- split_configs[length(split_configs) - 1] # update if you decide to remove ecs thresh
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
      feature_configs[[config_name]][[n_steps]][["clustering_importance"]] <- get_clustering_difference(
        graph_adjacency_matrix = feature_configs[[config_name]][[n_steps]]$adj_matrix,
        resolution = resolution_sequence,
        n_repetitions = n_repetitions,
        ecs_thresh = ecs_threshold,
        ncores = n_cores,
        algorithm = 1:4
      )

      ecc_clustering_methods <- lapply(
        feature_configs[[config_name]][[n_steps]]$clustering_importance$all,
        function(clust_method) {
          c(sapply(clust_method, c))
        }
      )

      best_clustering_method <- rank_configs(ecc_clustering_methods, rank_by = ranking_criterion)[1]
      feature_configs[[config_name]][[n_steps]]$stable_config[["clustering_method"]] <- algorithm_names[best_clustering_method]
    }
  }

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
      temp_list = list()

      for(k in names(feature_configs[[feature_config_name]][[n_steps]]$resolution_importance$split_by_k[[1]])) {
        temp_list[[k]] = feature_configs[[feature_config_name]][[n_steps]]$resolution_importance$split_by_k[[1]][[k]]$partitions[[1]]$mb
      }

      feature_configs[[feature_config_name]][[n_steps]][["stable_partitions"]] = temp_list
    }
  }

  return(feature_configs)
}
