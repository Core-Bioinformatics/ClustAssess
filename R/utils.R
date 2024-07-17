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

# source: https://github.com/satijalab/seurat/blob/763259d05991d40721dee99c9919ec6d4491d15e/R/objects.R#L2415
get_count_features <- function(counts_matrix) {
    return(list(
        nCount = Matrix::colSums(counts_matrix),
        nFeature = Matrix::colSums(counts_matrix > 0)
    ))
}

process_umap_arguments <- function(umap_arguments, ncells) {
    umap_arguments[["n_threads"]] <- 1
    umap_arguments[["n_sgd_threads"]] <- 1

    if (!("n_neighbors" %in% umap_arguments)) {
        umap_arguments[["n_neighbors"]] <- 15 # uwot's default
    }

    umap_arguments[["n_neighbors"]] <- min(umap_arguments[["n_neighbors"]], ncells - 1)
    return(umap_arguments)
}

process_clustering_arguments <- function(clustering_arguments, clustering_method) {
    num_starts <- 1
    num_iters <- 10
    if ("n.start" %in% names(clustering_arguments)) {
        num_starts <- clustering_arguments$n.start
        clustering_arguments$n.start <- NULL
    }

    if ("num_starts" %in% names(clustering_arguments)) {
        num_starts <- clustering_arguments$num_starts
        clustering_arguments$num_starts <- NULL
    }

    if ("n.iter" %in% names(clustering_arguments)) {
        num_iters <- clustering_arguments$n.iter
        clustering_arguments$n.iter <- NULL
    }

    if ("num_iter" %in% names(clustering_arguments)) {
        num_iters <- clustering_arguments$num_iter
        clustering_arguments$num_iter <- NULL
    }

    if (!("verbose" %in% names(clustering_arguments))) {
        clustering_arguments$verbose <- FALSE
    }

    clustering_arguments$num_starts <- num_starts
    clustering_arguments$num_iters <- num_iters

    clustering_arguments$algorithm <- clustering_method

    return(clustering_arguments)
}

#' Extract config-specific clusters from a ClustAssess object
#'
#' @description Given the output of the `automatic_stability_assessment`
#' function, extract the clusters that are specific to a particular
#' configuration of feature type, feature size, clustering method and,
#' optionally, the number of clusters.
#'
#' @param clustassess_object Output of the `automatic_stability_assessment`.
#' @param feature_type Type of feature used for dimensionality reduction.
#' @param feature_size Size of the feature set used for clustering.
#' @param clustering_method Clustering method used.
#' @param nclusters Number of clusters to extract. If NULL, all clusters are
#' returned.
#'
#' @return A list of clusters that are specific to the given configuration.
#' Each number of cluster will contain the list of partitions with that
#' specific k and the ECC value indicating the overall stability of k.
#'
#' @export
get_clusters_from_clustassess_object <- function(clustassess_object,
                                                 feature_type,
                                                 feature_size,
                                                 clustering_method,
                                                 nclusters = NULL) {

    if (!(feature_type %in% names(clustassess_object))) {
        available_options <- setdiff(names(clustassess_object), "feature_stability")
        stop(glue::glue("Feature type not found in clustassess object.\nAvailable options: {available_options}"))
    }

    clustassess_object <- clustassess_object[[feature_type]]
    if (!(feature_size %in% names(clustassess_object))) {
        available_options <- paste(setdiff(names(clustassess_object), "feature_list"), collapse = ", ")
        stop(glue::glue("Feature size not found in clustassess object.\nAvailable options: {available_options}"))
    }

    feature_size <- as.character(feature_size)

    clustassess_object <- clustassess_object[[feature_size]]$clustering_stability$split_by_k
    if (!(clustering_method %in% names(clustassess_object))) {
        available_options <- paste(names(clustassess_object), collapse = ", ")
        stop(glue::glue("Clustering method not found in clustassess object.\nAvailable options: {available_options}"))
    }

    clustassess_object <- clustassess_object[[clustering_method]]

    if (is.null(nclusters)) {
        return(clustassess_object)
    }

    nclusters <- as.character(nclusters)
    nclusters <- intersect(nclusters, names(clustassess_object))

    if (length(nclusters) == 0) {
        available_options <- paste(names(clustassess_object), collapse = ", ")
        stop(glue::glue("Number of clusters not found in clustassess object.\nAvailable options: {available_options}"))
    }

    return(clustassess_object[nclusters])
}

#' Choose stable clusters based on ECC and frequency
#'
#' @description Filter the list of clusters obtained by the automatic
#' ClustAssess pipeline using the ECC and frequency thresholds. The ECC
#' threshold is meant to filter out the partitions that are highly sensitive
#' to the change of the random seed, while the purpose of the frequency
#' threshold is to assure a statistical significance of the inferred stability.
#'
#' @param clusters_list List of clusters obtained from the
#' `get_clusters_from_clustassess_object` function.
#' @param ecc_threshold Minimum ECC value to consider a cluster as stable.
#' Default is 0.9.
#' @param freq_threshold Minimum total frequency of the partitions to consider.
#' Default is 30.
#' @param summary_function Function to summarize the ECC values. Default
#' is `mean`. To match the results from the ClustAssess Shiny App, use
#' `median`.
#'
#' @return A list of stable clusters that satisfy the ECC and frequency.
#'
#' @export
choose_stable_clusters <- function(clusters_list,
                                   ecc_threshold = 0.9,
                                   freq_threshold = 30,
                                   summary_function = mean) {

    if (!all(c("ecc", "partitions") %in% names(clusters_list[[1]]))) {
        stop("Invalid clusters list. Please run `get_clusters_from_clustassess_object` first on the ClustAssess object.")
    }

    ecc_vals <- sapply(clusters_list, function(x) summary_function(x$ecc))
    total_freq <- sapply(clusters_list, function(x) {
        sum(sapply(x$partitions, function(y) y$freq))
    })

    index_ecc <- which(ecc_vals >= ecc_threshold)
    index_freq <- which(total_freq >= freq_threshold)

    index_chosen <- intersect(index_ecc, index_freq)
   
    if (length(index_chosen) == 0) {
        stop("No stable clusters found. Consider lowering the `ecc_threshold` or `freq_threshold`.")
    }

    return(clusters_list[index_chosen])
}