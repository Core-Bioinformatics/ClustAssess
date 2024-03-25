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
