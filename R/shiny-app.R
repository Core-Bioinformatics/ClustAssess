

# stab_obj <- readRDS("stability_object2.rds")
# expression_matrix <- readRDS("expression_matrix.rds")
# metadata <- readRDS("metadata.rds")

# write_objects(stab_obj, expression_matrix, metadata, compression_level = 6)


write_objects <- function(clustassess_object,
                          expression_matrix,
                          metadata,
                          project_folder = ".",
                          compression_level = 6) {
    metadata_file_name <- file.path(project_folder, "metadata.rds")
    stability_file_name <- file.path(project_folder, "stability.h5")
    expr_file_name <- file.path(project_folder, "expression.h5")

    if (file.exists(expr_file_name)) {
        stop(glue::glue("File {expr_file_name} already exists!"))
    }

    if (file.exists(stability_file_name)) {
        stop(glue::glue("File {stability_file_name} already exists!"))
    }

    # if (file.exists(metadata_file_name)) {
    #     stop(glue::glue("File {metadata_file_name} already exists!"))
    # }

    # metadata file
    saveRDS(metadata, metadata_file_name)

    feature_ordering <- list(original = list(), stable = list(), resolution = list())

    for (ftype in names(stab_obj$feature_stability$by_steps)) {
        feature_ordering$original[[ftype]] <- names(stab_obj$feature_stability$by_steps[[ftype]])
        feature_ordering$resolution[[ftype]] <- list()
        for (fsize in names(stab_obj$feature_stability$by_steps[[ftype]])) {
            feature_ordering$resolution[[ftype]][[fsize]] <- names(stab_obj$feature_stability$by_steps[[ftype]][[fsize]])
            for (res in names(stab_obj$feature_stability$by_steps[[ftype]][[fsize]])) {
                stab_obj$feature_stability$by_steps[[ftype]][[fsize]][[res]]$mb <- as.integer(stab_obj$feature_stability$by_steps[[ftype]][[fsize]][[res]]$most_frequent_partition$mb)
                stab_obj$feature_stability$by_steps[[ftype]][[fsize]][[res]]$most_frequent_partition <- NULL
                stab_obj$feature_stability$by_steps[[ftype]][[fsize]][[res]]$n_different_partitions <- NULL
            }
        }
    }

    for (ftype in names(feature_ordering$original)) {
        nsteps <- length(stab_obj[[ftype]]) - 1
        feature_ordering$stable[[ftype]] <- names(stab_obj[[ftype]])[seq_len(nsteps)]
    }

    for (ftype in names(feature_ordering$stable)) {
        for (fsize in feature_ordering$stable[[ftype]]) {
            # remove the adjacency matrix, not needed for the shiny app
            stab_obj[[ftype]][[fsize]]$adj_matrix <- NULL

            # remove the list of unecessary partitions
            cl_methods <- names(stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k)

            for (cl_method in cl_methods) {
                # split by resolution
                for (res in names(stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]])) {
                    for (n_clusters in names(stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters)) {
                        stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$mb <- as.integer(stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions[[1]]$mb)
                        stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$first_freq <- stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions[[1]]$freq
                        stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$total_freq <- sum(
                            sapply(
                                stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions,
                                function(x) { x$freq }
                            )
                        )
                        stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions <- NULL
                    }
                }

                # split by k
                for (n_clusters in names(stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]])) {
                    stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$mb <- as.integer(stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions[[1]]$mb)
                    stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$first_freq <- stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions[[1]]$freq
                    stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$total_freq <- sum(
                        sapply(
                            stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions,
                            function(x) { x$freq }
                        )
                    )
                    stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$n_partitions <- length(stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions)
                    stab_obj[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions <- NULL
                }
            }
        }
    }

    stab_obj$feature_ordering <- feature_ordering
    feature_names <- names(feature_ordering$original)
    genes_of_interest <- c()
    for (feature_name in feature_names) {
        genes_of_interest <- union(genes_of_interest, stab_obj[[feature_name]]$feature_list)
    }

    others_sorted <- sort(rownames(expression_matrix)[!(rownames(expression_matrix) %in% genes_of_interest)])
    genes <- c(genes_of_interest, others_sorted)
    cells <- colnames(expression_matrix)

    # the object for expression matrix
    rhdf5::h5createFile(expr_file_name)

    rhdf5::h5createDataset(expr_file_name, "matrix_of_interest",
        level = compression_level,
        dims = c(length(genes_of_interest), ncol(expression_matrix)),
        chunk = c(length(genes_of_interest), 1)
    )
    rhdf5::h5createDataset(expr_file_name, "matrix_others",
        level = compression_level,
        dims = c(length(others_sorted), ncol(expression_matrix)),
        chunk = c(length(others_sorted), 1)
    )

    rhdf5::h5write(genes_of_interest, expr_file_name, "genes_of_interest")
    rhdf5::h5write(others_sorted, expr_file_name, "genes_others")
    rhdf5::h5write(cells, expr_file_name, "cells")
    rhdf5::h5write(as.matrix(expression_matrix[genes_of_interest, ]), expr_file_name, "matrix_of_interest")
    rhdf5::h5write(as.matrix(expression_matrix[others_sorted, ]), expr_file_name, "matrix_others")

    # the object for stability
    rhdf5::h5createFile(stability_file_name)
    rhdf5::h5write(stab_obj$feature_stability, stability_file_name, "feature_stability", level = compression_level)
    rhdf5::h5write(stab_obj$feature_ordering, stability_file_name, "feature_ordering", level = compression_level)

    for (feature_type in names(stab_obj$feature_ordering[[1]])) {
        rhdf5::h5write(stab_obj[[feature_type]], stability_file_name, feature_type, level = compression_level)
    }

    rhdf5::h5closeAll()
}
