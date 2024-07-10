
#' Create monocle object
#'
#' @description Use a normalized expression matrix and, potentially, an already
#' generated PCA / UMAP embedding, to create a Monocle object.
#'
#' @param normalized_expression_matrix The normalized expression matrix
#' having genes on rows and cells on columns.
#' @param count_matrix The count matrix having genes on rows and cells on
#' columns. If NULL, the normalized_expression_matrix will be used.
#' @param pca_emb The PCA embedding of the expression matrix. If NULL, the
#' pca will be created using the `monocle3` package (default parameters).
#' @param umap_emb The UMAP embedding of the expression matrix. If NULL, the
#' umap will be created using the `monocle3` package (default parameters).
#' @param metadata_df The metadata dataframe having the cell names as rownames.
#' If NULL, a dataframe with a single column named `identical_ident` will be
#' created.
#'
#' @return A Monocle object of the expression matrix, having the stable number
#' of clusters identified by ClustAssess.
#'
#' @export
#' @examples
#' set.seed(2024)
#' # create an already-transposed artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(20 * 10), runif(30 * 10, min = 3, max = 4)),
#'     nrow = 10, byrow = FALSE
#' )
#' colnames(expr_matrix) <- as.character(seq_len(ncol(expr_matrix)))
#' rownames(expr_matrix) <- paste("feature", seq_len(nrow(expr_matrix)))
#'
#' # uncomment to create the monocle object
#' mon_obj <- create_monocle_default(
#'     normalized_expression_matrix = expr_matrix,
#'     pca_emb = NULL,
#'     umap_emb = NULL,
#'     metadata_df = NULL
#' )
create_monocle_default <- function(normalized_expression_matrix,
                                   count_matrix = NULL,
                                   pca_emb = NULL,
                                   umap_emb = NULL,
                                   metadata_df = NULL) {
    if (is.null(count_matrix)) {
        count_matrix <- normalized_expression_matrix
    } else {
        count_matrix <- count_matrix[rownames(normalized_expression_matrix), colnames(normalized_expression_matrix)]
    }

    gene_names <- rownames(normalized_expression_matrix)
    cell_names <- colnames(normalized_expression_matrix)

    if (is.null(metadata_df)) {
        metadata_df <- data.frame(
            one_level = rep("one", ncol(normalized_expression_matrix)),
            row.names = cell_names
        )
    } else {
        if (nrow(metadata_df) != ncol(normalized_expression_matrix)) {
            stop("The number of rows in metadata_df should be equal to the number of cells in normalized_expression_matrix")
        }

        rownames(metadata_df) <- cell_names
    }


    monocle_cds <- monocle3::new_cell_data_set(
        expression_data = count_matrix,
        cell_metadata = metadata_df,
        gene_metadata = data.frame("gene_short_name" = gene_names, row.names = gene_names)
    )

    monocle_cds <- monocle3::preprocess_cds(
        cds = monocle_cds,
        num_dim = ifelse(is.null(pca_emb), 30, ncol(pca_emb)),
        norm_method = "none",
        scaling = FALSE
    )

    monocle_cds@assays@data$normalized_data <- normalized_expression_matrix

    if (!is.null(pca_emb)) {
        rownames(pca_embedding) <- cell_names
        colnames(pca_embedding) <- paste0("PC_", seq_len(ncol(pca_embedding)))
        monocle_cds@int_colData$reducedDims$PCA <- pca_embedding
    }

    if (!is.null(umap_emb)) {
        rownames(umap_embedding) <- cell_names
        colnames(umap_embedding) <- paste0("UMAP_", seq_len(ncol(umap_embedding)))
        monocle_cds@int_colData$reducedDims$UMAP <- umap_embedding
    }

    return(monocle_cds)
}

#' Create monocle object from a ClustAssess object
#'
#' @description Use the object generated using the ClustAssess
#' `automatic_stability_assessment` function to create a Monocle object
#' which has the stable number of clusters.
#'
#' @param normalized_expression_matrix The normalized expression matrix
#' having genes on rows and cells on columns.
#' @param count_matrix The count matrix having genes on rows and cells on
#' columns. If NULL, the normalized_expression_matrix will be used.
#' @param clustassess_object The output of the `automatic_stability_assessment`.
#' @param metadata_df The metadata dataframe having the cell names as rownames.
#' If NULL, a dataframe with a single column named `identical_ident` will be
#' created.
#' @param stable_feature_type The feature type which leads to stable clusters.
#' @param stable_feature_set_size The feature size which leads to stable
#' clusters.
#' @param stable_clustering_method The clustering method which leads to stable
#' clusters.
#' @param stable_n_clusters The number of clusters that are stable. If NULL,
#' all the clusters will be provided. Defaults to `NULL`.
#' @param use_all_genes A boolean value indicating if the expression matrix
#' should be truncated to the genes used in the stability assessment. Defaults
#' to `FALSE`.
#'
#' @return A Monocle object of the expression matrix, having the stable number
#' of clusters identified by ClustAssess.
#'
#' @export
#' @examples
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
#'         # from the default values to ensure that the function
#'         # will run under 5 seconds
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
#'
#' # uncomment to create the monocle object
#' # mon_obj <- create_monocle_from_clustassess(
#' #     normalized_expression_matrix = expr_matrix,
#' #     clustassess_object = autom_object,
#' #     metadata = NULL,
#' #     stable_feature_type = "set1",
#' #     stable_feature_set_size = "5",
#' #     stable_clustering_method = "Louvain"
#' # )
create_monocle_from_clustassess <- function(normalized_expression_matrix,
                                            count_matrix = NULL,
                                            clustassess_object,
                                            metadata_df,
                                            stable_feature_type,
                                            stable_feature_set_size,
                                            stable_clustering_method,
                                            stable_n_clusters = NULL,
                                            use_all_genes = FALSE) {
    gene_names <- rownames(normalized_expression_matrix)
    cell_names <- colnames(normalized_expression_matrix)

    if (is.null(gene_names)) {
        stop("The expression matrix should contain gene names as rownames.")
    }

    if (is.null(cell_names)) {
        stop("The expression matrix should contain gene names as colnames.")
    }

    if (!use_all_genes) {
        used_genes <- clustassess_object[[stable_feature_type]]$feature_list[seq_len(as.integer(stable_feature_set_size))]
        normalized_expression_matrix <- normalized_expression_matrix[used_genes, ]
        if (!is.null(count_matrix)) {
            count_matrix <- count_matrix[used_genes, ]
        }

        gene_names <- used_genes
    }

    clustassess_object <- clustassess_object[[stable_feature_type]][[as.character(stable_feature_set_size)]]
    gc()
    available_n_clusters <- names(clustassess_object$clustering_stability$split_by_k[[stable_clustering_method]])

    if (is.null(stable_n_clusters)) {
        stable_n_clusters <- available_n_clusters
    } else {
        stable_n_clusters <- as.character(stable_n_clusters)
        which_present <- which(available_n_clusters %in% stable_n_clusters)

        if (length(which_present) != length(available_n_clusters)) {
            warning(paste0("Only the n clusters of ", paste(available_n_clusters[which_present], collapse = " "), " are obtained for the provided configuration."))
        }

        stable_n_clusters <- available_n_clusters[which_present]
    }

    if (!is.data.frame(metadata_df) || nrow(metadata_df) == 0) {
        metadata_df <- data.frame(
            one_level = rep(1, length(cell_names)),
            row.names = cell_names
        )
    } else {
        metadata_values <- colnames(metadata_df)
        metadata_values[which(metadata_values == "sample_name")] <- "sample_names"
        colnames(metadata_df) <- metadata_values
    }

    for (k in stable_n_clusters) {
        metadata_df[[paste0("stable_", k, "_clusters")]] <- factor(clustassess_object$clustering_stability$split_by_k[[stable_clustering_method]][[k]]$partitions[[1]]$mb)
        metadata_df[[paste0("ecc_", k)]] <- clustassess_object$clustering_stability$split_by_k[[stable_clustering_method]][[k]]$ecc
    }

    return(create_monocle_default(
        normalized_expression_matrix,
        count_matrix,
        clustassess_object$pca,
        clustassess_object$umap,
        metadata_df
    ))
}

#' Create monocle object from a ClustAssess shiny app
#'
#' @description Use the files generated in the ClustAssess app to create a
#' Monocle object which has the stable number of clusters.
#'
#' @param app_folder Path pointing to the folder containing a ClustAssess app.
#' @param stable_feature_type The feature type which leads to stable clusters.
#' @param stable_feature_set_size The feature size which leads to stable
#' clusters.
#' @param stable_clustering_method The clustering method which leads to stable
#' clusters.
#' @param stable_n_clusters The number of clusters that are stable. If NULL,
#' all the clusters will be provided. Defaults to `NULL`.
#' @param use_all_genes A boolean value indicating if the expression matrix
#' should be truncated to the genes used in the stability assessment. Defaults
#' to `FALSE`.
#'
#' @return A Monocle object of the expression matrix, having the stable number
#' of clusters identified by ClustAssess.
#'
#' @export
create_monocle_from_clustassess_app <- function(app_folder,
                                                stable_feature_type,
                                                stable_feature_set_size,
                                                stable_clustering_method,
                                                stable_n_clusters = NULL,
                                                use_all_genes = FALSE) {
    if (!dir.exists(app_folder)) {
        stop(paste0("The provided app_folder: ", app_folder, " does not exist."))
    }

    expr_path <- file.path(app_folder, "expression.h5")
    stab_path <- file.path(app_folder, "stability.h5")
    gene_names <- rhdf5::h5read(expr_path, "genes")
    cell_names <- rhdf5::h5read(expr_path, "cells")

    expr_matrix <- rhdf5::h5read(expr_path, "expression")
    colnames(expr_matrix) <- cell_names
    rownames(expr_matrix) <- gene_names

    if (!use_all_genes) {
        used_genes <- rhdf5::h5read(stab_path, paste0(stable_feature_type, "/feature_list"))
        used_genes <- used_genes[seq_len(as.integer(stable_feature_set_size))]
        expr_matrix <- expr_matrix[used_genes, ]
        gene_names <- used_genes
    }

    prefix <- paste0(stable_feature_type, "/", stable_feature_set_size)
    cl_mb_prefix <- paste0(prefix, "/clustering_stability/split_by_k/mbs/", stable_clustering_method)
    cl_ecc_prefix <- paste0(prefix, "/clustering_stability/split_by_k/ecc/")
    metadata_df <- readRDS(file.path(app_folder, "metadata.rds"))$metadata

    available_n_clusters <- names(rhdf5::h5read(stab_path, cl_mb_prefix))

    if (is.null(stable_n_clusters)) {
        stable_n_clusters <- available_n_clusters
    } else {
        stable_n_clusters <- as.character(stable_n_clusters)
        which_present <- which(available_n_clusters %in% stable_n_clusters)

        if (length(which_present) != length(available_n_clusters)) {
            warning(paste0("Only the n clusters of ", paste(available_n_clusters[which_present], collapse = " "), " are obtained for the provided configuration."))
        }

        stable_n_clusters <- available_n_clusters[which_present]
    }

    metadata_values <- colnames(metadata_df)
    metadata_values[which(metadata_values == "sample_name")] <- "sample_names"
    colnames(metadata_df) <- metadata_values

    for (k in stable_n_clusters) {
        metadata_df[[paste0("stable_", k, "_clusters")]] <- rhdf5::h5read(stab_path, paste0(cl_mb_prefix, "/", k))
        metadata_df[[paste0("ecc_", k)]] <- rhdf5::hread(
            stab_path,
            paste0(cl_ecc_prefix, sprintf("%06d;%s", k, stable_clustering_method))
        )
    }

    return(create_monocle_default(
        expr_matrix,
        NULL,
        rhdf5::h5read(stab_path, paste0(prefix, "/pca")),
        rhdf5::h5read(stab_path, paste0(prefix, "/umap")),
        metadata_df
    ))
}

# TODO create this funciton
update_seurat_object <- function(original_seurat_object,
                                 clustassess_object,
                                 stable_feature_type,
                                 stable_feature_set_size,
                                 stable_clustering_method,
                                 stable_n_clusters,
                                 seurat_assay = "SCT") {


}
