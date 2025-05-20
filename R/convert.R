
get_dim_reduction_from_clustassess_app <- function(app_folder,
                                                   dim_type = c("pca", "umap"),
                                                   dim_key = "PC_",
                                                   stable_feature_type,
                                                   stable_feature_set_size) {
    dim_type <- dim_type[1]

    expr_path <- file.path(app_folder, "expression.h5")
    stab_path <- file.path(app_folder, "stability.h5")

    if (!file.exists(expr_path)) {
        stop(paste0("The expression.h5 file does not exist in the provided app_folder: ", app_folder))
    }

    if (!file.exists(stab_path)) {
        stop(paste0("The stability.h5 file does not exist in the provided app_folder: ", app_folder))
    }

    cell_names <- as.character(rhdf5::h5read(expr_path, "cells"))

    available_configs <- rhdf5::h5read(stab_path, "feature_ordering/stable")
    available_ftypes <- names(available_configs)

    if (!stable_feature_type %in% available_ftypes) {
        stop(paste0("The provided stable_feature_type: ", stable_feature_type, " is not available in the app.\nAvailable options: ", paste(available_ftypes, collapse = ", ")))
    }

    available_fsizes <- available_configs[[stable_feature_type]]

    if (!stable_feature_set_size %in% available_fsizes) {
        stop(paste0("The provided stable_feature_set_size: ", stable_feature_set_size, " is not available in the app.\nAvailable options: ", paste(available_fsizes, collapse = ", ")))
    }

    prefix <- paste0(stable_feature_type, "/", stable_feature_set_size)
    stab_structure <- rhdf5::h5ls(stab_path)
    avail_options <- stab_structure$name[stab_structure$group == paste0("/", prefix)]

    if (!(dim_type %in% avail_options)) {
        stop(paste0("The provided dim_type: ", dim_type, " is not available in the app.\nAvailable options: ", paste(avail_options, collapse = ", ")))
    }

    dim_emb <- rhdf5::h5read(stab_path, paste0(prefix, "/", dim_type))
    rownames(dim_emb) <- cell_names
    colnames(dim_emb) <- paste0(dim_key, seq_len(ncol(dim_emb)))

    return(dim_emb)
}

get_clusters_from_clustassess_app <- function(app_folder,
                                              stable_feature_type,
                                              stable_feature_set_size,
                                              stable_clustering_method) {
    stab_path <- file.path(app_folder, "stability.h5")
    if (!file.exists(stab_path)) {
        stop(paste0("The stability.h5 file does not exist in the provided app_folder: ", app_folder))
    }

    available_configs <- rhdf5::h5read(stab_path, "feature_ordering/stable")
    available_ftypes <- names(available_configs)

    if (!stable_feature_type %in% available_ftypes) {
        stop(paste0("The provided stable_feature_type: ", stable_feature_type, " is not available in the app.\nAvailable options: ", paste(available_ftypes, collapse = ", ")))
    }

    available_fsizes <- available_configs[[stable_feature_type]]

    if (!stable_feature_set_size %in% available_fsizes) {
        stop(paste0("The provided stable_feature_set_size: ", stable_feature_set_size, " is not available in the app.\nAvailable options: ", paste(available_fsizes, collapse = ", ")))
    }

    prefix <- paste0(stable_feature_type, "/", stable_feature_set_size)
    cl_mb_prefix <- paste0(prefix, "/clustering_stability/split_by_k/mbs/", stable_clustering_method)
    cl_ecc_prefix <- paste0(prefix, "/clustering_stability/split_by_k/ecc")

    tryCatch({
        available_n_clusters <- names(rhdf5::h5read(stab_path, cl_mb_prefix))
    }, error = function(e) {
        stop(paste0("The provided stable_clustering_method: ", stable_clustering_method, " is not available in the app."))
    })

    df_list <- lapply(available_n_clusters, function(k) {
        mb <- factor(as.integer(rhdf5::h5read(stab_path, paste0(cl_mb_prefix, "/", k))))
        ecc <- as.numeric(rhdf5::h5read(
            stab_path,
            paste0(cl_ecc_prefix, sprintf("/%06d;%s", as.integer(k), stable_clustering_method))
        ))
        ecc_order <- as.numeric(rhdf5::h5read(
            stab_path,
            paste0(cl_ecc_prefix, sprintf("_order/%06d;%s", as.integer(k), stable_clustering_method))
        ))

        temp_df <- data.frame(mb, ecc[ecc_order])
        colnames(temp_df) <- c(paste0("stable_", k, "_clusters"), paste0("ecc_", k))
        return(temp_df)
    })

    return(do.call(cbind, df_list))
}

get_metadata_from_clustassess_app <- function(app_folder) {
    metadata_df <- readRDS(file.path(app_folder, "metadata.rds"))$metadata
    metadata_values <- colnames(metadata_df)
    metadata_values[which(metadata_values == "sample_name")] <- "sample_names"
    colnames(metadata_df) <- metadata_values

    return(metadata_df)
}

get_expression_matrix_from_clustassess_app <- function(app_folder,
                                                       stable_feature_type,
                                                       stable_feature_set_size,
                                                       use_all_genes = FALSE) {
    expr_path <- file.path(app_folder, "expression.h5")
    stab_path <- file.path(app_folder, "stability.h5")
    if (!file.exists(expr_path)) {
        stop(paste0("The expression.h5 file does not exist in the provided app_folder: ", app_folder))
    }

    available_configs <- rhdf5::h5read(stab_path, "feature_ordering/stable")
    available_ftypes <- names(available_configs)

    if (!stable_feature_type %in% available_ftypes) {
        stop(paste0("The provided stable_feature_type: ", stable_feature_type, " is not available in the app.\nAvailable options: ", paste(available_ftypes, collapse = ", ")))
    }

    available_fsizes <- available_configs[[stable_feature_type]]

    if (!stable_feature_set_size %in% available_fsizes) {
        stop(paste0("The provided stable_feature_set_size: ", stable_feature_set_size, " is not available in the app.\nAvailable options: ", paste(available_fsizes, collapse = ", ")))
    }

    gene_names <- as.character(rhdf5::h5read(expr_path, "genes"))

    if (!use_all_genes) {
        used_genes <- as.character(rhdf5::h5read(stab_path, paste0(stable_feature_type, "/feature_list")))
        used_genes <- used_genes[seq_len(as.integer(stable_feature_set_size))]
        expr_matrix <- rhdf5::h5read(expr_path, "expression_matrix")
        rownames(expr_matrix) <- gene_names
        expr_matrix <- expr_matrix[used_genes, , drop = FALSE]
        gc()
        gene_names <- used_genes
    } else {
        expr_matrix <- rhdf5::h5read(expr_path, "expression_matrix")
    }

    rownames(expr_matrix) <- gene_names
    colnames(expr_matrix) <- as.character(rhdf5::h5read(expr_path, "cells"))

    return(expr_matrix)
}

get_info_from_clustassess_app <- function(app_folder,
                                          info_type = c("expression_matrix", "metadata", "dim_reduction", "clusters"),
                                          stable_feature_type,
                                          stable_feature_set_size,
                                          stable_clustering_method,
                                          stable_n_clusters = NULL,
                                          dim_type = c("pca", "umap"),
                                          dim_key = "PC_",
                                          use_all_genes = FALSE) {
    info_type <- info_type[1]
    if (!(info_type %in% c("expression_matrix", "metadata", "dim_reduction", "clusters"))) {
        stop(paste0("The provided info_type: ", info_type, " is not available."))
    }

    if (info_type == "expression_matrix") {
        return(get_expression_matrix_from_clustassess_app(
            app_folder,
            stable_feature_type,
            stable_feature_set_size,
            use_all_genes
        ))
    }

    if (info_type == "dim_reduction") {
        return(get_dim_reduction_from_clustassess_app(
            app_folder,
            dim_type,
            dim_key,
            stable_feature_type,
            stable_feature_set_size
        ))
    }

    clusters_df <- get_clusters_from_clustassess_app(
        app_folder,
        stable_feature_type,
        stable_feature_set_size,
        stable_clustering_method
    )

    if (info_type == "clusters") {
        return(clusters_df)
    }

    metadata_df <- get_metadata_from_clustassess_app(app_folder)

    if (is.null(stable_n_clusters)) {
        return(cbind(metadata_df, clusters_df))
    }

    for (k in stable_n_clusters) {
        if (!(paste0("stable_", k, "_clusters") %in% colnames(metadata_df))) {
            metadata_df[[paste0("stable_", k, "_clusters")]] <- clusters_df[[paste0("stable_", k, "_clusters")]]
            metadata_df[[paste0("ecc_", k)]] <- clusters_df[[paste0("ecc_", k)]]
        }
    }

    return(metadata_df)
}

#' Create monocle object
#'
#' @description Use a normalized expression matrix and, potentially, an already
#' generated PCA / UMAP embedding, to create a Monocle object.
#'
#' @param normalized_expression_matrix The normalized expression matrix
#' having genes on rows and cells on columns.
#' @param count_matrix The count matrix having genes on rows and cells on
#' columns. If NULL, the normalized_expression_matrix will be used.
#' @param pca_embedding The PCA embedding of the expression matrix. If NULL, the
#' pca will be created using the `monocle3` package (default parameters).
#' @param umap_embedding The UMAP embedding of the expression matrix. If NULL, the
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
#' # uncomment to create the monocle object
#' mon_obj <- create_monocle_default(
#'     normalized_expression_matrix = expr_matrix,
#'     pca_emb = NULL,
#'     umap_emb = NULL,
#'     metadata_df = NULL
#' )
#' }
create_monocle_default <- function(normalized_expression_matrix,
                                   count_matrix = NULL,
                                   pca_embedding = NULL,
                                   umap_embedding = NULL,
                                   metadata_df = NULL) {
    if (is.null(count_matrix)) {
        count_matrix <- normalized_expression_matrix
    } else {
        count_matrix <- count_matrix[rownames(normalized_expression_matrix), colnames(normalized_expression_matrix)]
    }

    # check if normalized expression matrix is sparse
    if (!inherits(normalized_expression_matrix, "dgCMatrix")) {
        normalized_expression_matrix <- methods::as(normalized_expression_matrix, "dgCMatrix")
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
        num_dim = ifelse(is.null(pca_embedding), 30, ncol(pca_embedding)),
        norm_method = "none",
        scaling = FALSE
    )

    monocle_cds@assays@data$normalized_data <- normalized_expression_matrix

    if (!is.null(pca_embedding)) {
        rownames(pca_embedding) <- cell_names
        colnames(pca_embedding) <- paste0("PC_", seq_len(ncol(pca_embedding)))
        monocle_cds@int_colData$reducedDims$PCA <- pca_embedding
    }

    if (!is.null(umap_embedding)) {
        rownames(umap_embedding) <- cell_names
        colnames(umap_embedding) <- paste0("UMAP_", seq_len(ncol(umap_embedding)))
        monocle_cds@int_colData$reducedDims$UMAP <- umap_embedding
    } else {
        monocle_cds <- monocle3::reduce_dimension(
            cds = monocle_cds,
            reduction_method = "UMAP",
            umap.n_neighbors = 25,
            umap.min_dist = 0.3,
            umap.metric = "cosine"
        )
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
#' }
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
    gene_names <- as.character(rhdf5::h5read(expr_path, "genes"))
    cell_names <- as.character(rhdf5::h5read(expr_path, "cells"))

    # expr_matrix <- Matrix::Matrix(
    #     rhdf5::h5read(expr_path, "expression_matrix"),
    #     sparse = TRUE,
    #     dimnames = list(gene_names, cell_names)
    # )
    expr_matrix <- rhdf5::h5read(expr_path, "expression_matrix")
    rownames(expr_matrix) <- gene_names
    colnames(expr_matrix) <- cell_names

    available_configs <- rhdf5::h5read(stab_path, "feature_ordering/stable")
    available_ftypes <- names(available_configs)

    if (!stable_feature_type %in% available_ftypes) {
        stop(paste0("The provided stable_feature_type: ", stable_feature_type, " is not available in the app.\nAvailable options: ", paste(available_ftypes, collapse = ", ")))
    }

    available_fsizes <- available_configs[[stable_feature_type]]

    if (!stable_feature_set_size %in% available_fsizes) {
        stop(paste0("The provided stable_feature_set_size: ", stable_feature_set_size, " is not available in the app.\nAvailable options: ", paste(available_fsizes, collapse = ", ")))
    }

    if (!use_all_genes) {
        used_genes <- rhdf5::h5read(stab_path, paste0(stable_feature_type, "/feature_list"))
        used_genes <- used_genes[seq_len(as.integer(stable_feature_set_size))]
        expr_matrix <- expr_matrix[used_genes, , drop = FALSE]
        gc()
        gene_names <- used_genes
    }

    prefix <- paste0(stable_feature_type, "/", stable_feature_set_size)
    cl_mb_prefix <- paste0(prefix, "/clustering_stability/split_by_k/mbs/", stable_clustering_method)
    cl_ecc_prefix <- paste0(prefix, "/clustering_stability/split_by_k/ecc")
    metadata_df <- readRDS(file.path(app_folder, "metadata.rds"))$metadata

    tryCatch({
        available_n_clusters <- names(rhdf5::h5read(stab_path, cl_mb_prefix))

    }, error = function(e) {
        stop(paste0("The provided stable_clustering_method: ", stable_clustering_method, " is not available in the app."))
    })

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
        metadata_df[[paste0("stable_", k, "_clusters")]] <- factor(as.integer(rhdf5::h5read(stab_path, paste0(cl_mb_prefix, "/", k))))
        ecc <- as.numeric(rhdf5::h5read(
            stab_path,
            paste0(cl_ecc_prefix, sprintf("/%06d;%s", as.integer(k), stable_clustering_method))
        ))
        ecc_order <- as.numeric(rhdf5::h5read(
            stab_path,
            paste0(cl_ecc_prefix, sprintf("_order/%06d;%s", as.integer(k), stable_clustering_method))
        ))
        metadata_df[[paste0("ecc_", k)]] <- ecc[ecc_order]
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

#' Create Seurat object
#'
#' @description Use a normalized expression matrix and, potentially, an already
#' generated PCA / UMAP embedding, to create a Seurat object.
#'
#' @param normalized_expression_matrix The normalized expression matrix
#' having genes on rows and cells on columns.
#' @param count_matrix The count matrix having genes on rows and cells on
#' columns. If NULL, the normalized_expression_matrix will be used.
#' @param pca_embedding The PCA embedding of the expression matrix. If NULL, the
#' pca will be created using the `Seurat` package (default parameters).
#' @param umap_embedding The UMAP embedding of the expression matrix. If NULL, the
#' umap will be created using the `Seurat` package (default parameters).
#' @param metadata_df The metadata dataframe having the cell names as rownames.
#' If NULL, a dataframe with a single column named `identical_ident` will be
#' created.
#'
#' @return A Seurat object of the expression matrix, having the stable number
#' of clusters identified by ClustAssess.
#'
#' @export
create_seurat_object_default <- function(normalized_expression_matrix,
                                         count_matrix = NULL,
                                         pca_embedding = NULL,
                                         umap_embedding = NULL,
                                         metadata_df = NULL) {


    if (is.null(count_matrix)) {
        count_matrix <- normalized_expression_matrix
    } else {
        count_matrix <- count_matrix[rownames(normalized_expression_matrix), colnames(normalized_expression_matrix)]
    }

    if (!inherits(normalized_expression_matrix, "dgCMatrix")) {
        normalized_expression_matrix <- Matrix::Matrix(normalized_expression_matrix, sparse = TRUE)
    }

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

    seurat_obj <- Seurat::CreateSeuratObject(
        counts = count_matrix,
        meta.data = metadata_df
    )

    seurat_obj <- Seurat::NormalizeData(seurat_obj)
    seurat_obj@assays$RNA@layers$data <- normalized_expression_matrix
    seurat_obj <- Seurat::ScaleData(seurat_obj)
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- Seurat::RunPCA(seurat_obj)

    if (!is.null(pca_embedding)) {
        rownames(pca_embedding) <- cell_names
        colnames(pca_embedding) <- paste0("PC_", seq_len(ncol(pca_embedding)))
        seurat_obj$pca <- SeuratObject::CreateDimReducObject(
            embedding = pca_embedding,
            key = "PCA"
        )
    }

    seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)

    if (!is.null(umap_embedding)) {
        rownames(umap_embedding) <- cell_names
        colnames(umap_embedding) <- paste0("UMAP_", seq_len(ncol(umap_embedding)))
        seurat_obj$umap <- SeuratObject::CreateDimReducObject(
            embedding = umap_embedding,
            key = "UMAP"
        )
    }

    return(seurat_obj)
}

#' Create Seurat object from a ClustAssess shiny app
#'
#' @description Use the files generated in the ClustAssess app to create a
#' Seurat object which has the stable number of clusters.
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
#' @return A Seurat object of the expression matrix, having the stable number
#' of clusters identified by ClustAssess.
#'
#' @export
create_seurat_object_from_clustassess_app <- function(app_folder,
                                                      stable_feature_type,
                                                      stable_feature_set_size,
                                                      stable_clustering_method,
                                                      stable_n_clusters = NULL,
                                                      use_all_genes = FALSE) {
    create_seurat_object_default(
        normalized_expression_matrix = get_info_from_clustassess_app(
            app_folder,
            info_type = "expression_matrix",
            stable_feature_type = stable_feature_type,
            stable_feature_set_size = stable_feature_set_size,
            stable_clustering_method = stable_clustering_method,
            stable_n_clusters = stable_n_clusters,
            use_all_genes = use_all_genes
        ),
        pca_embedding = get_info_from_clustassess_app(
            app_folder,
            info_type = "dim_reduction",
            dim_type = "pca",
            dim_key = "PC_",
            stable_feature_type = stable_feature_type,
            stable_feature_set_size = stable_feature_set_size,
            stable_clustering_method = stable_clustering_method,
            stable_n_clusters = stable_n_clusters
        ),
        umap_embedding = get_info_from_clustassess_app(
            app_folder,
            info_type = "dim_reduction",
            dim_type = "umap",
            dim_key = "UMAP_",
            stable_feature_type = stable_feature_type,
            stable_feature_set_size = stable_feature_set_size,
            stable_clustering_method = stable_clustering_method,
            stable_n_clusters = stable_n_clusters
        ),
        metadata_df = get_info_from_clustassess_app(
            app_folder,
            info_type = "metadata",
            stable_feature_type = stable_feature_type,
            stable_feature_set_size = stable_feature_set_size,
            stable_clustering_method = stable_clustering_method,
            stable_n_clusters = stable_n_clusters
        )
    )
}
