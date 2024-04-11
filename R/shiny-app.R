# BUG the app doesn't work fully on Windows - check the problems
ppi <- 72
single_color <- "#025147"
generate_colours <- function(n_unique_values, qualpalr_colorspace, single_color = "#017c6b") {
    if (n_unique_values > 99) {
        return(sample(grDevices::colors(), n_unique_values))
    }

    if (n_unique_values > 1) {
        return(qualpalr::qualpal(n_unique_values, colorspace = qualpalr_colorspace)$hex)
    }

    single_color
}

#' Write the objects for the ClustAssess ShinyApp
#'
#' @description Given the output of the ClustAssess pipeline, the expression matrix
#' and the metadata, this function creates the files needed for the ClustAssess
#' ShinyApp. The files are written in the project_folder and are the following:
#' - metadata.rds: the metadata file
#' - stability.h5: contains the stability results
#' - expression.h5: contains the expression matrix and the rank matrix
#'
#' @param clustassess_object The output of the ClustAssess automatic pipeline
#' @param expression_matrix The expression matrix
#' @param metadata The metadata
#' @param project_folder The folder where the files will be written
#' @param compression_level The compression level for the h5 files (See `rhdf5::h5createFile`` for more details)
#' @param chunk_size The chunk size for the rank matrix (See `rhdf5::h5createDataset` for more details)
#' @param gene_variance_threshold The threshold for the gene variance; genes with variance below this threshold will be removed
#' @param summary_function The function used for summarizing the stability values; the default is `median`
#' @param qualpalr_colorspace The colorspace used for generating the colors; the default is `pretty`
#'
#' @return NULL (the files are written in the project_folder)
#'
#' @export
write_objects <- function(clustassess_object,
                          expression_matrix,
                          metadata,
                          project_folder = ".",
                          compression_level = 6,
                          chunk_size = 100,
                          gene_variance_threshold = 0,
                          summary_function = stats::median,
                          qualpalr_colorspace = "pretty") {
    metadata_file_name <- file.path(project_folder, "metadata.rds")
    stability_file_name <- file.path(project_folder, "stability.h5")
    expr_file_name <- file.path(project_folder, "expression.h5")

    if (file.exists(expr_file_name)) {
        stop(glue::glue("File {expr_file_name} already exists!"))
    }

    if (file.exists(stability_file_name)) {
        stop(glue::glue("File {stability_file_name} already exists!"))
    }

    if (file.exists(metadata_file_name)) {
        stop(glue::glue("File {metadata_file_name} already exists!"))
    }

    # metadata file
    print(glue::glue("[{Sys.time()}] Writing the metadata"))
    metadata$one_level <- factor(rep("one_level", nrow(metadata)))
    metadata_columns <- colnames(metadata)
    metadata_colors <- list()
    metadata_unique <- list()
    for (mtd_col in metadata_columns) {
        if (is.factor(metadata[, mtd_col])) {
            metadata[, mtd_col] <- droplevels(metadata[, mtd_col])
            metadata_unique[[mtd_col]] <- levels(metadata[, mtd_col])
            if (length(metadata_unique[[mtd_col]]) > 502) {
                next
            }


            metadata_colors[[mtd_col]] <- generate_colours(length(metadata_unique[[mtd_col]]), qualpalr_colorspace)
        } else if (is.character(metadata[, mtd_col])) {
            metadata[, mtd_col] <- factor(metadata[, mtd_col])
            metadata_unique[[mtd_col]] <- levels(metadata[, mtd_col])
            if (length(metadata_unique[[mtd_col]]) > 502) {
                next
            }

            metadata_colors[[mtd_col]] <- generate_colours(length(metadata_unique[[mtd_col]]), qualpalr_colorspace)
        } else if (is.logical(metadata[, mtd_col])) {
            metadata_unique[[mtd_col]] <- c(FALSE, TRUE)
            metadata_colors[[mtd_col]] <- qualpalr::qualpal(2, colorspace = qualpalr_colorspace)$hex
        }
    }

    saveRDS(
        list(
            metadata = metadata,
            metadata_colors = metadata_colors,
            metadata_unique = metadata_unique
        ),
        metadata_file_name
    )

    # establish the feature ordering (original and stable) and convert to data.table
    feature_ordering <- list(original = list(), stable = list(), original_incremental = list())
    resolution_values <- c()

    ftype_index <- 1
    nftypes <- length(clustassess_object$feature_stability$by_steps)
    unique_n_colors <- c(nftypes)

    clustassess_object$feature_stability$colours <- generate_colours(nftypes, qualpalr_colorspace)

    for (ftype in names(clustassess_object$feature_stability$by_steps)) {
        feature_ordering$original[[ftype]] <- names(clustassess_object$feature_stability$by_steps[[ftype]])
        fsize_index <- 1
        upper_limit <- length(clustassess_object$feature_stability$incremental[[ftype]])

        for (fsize in names(clustassess_object$feature_stability$by_steps[[ftype]])) {
            resolution_values <- union(resolution_values, names(clustassess_object$feature_stability$by_steps[[ftype]][[fsize]]))

            summary_values_by_step <- rep(0, length(clustassess_object$feature_stability$by_steps[[ftype]][[fsize_index]]))
            if (fsize_index <= upper_limit) {
                summary_values_incremental <- rep(0, length(clustassess_object$feature_stability$incremental[[ftype]][[fsize_index]]))
            }

            res_index <- 1

            for (resval in names(clustassess_object$feature_stability$by_steps[[ftype]][[fsize]])) {
                # mb <- original_stab_obj$feature_stability$by_steps[[ftype]][[fsize_index]][[resval]]$most_frequent_partition$mb # uncomment if you want to include the mb vectors aswell

                ecc_by_step <- clustassess_object$feature_stability$by_steps[[ftype]][[fsize_index]][[resval]]$ecc
                summary_values_by_step[res_index] <- summary_function(ecc_by_step)
                res_dt_by_step <- data.frame(ecc = ecc_by_step, res = resval) #  mb = mb if you want to also include the mb vector
                if (res_index == 1) {
                    size_dt_by_step <- res_dt_by_step
                } else {
                    size_dt_by_step <- rbind(size_dt_by_step, res_dt_by_step)
                }

                if (fsize_index <= upper_limit) {
                    ecc_incremental <- clustassess_object$feature_stability$incremental[[ftype]][[fsize_index]][[resval]]

                    summary_values_incremental[res_index] <- summary_function(ecc_incremental)
                    res_dt_incremental <- data.frame(ecc = ecc_incremental, res = resval) #  mb = mb if you want to also include the mb vector
                    if (res_index == 1) {
                        size_dt_incremental <- res_dt_incremental
                    } else {
                        size_dt_incremental <- rbind(size_dt_incremental, res_dt_incremental)
                    }
                }

                res_index <- res_index + 1
            }

            size_dt_by_step[, "fsize"] <- fsize_index
            summary_values_by_step <- data.frame(ecc = summary_values_by_step, fsize = fsize_index)

            if (fsize_index == 1) {
                ftype_dt_by_step <- size_dt_by_step
                ftype_summary_dt_by_step <- summary_values_by_step
            } else {
                ftype_dt_by_step <- rbind(ftype_dt_by_step, size_dt_by_step)
                ftype_summary_dt_by_step <- rbind(ftype_summary_dt_by_step, summary_values_by_step)
            }

            if (fsize_index > upper_limit) {
                next
            }

            size_dt_incremental[, "fsize"] <- fsize_index
            summary_values_incremental <- data.frame(ecc = summary_values_incremental, fsize = fsize_index)

            if (fsize_index == 1) {
                ftype_dt_incremental <- size_dt_incremental
                ftype_summary_dt_incremental <- summary_values_incremental
            } else {
                ftype_dt_incremental <- rbind(ftype_dt_incremental, size_dt_incremental)
                ftype_summary_dt_incremental <- rbind(ftype_summary_dt_incremental, summary_values_incremental)
            }

            fsize_index <- fsize_index + 1
        }

        ftype_dt_by_step[, "ftype"] <- ftype
        ftype_dt_incremental[, "ftype"] <- ftype

        ftype_summary_dt_by_step[, "ftype"] <- ftype
        ftype_summary_dt_incremental[, "ftype"] <- ftype
        if (ftype_index == 1) {
            overall_dtable_by_step <- ftype_dt_by_step
            overall_dtable_incremental <- ftype_dt_incremental
            overall_summary_dt_by_step <- ftype_summary_dt_by_step
            overall_summary_dt_incremental <- ftype_summary_dt_incremental
        } else {
            overall_dtable_by_step <- rbind(overall_dtable_by_step, ftype_dt_by_step)
            overall_dtable_incremental <- rbind(overall_dtable_incremental, ftype_dt_incremental)
            overall_summary_dt_by_step <- rbind(overall_summary_dt_by_step, ftype_summary_dt_by_step)
            overall_summary_dt_incremental <- rbind(overall_summary_dt_incremental, ftype_summary_dt_incremental)
        }

        ftype_index <- ftype_index + 1
    }

    feature_ordering$resolution <- stringr::str_sort(resolution_values, numeric = TRUE)
    # split the data tables by resolution
    clustassess_object$feature_stability$by_steps <- lapply(feature_ordering$resolution, function(resval) {
        subdt <- overall_dtable_by_step %>% dplyr::filter(.data$res == resval) # %>% dplyr::arrange(order(.data$ecc))
        subdt$fsize <- factor(subdt$fsize)
        subdt$resval <- NULL

        return(subdt)
    })
    names(clustassess_object$feature_stability$by_steps) <- feature_ordering$resolution

    clustassess_object$feature_stability$incremental <- lapply(feature_ordering$resolution, function(resval) {
        subdt <- overall_dtable_incremental %>%
            dplyr::filter(.data$res == resval) %>%
            dplyr::arrange(order(.data$ecc))
        subdt$fsize <- factor(subdt$fsize)

        return(subdt)
    })
    names(clustassess_object$feature_stability$incremental) <- feature_ordering$resolution

    clustassess_object$feature_stability$overall <- list(
        by_step = overall_summary_dt_by_step,
        incremental = overall_summary_dt_incremental
    )

    # store the names and ordering of the stable feature sizes
    for (ftype in names(feature_ordering$original)) {
        nsteps <- length(clustassess_object[[ftype]]) - 1
        feature_ordering$stable[[ftype]] <- names(clustassess_object[[ftype]])[seq_len(nsteps)]

        feature_ordering$original_incremental[[ftype]] <- rep(0, length(feature_ordering$original[[ftype]]) - 1)

        for (fsize_index in seq_len(length(feature_ordering$original[[ftype]]) - 1)) {
            feature_ordering$original_incremental[[ftype]][fsize_index] <- paste(
                feature_ordering$original[[ftype]][fsize_index],
                feature_ordering$original[[ftype]][fsize_index + 1],
                sep = "-"
            )
        }
    }

    # ---
    # create the clustering data tables
    # ecc_by_k_dt - columns: clustering method, k and ecc
    # ecc_by_res_dt - columns: clustering method, res and ecc
    for (ftype in names(feature_ordering$stable)) {
        for (fsize in feature_ordering$stable[[ftype]]) {
            # remove the adjacency matrix, not needed for the shiny app
            clustassess_object[[ftype]][[fsize]]$adj_matrix <- NULL
            mbs_list <- list()
            ecc_list <- list()
            ecc_list_by_res <- list()
            ecc_order_list <- list()
            ecc_order_list_by_res <- list()
            structure_list <- list()

            # remove the list of unecessary partitions
            cl_methods <- names(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k)
            first_dt_by_k <- TRUE
            first_dt_by_res <- TRUE
            # algorithm_names_mapping

            for (cl_method in cl_methods) {
                # mbs_list[[cl_method]] <- list()
                # ecc_list[[cl_method]] <- list()
                # ecc_order_list[[cl_method]] <- list()
                # ecc_order_list_by_res[[cl_method]] <- list()
                # ecc_list_by_res[[cl_method]] <- list()

                cl_method_index <- algorithm_names_mapping[[cl_method]]
                # --- split by resolution ---
                for (res in names(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]])) {
                    # update the ecc lists for the boxplots
                    ecc <- clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$ecc
                    ecc_order <- order(ecc)
                    res_format <- sprintf("%.6f", as.numeric(res))
                    ecc_list_by_res[[paste(res_format, cl_method, sep = ";")]] <- ecc[ecc_order]
                    ecc_order_list_by_res[[paste(res_format, cl_method, sep = ";")]] <- Matrix::invPerm(ecc_order)

                    # create the summary table for the complex plot
                    for (n_clusters in names(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters)) {
                        # clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$first_freq <- clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions[[1]]$freq
                        # clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$total_freq <- sum(
                        #     sapply(
                        #         clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions,
                        #         function(x) { x$freq }
                        #     )
                        # )
                        # clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions <- NULL

                        temp_by_res_summary <- data.frame(
                            k = as.integer(n_clusters),
                            res = as.numeric(res),
                            cl_method = cl_method,
                            first_freq = clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions[[1]]$freq,
                            total_freq = sum(
                                sapply(
                                    clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$partitions,
                                    function(x) {
                                        x$freq
                                    }
                                )
                            ),
                            ecc = summary_function(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution[[cl_method]][[res]]$clusters[[n_clusters]]$ecc)
                        )

                        if (first_dt_by_res) {
                            first_dt_by_res <- FALSE
                            by_res_summary <- temp_by_res_summary
                        } else {
                            by_res_summary <- rbind(by_res_summary, temp_by_res_summary)
                        }
                    }
                }

                # --- split by k ---
                structure_list[[cl_method]] <- names(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]])
                for (n_clusters in structure_list[[cl_method]]) {
                    unique_n_colors <- union(unique_n_colors, n_clusters)
                    mbs_list[[cl_method]][[n_clusters]] <- as.integer(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions[[1]]$mb)

                    # clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$mb <- as.integer(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions[[1]]$mb)
                    # clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$first_freq <- clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions[[1]]$freq
                    # clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$total_freq <- sum(
                    #     sapply(
                    #         clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions,
                    #         function(x) { x$freq }
                    #     )
                    # )
                    # clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$n_partitions <- length(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions)
                    # clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions <- NULL

                    # temp_n_cl_ecc_dt <- data.table::data.table(
                    #     ecc = clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$ecc,
                    #     k = n_clusters,
                    #     cl_method = cl_method_index
                    # )
                    ecc <- clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$ecc
                    ecc_order <- order(ecc)
                    n_clusters_format <- sprintf("%06d", as.integer(n_clusters))
                    ecc_list[[paste(n_clusters_format, cl_method, sep = ";")]] <- ecc[ecc_order]
                    ecc_order_list[[paste(n_clusters_format, cl_method, sep = ";")]] <- Matrix::invPerm(ecc_order)

                    temp_by_k_summary <- data.frame(
                        k = as.integer(n_clusters),
                        cl_method = cl_method,
                        n_partitions = length(clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions),
                        total_freq = sum(
                            sapply(
                                clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions,
                                function(x) {
                                    x$freq
                                }
                            )
                        ),
                        first_freq = clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[cl_method]][[n_clusters]]$partitions[[1]]$freq,
                        ecc = summary_function(ecc)
                    )


                    if (first_dt_by_k) {
                        first_dt_by_k <- FALSE
                        by_k_summary <- temp_by_k_summary
                    } else {
                        by_k_summary <- rbind(by_k_summary, temp_by_k_summary)
                    }
                }
            }

            clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_k <- list(
                mbs = mbs_list,
                summary = by_k_summary %>% dplyr::arrange(.data$k, .data$cl_method),
                ecc = ecc_list,
                ecc_order = ecc_order_list,
                structure_list = structure_list
            )

            clustassess_object[[ftype]][[fsize]]$clustering_stability$split_by_resolution <- list(
                summary = by_res_summary %>% dplyr::arrange(.data$res, .data$cl_method),
                ecc = ecc_list_by_res,
                ecc_order = ecc_order_list_by_res
            )
        }
    }
    unique_n_colors <- union(unique_n_colors, 1:4)
    unique_colors <- lapply(unique_n_colors, function(n) {
        n <- as.integer(n)
        if (n == 1) {
            return(single_color)
        }
        generate_colours(n, qualpalr_colorspace)
    })
    names(unique_colors) <- as.character(unique_n_colors)

    clustassess_object$feature_ordering <- feature_ordering
    feature_names <- names(feature_ordering$original)
    # genes_of_interest <- c()
    # for (feature_name in feature_names) {
    #     genes_of_interest <- union(genes_of_interest, clustassess_object[[feature_name]]$feature_list)
    # }

    # should sparse matrix be used for storing?
    # advantages: lower required space
    # disadvantages: would the rows be as easily accesible as in a normal matrix?
    print(glue::glue("[{Sys.time()}] Removing the genes from the expression matrix with low variance"))
    num_rows <- nrow(expression_matrix)
    num_chunks <- ceiling(num_rows / chunk_size)
    dense_chunks <- lapply(seq_len(num_chunks), function(i) {
        start_row <- (i - 1) * chunk_size + 1
        end_row <- min(i * chunk_size, num_rows)
        chunk_matrix <- as.matrix(expression_matrix[start_row:end_row, ])
        gene_vars <- matrixStats::rowVars(chunk_matrix)
        # filter by var threshold
        chunk_matrix[gene_vars >= gene_variance_threshold, ]
    })

    print("merging chunks")
    expression_matrix <- do.call(rbind, dense_chunks)
    print("chunks merged")

    # others_sorted <- sort(rownames(expression_matrix)[!(rownames(expression_matrix) %in% genes_of_interest)])
    # genes <- rownames(expression_matrix) #c(genes_of_interest, others_sorted)
    # cells <- colnames(expression_matrix)

    gene_avg_expression <- matrixStats::rowMeans2(expression_matrix)
    order_expression <- order(gene_avg_expression, decreasing = TRUE)
    # order the genes based on their average expression
    expression_matrix <- expression_matrix[order_expression, ]
    gene_avg_expression <- gene_avg_expression[order_expression]

    print(glue::glue("[{Sys.time()}] Writing the gene matrix"))

    rhdf5::h5createFile(expr_file_name)

    # rhdf5::h5createDataset(expr_file_name, "matrix_of_interest",
    #     level = compression_level,
    #     dims = c(length(genes_of_interest), ncol(expression_matrix)),
    #     storage.mode = "double",
    #     # chunk = c(length(genes_of_interest), 1)
    #     chunk = c(1, ncol(expression_matrix))
    # )

    rhdf5::h5createDataset(expr_file_name, "rank_matrix",
        level = compression_level,
        dims = c(nrow(expression_matrix), ncol(expression_matrix)),
        storage.mode = "integer",
        chunk = c(chunk_size, ncol(expression_matrix))
    )

    rhdf5::h5createDataset(expr_file_name, "expression_matrix",
        level = compression_level,
        dims = c(nrow(expression_matrix), ncol(expression_matrix)),
        maxdims = c(nrow(expression_matrix), ncol(expression_matrix)),
        storage.mode = "double",
        chunk = c(1, ncol(expression_matrix))
    )

    # rhdf5::h5write(genes_of_interest, expr_file_name, "genes_of_interest")
    # rhdf5::h5write(others_sorted, expr_file_name, "genes_others")
    rhdf5::h5write(rownames(expression_matrix), expr_file_name, "genes")
    rhdf5::h5write(colnames(expression_matrix), expr_file_name, "cells")
    rhdf5::h5write(gene_avg_expression, expr_file_name, "average_expression")
    rhdf5::h5write(chunk_size, expr_file_name, "chunk_size")
    rhdf5::h5write(expression_matrix, expr_file_name, "expression_matrix")
    # rhdf5::h5write(expression_matrix[genes_of_interest, ], expr_file_name, "matrix_of_interest")
    # rhdf5::h5write(expression_matrix[others_sorted, ], expr_file_name, "matrix_others")

    print(glue::glue("[{Sys.time()}] Writing the rank matrix"))
    nchunks <- ceiling(nrow(expression_matrix) / chunk_size)
    for (i in seq_len(nchunks - 1)) {
        index_list <- seq(from = chunk_size * (i - 1) + 1, by = 1, length.out = chunk_size)
        rhdf5::h5write(
            t(apply(
                expression_matrix[index_list, ],
                1,
                function(x) {
                    rank(x, ties.method = "min")
                }
            )),
            expr_file_name,
            "rank_matrix",
            index = list(index_list, NULL)
        )
    }

    index_list <- seq(from = chunk_size * (nchunks - 1) + 1, by = 1, to = nrow(expression_matrix))
    if (length(index_list) == 1) {
        rhdf5::h5write(
            rank(expression_matrix[index_list, ], ties.method = "min"),
            expr_file_name,
            "rank_matrix",
            index = list(index_list, NULL)
        )
    } else {
        rhdf5::h5write(
            t(apply(
                expression_matrix[index_list, ],
                1,
                function(x) {
                    rank(x, ties.method = "min")
                }
            )),
            expr_file_name,
            "rank_matrix",
            index = list(index_list, NULL)
        )
    }

    # rank_matrix <- matrix(nrow = length(genes_of_interest), ncol = ncol(expression_matrix))

    # for (i in seq_len(nrow(rank_matrix))) {
    #   rank_matrix[i, ] <- rank(expression_matrix[genes_of_interest[i], ], ties.method = "min")
    # }
    # rhdf5::h5write(rank_matrix, expr_file_name, "rank_of_interest")

    print(glue::glue("[{Sys.time()}] Writing the feature stability"))

    # the object for stability
    rhdf5::h5createFile(stability_file_name)
    rhdf5::h5write(
        clustassess_object$feature_stability,
        stability_file_name,
        "feature_stability",
        level = compression_level
    )
    rhdf5::h5write(clustassess_object$feature_ordering, stability_file_name, "feature_ordering", level = compression_level)
    print(glue::glue("[{Sys.time()}] Writing the stability object"))

    for (feature_type in names(clustassess_object$feature_ordering[[1]])) {
        rhdf5::h5write(clustassess_object[[feature_type]], stability_file_name, feature_type, level = compression_level)
    }
    rhdf5::h5write(unique_colors, stability_file_name, "colors", level = compression_level)
    print(glue::glue("[{Sys.time()}] Done"))

    rhdf5::h5closeAll()
}

#' @rdname write_shiny_app
#' @export
write_shiny_app.Seurat <- function(object,
                                   metadata = NULL,
                                   assay_name,
                                   clustassess_object,
                                   project_folder,
                                   compression_level = 6,
                                   summary_function = stats::median,
                                   shiny_app_title = "",
                                   organism_enrichment = "hsapiens",
                                   height_ratio = 0.6,
                                   qualpalr_colorspace = "pretty") {
    write_shiny_app(
        object = object@assays[[assay_name]]@data,
        metadata = object@meta.data,
        clustassess_object = clustassess_object,
        project_folder = project_folder,
        compression_level = compression_level,
        summary_function = summary_function,
        shiny_app_title = shiny_app_title,
        organism_enrichment = organism_enrichment,
        height_ratio = height_ratio,
        qualpalr_colorspace = qualpalr_colorspace
    )
}


#' @rdname write_shiny_app
#' @export
write_shiny_app.default <- function(object,
                                    metadata,
                                    assay_name = NULL,
                                    clustassess_object,
                                    project_folder,
                                    compression_level = 6,
                                    summary_function = stats::median,
                                    shiny_app_title = "",
                                    organism_enrichment = "hsapiens",
                                    height_ratio = 0.6,
                                    qualpalr_colorspace = "pretty") {
    nFeature <- DelayedMatrixStats::colSums2(object > 0)
    warning_message <- ""
    shiny_app_title <- paste("ClustAssess ShinyApp -", shiny_app_title)

    for (ftype in names(clustassess_object$feature_stability$by_steps)) {
        fsizes <- as.integer(names(clustassess_object$feature_stability$by_steps[[ftype]]))
        fsizes <- fsizes[which(fsizes > stats::median(nFeature))]
        if (length(fsizes) > 0) {
            warning_message <- paste0(warning_message, ftype, " - ", paste(fsizes, collapse = ", "), "; ")
        }
    }

    if (stringr::str_length(warning_message) > 0) {
        while (TRUE) {
            warning(glue::glue("WARNING: The following configurations -- {warning_message} have a size above the average number of nFeatures per cell - {median(nFeature)}.\nIncreasing the number of features above the average will lead to introduction of noise in the data, thus we recommed re-running ClustAssess with lower values for the feature steps.\nPlease type `yes` or `no` if you want to continue creating the ClustAssess shiny app."), immediate. = TRUE)
            user_input <- readline()
            if (user_input == "no") {
                return()
            }

            if (user_input == "yes") {
                break
            }
        }
    }

    if (!dir.exists(project_folder)) {
        dir.create(project_folder)
    }

    write_objects(
        clustassess_object = clustassess_object,
        expression_matrix = object,
        metadata = metadata,
        project_folder = project_folder,
        compression_level = compression_level,
        summary_function = summary_function,
        qualpalr_colorspace = qualpalr_colorspace
    )

    file_content <- paste0("
        library(ggplot2)
        library(shiny)
        library(shinyjs)
        library(rhdf5)
        library(ClustAssess)
        library(dplyr)

        tabs_numbers <- seq_len(6)
        names(tabs_numbers) <- c(
            \"Home\",
            \"Dimensionality Reduction\",
            \"Graph Construction\",
            \"Graph Clustering\",
            \"Comparison\",
            \"Sandbox\"
        )

        ui <- fluidPage(
        tags$head(
            tags$style(
            HTML(
                \".first-element-tab {
                margin-top: 150px;
                }\",
                \".page-info {
                top: 0px;
                margin-top: 55px;
                position: fixed;
                z-index: 100;
                }\",
                \"#tabset_id {
                            position: fixed;
                            width: 100%;
                            background-color: white;
                            top: 0;
                            z-index: 100;
                            font-size: 20px;
                            margin-left: 40px;
                            margin-top: 50px;
                            }\",
                \".show_config {
                            z-index: 10000;
                            position: fixed;
                            top: 10px;
                            font-size: 15px;
                            right: 25px;
                        }\",
                \".shiny-split-layout > div {
                            overflow: visible;
                            white-space: normal;
                        }\",
                \"div.vertical-line{
                        width: 1px; /* Line width */
                        background-color: black; /* Line color */
                        height: 100%; /* Override in-line if you want specific height. */
                        float: left; /* Causes the line to float to left of content.
                            You can instead use position:absolute or display:inline-block
                            if this fits better with your design */
                        }\"
            )
            )
        ),
        tags$head(tags$script(HTML('
                    var dimension = [0, 0];
                    var resizeId;
                    $(document).on(\"shiny:connected\", function(e) {
                        dimension[0] = Math.max(window.innerWidth - 20, 400);
                        dimension[1] = Math.max(window.innerHeight - 30, 400);
                        Shiny.onInputChange(\"dimension\", dimension);
                    });

                    function transferWindowSize() {
                        console.log(dimension);
                        console.log(window.innerHeight);

                        let dif_width = Math.abs(Math.max(window.innerWidth - 20, 400) - dimension[0]);
                        let dif_height = Math.abs(Math.max(window.innerHeight - 30, 400) - dimension[1]);
                        console.log(dif_height);

                        if (dif_width >= 200 || dif_height >= 200) {
                        console.log(\"Changed\")
                        dimension[0] = Math.max(window.innerWidth - 20, 400);
                        dimension[1] = Math.max(window.innerHeight - 30, 400);
                        Shiny.onInputChange(\"dimension\", dimension);
                        }
                    }

                    $(window).resize(function() {
                        clearTimeout(resizeId);
                        resizeId = setTimeout(transferWindowSize, 500);
                    });
                    '))),
        tags$head(
            tags$link(rel = \"stylesheet\", type = \"text/css\", href = \"www/shiny.css\")
        ),
        useShinyjs(),
        navbarPage(\"", shiny_app_title, "\", windowTitle = \"", shiny_app_title, "\", position = \"fixed-top\", inverse = TRUE),
        tabsetPanel(
            id = \"tabset_id\",
            ui_landing_page(\"landing_page\"),
            ui_dimensionality_reduction(\"dim_reduc\"),
            ui_graph_construction(\"graph_constr\"),
            ui_graph_clustering(\"graph_clust\"),
            ui_comparisons(\"comparison\"),
            ui_sandbox(\"sandbox\")
        )
        )

        server <- function(input, output, session) {
        hideTab(\"tabset_id\", \"Dimensionality Reduction\")
        hideTab(\"tabset_id\", \"Graph Construction\")
        hideTab(\"tabset_id\", \"Graph Clustering\")
        hideTab(\"tabset_id\", \"Comparison\")
        hideTab(\"tabset_id\", \"Sandbox\")
        fchoice <- reactiveVal(
            list(
            chosen_feature_type = \"Highly_Variable\",
            chosen_set_size = \"1500\"
            )
        )
        cchoice <- reactiveVal()
        height_ratio <- ", height_ratio, " # used to control the height of the plot proportional to the window's height
        observe({
            runjs(\"window.scrollTo(0, 0)\")
            tab_number <- as.integer(tabs_numbers[input$tabset_id])
            print(tab_number)

            if (tab_number == 1) {
            server_landing_page(\"landing_page\", height_ratio, reactive(input$dimension), session, \"", organism_enrichment, "\")
            }

            if (tab_number == 2) {
            fchoice(server_dimensionality_reduction(\"dim_reduc\", session))
            }

            if (tab_number == 3) {
            server_graph_construction(\"graph_constr\", fchoice())
            }

            if (tab_number == 4) {
            cchoice(server_graph_clustering(\"graph_clust\", fchoice(), session))
            }

            if (tab_number == 5) {
            server_comparisons(\"comparison\", fchoice(), cchoice())
            }

            if (tab_number == 6) {
            server_sandbox(\"sandbox\")
            }
        }) %>% bindEvent(input$tabset_id)
        }

        shinyApp(ui, server)
    ")


    write(file_content, file.path(project_folder, "app.R"))
    styler::style_file(file.path(project_folder, "app.R"))
}

#' Add metadata to ClustAssess ShinyApp
#'
#' @description Adds new metadata into the ClustAssess ShinyApp without having
#' to update the object and re-create the app.
#'
#' @param app_folder The folder containing the ClustAssess ShinyApp
#' @param metadata The new metadata to be added
#' @param qualpalr_colorspace The colorspace to be used for the metadata
#'
#' @return NULL - the metadata object is updated in the app folder
#'
#' @export
add_metadata <- function(app_folder,
                         metadata,
                         qualpalr_colorspace = "pretty") {
    metadata_file_name <- file.path(app_folder, "metadata.rds")

    if (!file.exists(metadata_file_name)) {
        quit(glue::glue("The app folder {app_folder} should contain a file named `metadata.rds`!"))
    }

    if (!is.data.frame(app_folder)) {
        quit("The `metadata` object should be a dataframe!")
    }

    existing_metadata <- readRDS(metadata_file_name)
    metadata_columns <- colnames(metadata)
    for (mtd_col in metadata_columns) {
        if (is.factor(metadata[, mtd_col])) {
            metadata[, mtd_col] <- droplevels(metadata[, mtd_col])
            existing_metadata$metadata_unique[[mtd_col]] <- levels(metadata[, mtd_col])
            if (length(existing_metadata$metadata_unique[[mtd_col]]) > 502) {
                next
            }


            existing_metadata$metadata_colors[[mtd_col]] <- generate_colours(length(existing_metadata$metadata_unique[[mtd_col]]), qualpalr_colorspace)
        } else if (is.character(metadata[, mtd_col])) {
            metadata[, mtd_col] <- factor(metadata[, mtd_col])
            existing_metadata$metadata_unique[[mtd_col]] <- levels(metadata[, mtd_col])
            if (length(existing_metadata$metadata_unique[[mtd_col]]) > 502) {
                next
            }

            existing_metadata$metadata_colors[[mtd_col]] <- generate_colours(length(existing_metadata$metadata_unique[[mtd_col]]), qualpalr_colorspace)
        } else if (is.logical(metadata[, mtd_col])) {
            existing_metadata$metadata_unique[[mtd_col]] <- c(FALSE, TRUE)
            existing_metadata$metadata_colors[[mtd_col]] <- qualpalr::qualpal(2, colorspace = qualpalr_colorspace)$hex
        }

        existing_metadata$metadata[[mtd_col]] <- metadata[, mtd_col]
    }

    saveRDS(existing_metadata, metadata_file_name)
}
