# devtools::load_all()

# set.seed(2024)
# # create an artificial PCA embedding
# pca_embedding <- matrix(runif(100 * 30), nrow = 100)
# rownames(pca_embedding) <- paste0("cell_", seq_len(nrow(pca_embedding)))
# colnames(pca_embedding) <- paste0("PC_", 1:30)


# adj_matrix <- getNNmatrix(
#     RANN::nn2(pca_embedding, k = 10)$nn.idx,
#     10,
#     0,
#     -1
# )$nn
# rownames(adj_matrix) <- paste0("cell_", seq_len(nrow(adj_matrix)))
# colnames(adj_matrix) <- paste0("cell_", seq_len(ncol(adj_matrix)))

# alternatively, the adj_matrix can be calculated
# using the `Seurat::FindNeighbors` function.

clust_diff_obj <- assess_clustering_stability(
    graph_adjacency_matrix = adj_matrix,
    resolution = seq(from = 0.1, to = 2, by = 0.1),
    n_repetitions = 10,
    clustering_algorithm = 1:2,
    verbose = TRUE
)

# plot_k_n_partitions(clust_object = clust_diff_obj)


clust_hierplot_get_y_mapping <- function(df) {
    cl_names <- colnames(df)
    y_mapping <- list()
    n_max_cl <- 0

    for (i in cl_names) {
        # unique_cl <- stringr::str_sort(unique(df[, i]), numeric = TRUE)
        if (is.factor(df[, i])) {
            unique_cl <- levels(df[, i])
        } else {
            unique_cl <- unique(df[, i])
        }
        y_mapping[[i]] <- lapply(unique_cl, I)
        names(y_mapping[[i]]) <- unique_cl

        n_max_cl <- max(n_max_cl, length(unique_cl))
    }

    for (i in cl_names) {
        y_fact <- n_max_cl / length(y_mapping[[i]])
        for (j in seq_along(y_mapping[[i]])) {
            y_mapping[[i]][[j]] <- y_fact * j - y_fact / 2
        }
    }

    return(y_mapping)
}

clust_hierplot_create_node_df <- function(df, consistency_list, y_mapping) {
    nodes_df <- NA
    for (i in seq(from = 1, to = length(y_mapping))) {
        cl <- df[, i]
        unique_cl <- names(y_mapping[[names(y_mapping)[i]]])

        for (cl_name in unique_cl) {
            temp_df <- data.frame(
                part_name = names(y_mapping)[i],
                clust_name = cl_name,
                ecc = mean(consistency_list[[names(y_mapping)[i]]][cl == cl_name]),
                cluster_size = sum(cl == cl_name),
                actual_y = y_mapping[[names(y_mapping)[i]]][[cl_name]]
            )

            if (!is.data.frame(nodes_df) && is.na(nodes_df)) {
                nodes_df <- temp_df
            } else {
                nodes_df <- rbind(nodes_df, temp_df)
            }
        }
    }
    nodes_df$part_name <- factor(nodes_df$part_name, levels = names(y_mapping))
    nodes_df$clust_name <- as.character(nodes_df$clust_name)

    return(nodes_df)
}

clust_hierplot_create_edge_df <- function(df, y_mapping) {
    edges_df <- NA
    for (i in seq(from = 1, to = length(y_mapping) - 1)) {
        cl1 <- df[, i]
        cl2 <- df[, i + 1]

        ctg_table <- table(cl1, cl2)
        cl_sizes1 <- rowSums(ctg_table)
        cl_sizes2 <- colSums(ctg_table)

        for (j in rownames(ctg_table)) {
            for (k in colnames(ctg_table)) {
                ecs_val <- ctg_table[j, k] / max(cl_sizes1[j], cl_sizes2[k])

                if (ctg_table[j, k] == 0) {
                    next
                }

                temp_df <- data.frame(
                    part_name_from = names(y_mapping)[i],
                    part_name_to = names(y_mapping)[i + 1],
                    clust_name_from = y_mapping[[names(y_mapping)[i]]][[j]],
                    clust_name_to = y_mapping[[names(y_mapping)[i + 1]]][[k]],
                    ecs = ecs_val,
                    intersect_size = ctg_table[j, k]
                )

                if (!is.data.frame(edges_df) && is.na(edges_df)) {
                    edges_df <- temp_df
                } else {
                    edges_df <- rbind(edges_df, temp_df)
                }
            }
        }
    }
    edges_df$part_name_from <- factor(edges_df$part_name_from, levels = names(y_mapping))
    edges_df$part_name_to <- factor(edges_df$part_name_to, levels = names(y_mapping))

    return(edges_df %>% dplyr::arrange(.data$intersect_size))
}

clust_hierplot_create_dfs <- function(df, consistency_list, edge_threshold = 0.3) {
    y_mapping <- clust_hierplot_get_y_mapping(df)
    node_df <- clust_hierplot_create_node_df(df, consistency_list, y_mapping)
    edge_df <- clust_hierplot_create_edge_df(df, y_mapping)
    quantile_thresh <- stats::quantile(edge_df$intersect_size, edge_threshold)
    edge_df <- edge_df %>% filter(.data$intersect_size >= quantile_thresh)
    return(list(
        node_df,
        edge_df
    ))
}

clust_hierarchical_plot_default <- function(clustering_df,
                                            consistency_list,
                                            edge_threshold = 0.3,
                                            range_point_size = c(1, 6),
                                            range_edge_width = c(0.01, 3),
                                            edge_palette_name = "RColorBrewer::Greys",
                                            edge_palette_inverse = FALSE,
                                            node_palette_name = "viridis::rocket",
                                            node_palette_inverse = TRUE) {
    node_edge_list <- clust_hierplot_create_dfs(clustering_df, consistency_list, edge_threshold)
    edge_pal_colours <- get_colour_vector_from_palette(edge_palette_name, edge_palette_inverse, "RColorBrewer::Greys")
    node_pal_colours <- get_colour_vector_from_palette(node_palette_name, node_palette_inverse, "viridis::viridis")
    
    gplot_obj <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = node_edge_list[[2]], ggplot2::aes(x = .data$part_name_from, y = .data$clust_name_from, xend = .data$part_name_to, yend = .data$clust_name_to, colour = .data$ecs, linewidth = .data$intersect_size), show.legend = TRUE) +
        ggplot2::scale_colour_gradientn(colours = edge_pal_colours, name = "ECS") +
        ggplot2::scale_linewidth_continuous(range = range_edge_width, name = "Intersection size") +
        ggnewscale::new_scale("colour") +
        ggplot2::geom_point(data = node_edge_list[[1]], aes(x = .data$part_name, y = .data$actual_y, colour = .data$ecc, size = .data$cluster_size), show.legend = TRUE) +
        ggplot2::scale_colour_gradientn(colours = node_pal_colours, name = "ECC") +
        # ggplot2::scale_colour_viridis_c(option = "F", direction = -1) +
        ggplot2::scale_size_continuous(range = range_point_size, name = "Cluster size") +
        ggplot2::scale_y_continuous(breaks = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank()
        )

    return(gplot_obj)
}

clust_hierarchical_plot <- function(clustering_assessment,
                                    clustering_method = NULL,
                                    k = NULL,
                                    edge_threshold = 0.3,
                                    range_point_size = c(1, 6),
                                    range_edge_width = c(0.01, 3),
                                    edge_palette_name = "RColorBrewer::Greys",
                                    edge_palette_inverse = FALSE,
                                    node_palette_name = "viridis::rocket",
                                    node_palette_inverse = TRUE) {
    clustering_assessment <- clustering_assessment$split_by_k
    cl_methods <- names(clustering_assessment)
    if (is.null(clustering_method) || !(clustering_method %in% cl_methods)) {
        clustering_method <- cl_methods[1]
    }
    clustering_assessment <- clustering_assessment[[clustering_method]]
    available_k <- names(clustering_assessment)
    if (is.null(k)) {
        k <- available_k
    } else {
        k <- intersect(k, available_k)
        if (length(k) == 0) {
            k <- available_k
        }
    }
    clustering_assessment <- clustering_assessment[k]
    ecc_list <- lapply(k, function(x) clustering_assessment[[x]]$ecc)
    names(ecc_list) <- k
    clustering_assessment <- data.frame(lapply(k, function(x) clustering_assessment[[x]]$partitions[[1]]$mb))
    colnames(clustering_assessment) <- k

    clust_hierarchical_plot_default(
        clustering_df = clustering_assessment,
        consistency_list = ecc_list,
        edge_threshold = edge_threshold,
        range_point_size = range_point_size,
        range_edge_width = range_edge_width,
        edge_palette_name = edge_palette_name,
        edge_palette_inverse = edge_palette_inverse,
        node_palette_name = node_palette_name,
        node_palette_inverse = node_palette_inverse
    )

}

clust_hierarchical_plot(clust_diff_obj, edge_palette_inverse =  TRUE)
