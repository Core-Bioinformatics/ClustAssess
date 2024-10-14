# expression_matrix <- NA
# metadata <- NA
# genes <- NA
# stab_obj <- NA
# feature_types <- NA

# pkg_env <- as.environment("namsespace:ClustAsssess")#new.env(parent = emptyenv())
# TODO check https://stackoverflow.com/questions/41954302/where-to-create-package-environment-variables   and   https://github.com/tidyverse/dplyr/blob/bbcfe99e29fe737d456b0d7adc33d3c445a32d9d/R/zzz.r   and  https://adv-r.hadley.nz/environments.html   http://adv-r.had.co.nz/Environments.html
pkg_env <- .GlobalEnv # new.env(parent = baseenv())
filetypes <- list(
    "PDF" = pdf,
    "PNG" = function(filename, width, height) {
        ragg::agg_png(filename, width, height, units = "in", res = 300)
    },
    "JPEG" = function(filename, width, height) {
        ragg::agg_jpeg(filename, width, height, units = "in", res = 300)
    },
    "SVG" = svg
)


add_env_variable <- function(variable_name, value) {
    assign(variable_name, value, envir = pkg_env)
}

get_env_var <- function(variable_name) {
    return(get(variable_name, envir = pkg_env))
}

mean_functions <- list(
    logNormalize = function(x, pseudocount_use = 1, base = 2) {
        return(log(x = Matrix::rowMeans(x = expm1(x = x)) + pseudocount_use, base = base))
    },
    default = function(x, pseudocount_use = 1, base = 2) {
        return(log(x = Matrix::rowMeans(x = x) + pseudocount_use, base = base))
    }
)

cList <- list(
    c(
        "grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
        "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"
    ),
    c(
        "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
        "#FEE090", "#FDAE61", "#F46D43", "#D73027"
    )[c(1, 1:9, 9)],
    c(
        "#FDE725", "#AADC32", "#5DC863", "#27AD81", "#21908C",
        "#2C728E", "#3B528B", "#472D7B", "#440154"
    )
)
names(cList) <- c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")

tags_style <- shiny::tags$head(
    shiny::tags$style(
        "#tabset_id {
        position: fixed;
        width: 100%;
        background-color: white;
        top: 0;
        z-index: 999;
        font-size: 25px;
        }",
        ".tab-content  {
        margin-top: 72px;
        }"
    )
)

empty_ggplot <- function() {
    ggplot2::ggplot() +
        ggplot2::theme_void()
}

grouped_boxplot_dataframe <- function(dataframe,
                                      y_column,
                                      x_column,
                                      plt_height,
                                      plt_width,
                                      boxplot_width = 0.5,
                                      xlabel = NULL,
                                      ylabel = NULL,
                                      plot_title = "",
                                      title_size = 1,
                                      axis_size = 1) {
    unique_groups <- unique(dataframe[, x_column])
    n_groups <- length(unique_groups)
    # groups_colours <- rhdf5::h5read("stability.h5", paste0("colors/", n_groups))
    groups_colours <- pkg_env$discrete_colors[[as.character(n_groups)]]

    if (is.null(xlabel)) {
        xlabel <- x_column
    }

    if (is.null(ylabel)) {
        ylabel <- y_column
    }

    if (axis_size > 1.5) {
        current_margins <- graphics::par("mai")
        current_margins[2] <- current_margins[2] + 0.5 * graphics::strheight(ylabel, units = "inches", cex = axis_size)
        graphics::par(mai = current_margins)
    }

    graphics::boxplot(
        formula = stats::as.formula(paste0(y_column, " ~ ", x_column)),
        data = dataframe,
        col = groups_colours,
        cex.main = title_size,
        cex.axis = axis_size,
        cex.lab = axis_size,
        boxwex = boxplot_width,
        main = plot_title,
        xlab = xlabel,
        ylab = ylabel
    )
}

grouped_boxplot_list <- function(groups_list,
                                 groups_values,
                                 x_values,
                                 plt_height,
                                 plt_width,
                                 xlab_text,
                                 boxplot_width = 0.5,
                                 space_inter = 1,
                                 space_intra = 1,
                                 display_legend = TRUE,
                                 text_size = 1) {
    unique_groups_values <- unique(groups_values)
    unique_x_values <- unique(x_values)

    # print(plt_width)
    # convert pixels to inches
    plt_height <- plt_height / ppi

    plt_width <- plt_width / ppi

    n_groups <- length(unique_groups_values)
    n_boxplots <- length(groups_values)

    # groups_colours <- rhdf5::h5read("stability.h5", paste0("colors/", n_groups))
    groups_values <- match(groups_values, unique_groups_values)
    groups_colours <- pkg_env$discrete_colors[[as.character(n_groups)]]
    at_values <- rep(0, length(groups_values))
    text_coords <- rep(0, length(unique_x_values))
    abline_coords <- rep(0, length(unique_x_values) - 1)

    count_diff <- -1
    prev_value <- -1
    updated_count <- FALSE
    start_index <- 1

    for (i in seq_along(groups_values)) {
        if (x_values[i] != prev_value) {
            if (prev_value != -1) {
                updated_count <- TRUE
            }
            prev_value <- x_values[i]
            count_diff <- count_diff + 1
        }

        at_values[i] <- count_diff * (n_groups * space_intra + space_inter) + groups_values[i] + (groups_values[i] - 1) * (space_intra - 1)

        if (updated_count) {
            # abline_coords[count_diff] <- (at_values[i] + at_values[i-1]) / 2
            abline_coords[count_diff] <- (count_diff - 1) * (n_groups * space_intra + space_inter) + n_groups + (n_groups - 1) * (space_intra - 1) + (space_inter + space_intra) / 2
            updated_count <- FALSE
            text_coords[count_diff] <- mean(at_values[start_index:(i - 1)])
            start_index <- i
        }
    }

    # print(at_values)
    # print(abline_coords)
    # print(groups_values)

    text_coords[count_diff + 1] <- mean(at_values[start_index:n_boxplots])

    if (display_legend) {
        predicted_width <- graphics::strwidth(c(" ", unique_groups_values), units = "inches", cex = text_size)
        space_width <- predicted_width[1]
        predicted_width <- predicted_width[2:length(predicted_width)]
        number_columns <- min(
            max(
                plt_width %/% (6 * space_width + max(predicted_width)),
                1
            ),
            length(unique_groups_values)
        )
        number_rows <- ceiling(length(unique_groups_values) / number_columns)

        text_height <- graphics::strheight(
            paste(
                rep("TEXT", number_rows + 1),
                collapse = "\n"
            ),
            units = "inches",
            cex = text_size
        ) + 0.1

        graphics::layout(
            matrix(c(1, 2), nrow = 2),
            heights = c(
                graphics::lcm((plt_height - text_height - 0.05) * 2.54),
                graphics::lcm(text_height * 2.54)
            )
        )
    }

    graphics::boxplot(
        groups_list,
        outline = FALSE,
        at = at_values,
        col = groups_colours[groups_values],
        boxwex = boxplot_width * (n_groups + (space_intra - 1) * (n_groups - 1)) / n_groups,
        xaxt = "n",
        xlab = NA,
        cex.axis = text_size,
        cex.lab = text_size
    )
    graphics::abline(v = abline_coords, lty = "dashed", col = "grey")
    graphics::title(xlab = xlab_text, ylab = "ecc", cex.lab = text_size)
    graphics::axis(side = 1, at = text_coords, labels = unique_x_values, las = 2, cex.axis = text_size)

    if (!display_legend) {
        return()
    }

    graphics::par(mar = c(0, 0, 0, 0))
    plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
    graphics::legend(
        "topleft",
        legend = unique_groups_values,
        col = groups_colours,
        pch = 15,
        cex = text_size,
        pt.cex = text_size * 2,
        bty = "n",
        # horizontal = TRUE,
        horiz = TRUE,
        # text.width = NA,
        xpd = TRUE
    )
}

expression_ggplot <- function(embedding, expression, threshold = 0) {
    df <- data.frame(embedding)
    colnames(df) <- paste("UMAP", 1:2, sep = "_")

    if (threshold > 0) {
        mask <- (expression > threshold)

        return(mask_ggplot(df, mask))
    }
    ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2, color = as.numeric(expression))
    ) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::scale_color_gradientn("", colors = cList[[1]]) +
        ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 15)) +
        ggplot2::theme(legend.position = "bottom")
}

color_ggplot <- function(embedding,
                         color_info,
                         cell_mask = NULL,
                         sort_cells = c("original", "highest", "lowest"),
                         labels = FALSE,
                         text_size = 20,
                         axis_text_size = 20,
                         legend_text_size = 20,
                         pt_size = 1) {
    df <- data.frame(embedding)
    colnames(df) <- paste("UMAP", 1:2, sep = "_")

    if (!is.null(cell_mask)) {
        color_info[!cell_mask] <- NA
    }
    is_continuous <- !(is.factor(color_info) || is.character(color_info))
    df$cell_colour <- color_info

    cell_ordering <- switch(sort_cells[1],
        "original" = seq_len(nrow(embedding)),
        "highest" = order(color_info, decreasing = FALSE, na.last = FALSE),
        "lowest" = order(color_info, decreasing = TRUE, na.last = FALSE)
    )

    df <- df[cell_ordering, ]

    ggplot_obj <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = .data$UMAP_1,
            y = .data$UMAP_2,
        )
    ) +
        ggplot2::geom_point(ggplot2::aes(color = .data$cell_colour), size = pt_size) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            legend.position = "bottom",
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = legend_text_size),
            axis.text = ggplot2::element_text(size = axis_text_size),
            axis.title = ggplot2::element_text(size = axis_text_size),
            plot.title = ggtext::element_textbox_simple(hjust = 0.5, size = axis_text_size * 1.5),
            aspect.ratio = 1
        )


    if (!labels || is_continuous) {
        return(ggplot_obj)
    }

    if (is.character(color_info)) {
        unique_values <- unique(color_info)
    } else {
        unique_values <- levels(droplevels(color_info))
    }

    medians_values <- t(sapply(
        seq_along(unique_values),
        function(i) {
            mask <- which(color_info == unique_values[i])
            emb <- embedding[mask, ]
            if (is.null(nrow(emb)) || nrow(emb) == 1) {
                return(emb)
            }

            gm_res <- Gmedian::Gmedian(emb)
            gm_res
        }
    ))
    medians_values <- data.frame(medians_values)
    colnames(medians_values) <- c("median_x", "median_y")
    medians_values$labels <- unique_values

    ggplot_obj + ggrepel::geom_text_repel(
        data = medians_values,
        ggplot2::aes(x = .data$median_x, y = .data$median_y, label = .data$labels),
        size = text_size,
        color = "white",
        bg.color = "black",
        bg.r = 0.15,
        labels.padding = 0.15
    )
}

metadata_plot <- function(embedding,
                          metadata_name,
                          plt_height,
                          plt_width,
                          groups_highlight = NULL,
                          pt_size = 1,
                          pch = ".",
                          text_size = 1,
                          axis_size = 1,
                          sort_cells = c("original", "highest", "lowest"),
                          display_legend = FALSE,
                          predicted_height = NULL,
                          labels = FALSE) {
    if (is.null(groups_highlight)) {
        metadata_mask <- rep(TRUE, nrow(embedding))
    } else {
        metadata_mask <- pkg_env$metadata[[metadata_name]] %in% groups_highlight
    }

    mtd_unique <- pkg_env$metadata_unique[[metadata_name]]

    color_plot2(
        embedding = embedding,
        color_info = pkg_env$metadata[[metadata_name]],
        # color_values = pkg_env$metadata_colors[[metadata_name]],
        color_values = pkg_env$discrete_colors[[as.character(length(mtd_unique))]],
        unique_values = mtd_unique,
        plt_height = plt_height,
        plt_width = plt_width,
        cell_mask = metadata_mask,
        pt_size = pt_size,
        pch = pch,
        text_size = text_size,
        axis_size = axis_size,
        sort_cells = sort_cells,
        display_legend = display_legend,
        predicted_height = predicted_height,
        labels = labels
    )
}

only_legend_plot <- function(unique_values,
                             color_values,
                             color_info,
                             plt_width,
                             text_size = 1) {
    is_continuous <- is.null(unique_values)
    plt_width <- plt_width / ppi

    if (is_continuous) {
        graphics::par(mai = c(0.1, 0, 0.1, 0))
        unique_values <- c(min(color_info), max(color_info))
    } else {
        graphics::par(mar = c(0, 0, 0, 0))
        # calculate space needed for the legend
        predicted_width <- graphics::strwidth(c(" ", unique_values), units = "inches", cex = text_size)
        space_width <- predicted_width[1]
        predicted_width <- predicted_width[2:length(predicted_width)]
        number_columns <- min(
            max(
                plt_width %/% (6 * space_width + max(predicted_width)),
                1
            ),
            length(unique_values)
        )
        number_rows <- ceiling(length(unique_values) / number_columns)

        # print(strheight(paste(
        #           rep("TEXT", number_rows + 1),
        #           collapse = "\n"
        #           ),
        #           units = "inches",
        #           cex = text_size) * ppi)
    }

    if (is.null(color_values)) {
        color_values <- paletteer::paletteer_c(palette = "viridis::viridis", n = 50)
    }

    if (is.function(color_values)) {
        ncolors <- ifelse(is_continuous, 50, length(unique_values))
        color_values <- color_values(ncolors)
    }

    if (is_continuous) {
        legend_image <- grDevices::as.raster(matrix(color_values, nrow = 1))
        plot(c(0, 1), c(-1, 1), type = "n", axes = FALSE, bty = "n", ylab = "", xlab = "")
        graphics::text(y = -0.5, x = seq(0, 1, l = 5), labels = round(seq(from = unique_values[1], to = unique_values[2], length.out = 5), digits = 3), cex = text_size)
        graphics::rasterImage(legend_image, 0, 0, 1, 1)
    } else {
        plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
        graphics::legend(
            "topleft",
            legend = unique_values,
            col = color_values,
            pch = 15,
            cex = text_size,
            pt.cex = text_size * 2,
            bty = "n",
            ncol = number_columns,
            # text.width = NA,
            xpd = TRUE
        )
    }
    # graphics::par(old_par)

    # grDevices::dev.off()
}


only_legend_metadata_plot <- function(metadata_name,
                                      groups = NULL,
                                      text_size = 1,
                                      plt_width) {
    unique_values <- pkg_env$metadata_unique[[metadata_name]]
    color_values <- pkg_env$discrete_colors[[metadata_name]]
    color_info <- pkg_env$metadata[, metadata_name]

    if (!is.null(unique_values) && !is.null(groups)) {
        color_values <- color_values[which(unique_values %in% groups)]
        unique_values <- groups
    }

    only_legend_plot(
        unique_values = unique_values,
        color_values = color_values,
        color_info = color_info,
        text_size = text_size,
        plt_width = plt_width
    )
}

color_plot2 <- function(embedding,
                        color_info,
                        color_values,
                        plt_height,
                        plt_width,
                        cell_mask = NULL,
                        groups_highlight = NULL,
                        unique_values = NULL,
                        pt_size = 1,
                        pch = ".",
                        text_size = 1,
                        legend_text_size = 1,
                        axis_size = 1,
                        sort_cells = c("original", "highest", "lowest"),
                        display_legend = FALSE,
                        predicted_height = NULL,
                        labels = FALSE) {
    # xlim <- c(min(embedding[, 1]), max(embedding[, 1]))
    # ylim <- c(min(embedding[, 2]), max(embedding[, 2]))

    # convert pixels to inches
    plt_height <- plt_height / ppi
    plt_width <- plt_width / ppi

    is_continuous <- is.null(unique_values)
    is_logical <- is.logical(color_info)

    if (!is.null(cell_mask)) {
        color_info[!cell_mask] <- NA

        if (!is_continuous && !is_logical) {
            unique_values <- unique(color_info[cell_mask])
        }
    }

    if (!is_continuous && !is_logical) {
        sort_cells <- ifelse(is.null(cell_mask), "original", "mask")
    }

    cell_ordering <- switch(sort_cells[1],
        "original" = seq_len(nrow(embedding)),
        "highest" = order(color_info, decreasing = FALSE, na.last = FALSE),
        "lowest" = order(color_info, decreasing = TRUE, na.last = FALSE),
        "mask" = c(which(!cell_mask), which(cell_mask))
    )

    if (is_logical) {
        color_info <- as.character(color_info)
    }

    if (is.null(color_values)) {
        color_values <- paletteer::paletteer_c(palette = "viridis::viridis", n = 50)
    }

    if (is.function(color_values)) {
        ncolors <- ifelse(is_continuous, 50, length(unique_values))
        color_values <- color_values(ncolors)
    }

    if (display_legend) {
        if (is_continuous) {
            unique_values <- c(min(color_info), max(color_info))
            number_rows <- 2

            if (is.null(predicted_height)) {
                predicted_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = legend_text_size)
            }
        } else {
            # calculate space needed for the legend
            predicted_width <- graphics::strwidth(unique_values, units = "inches", cex = legend_text_size)
            space_width <- graphics::strwidth(" ", units = "inches", cex = legend_text_size)
            number_columns <- min(
                max(
                    plt_width %/% (6 * space_width + max(predicted_width)),
                    1
                ),
                length(unique_values)
            )


            if (is.null(predicted_height)) {
                number_rows <- ceiling(length(unique_values) / number_columns)
                predicted_height <- graphics::strheight(
                    paste(
                        rep("TEXT", number_rows + 1),
                        collapse = "\n"
                    ),
                    units = "inches",
                    cex = legend_text_size
                )
            }
        }
    }

    if (is_continuous) {
        color_info <- cut(color_info, breaks = 50)
    }

    colrs <- color_values[color_info]
    if (!is.null(cell_mask)) {
        colrs[!cell_mask] <- "ivory4" 
    }

    if (display_legend) {
        graphics::layout(
            matrix(c(1, 2), nrow = 2),
            heights = c(
                graphics::lcm((plt_height - predicted_height) * 2.5),
                graphics::lcm(predicted_height * 2.5)
            )
        )
    }

    plot(
        embedding[cell_ordering, ],
        pch = pch,
        col = colrs[cell_ordering],
        cex = pt_size,
        xlab = "UMAP_1",
        ylab = "UMAP_2",
        cex.axis = axis_size,
        cex.lab = axis_size
    )

    if (!is_continuous && labels) {
        for (unique_val in unique_values) {
            mask <- which(color_info == unique_val)
            emb <- embedding[mask, ]

            if (is.null(nrow(emb)) || nrow(emb) == 1) {
                geometric_median <- emb
            } else {
                geometric_median <- Gmedian::Gmedian(embedding[mask, ])
            }

            shadowtext(
                geometric_median[1], #+ text_width * 0.6,
                geometric_median[2], #+ text_size * 0.75,
                unique_val,
                col = "black",
                bg = "white",
                cex = text_size,
                xpd = TRUE
            )
        }
    }


    if (!display_legend) {
        return()
    }

    if (is_continuous) {
        graphics::par(mai = c(0.1, 0, 0.1, 0))
        legend_image <- grDevices::as.raster(matrix(color_values, nrow = 1))
        plot(c(0, 1), c(-1, 1), type = "n", axes = F, bty = "n", ylab = "", xlab = "")
        graphics::text(y = -0.5, x = seq(0, 1, l = 5), labels = round(seq(from = unique_values[1], to = unique_values[2], length.out = 5), digits = 3), cex = legend_text_size)
        graphics::rasterImage(legend_image, 0, 0, 1, 1)
    } else {
        graphics::par(mar = c(0, 0, 0, 0))
        plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
        graphics::legend(
            "topleft",
            legend = unique_values,
            col = color_values,
            pch = 15,
            cex = legend_text_size,
            pt.cex = legend_text_size * 2,
            bty = "n",
            ncol = number_columns,
            # text.width = NA,
            xpd = TRUE
        )
    }

    # grDevices::dev.off()
}

color_plot <- function(embedding,
                       color_info,
                       pt_size = 1,
                       pch = ".",
                       color_palette = NULL,
                       labels = FALSE,
                       ncolors = NULL,
                       text_size = 1) {
    if (is.null(color_palette)) {
        color_palette <- ifelse(is.null(ncolors), "palettesForR::Named", "viridis::viridis")
    }
    paletteer_function <- paletteer::paletteer_c
    if (is.null(ncolors)) {
        if (is.factor(color_info)) {
            unique_values <- levels(color_info)
        } else {
            unique_values <- unique(color_info)
        }
        ncolors <- length(unique_values)
        paletteer_function <- paletteer::paletteer_d
    } else {
        color_info <- cut(color_info, breaks = ncolors)
    }
    colors <- paletteer_function(
        palette = color_palette,
        n = ncolors
    )

    return({
        plot(
            embedding[, 1],
            embedding[, 2],
            pch = pch,
            col = colors[color_info],
            cex = pt_size,
            panel.first = {
                graphics::axis(1, tck = 1, lty = 2, col = "gray")
                graphics::axis(2, tck = 1, lty = 2, col = "gray")
            },
            xlab = "UMAP_1",
            ylab = "UMAP_2"
        )

        if (labels) {
            for (unique_val in unique_values) {
                geometric_median <- Gmedian::Gmedian(embedding[color_info == unique_val, ])

                text_width <- graphics::strwidth(unique_val, cex = text_size)
                graphics::rect(
                    geometric_median[1],
                    geometric_median[2],
                    geometric_median[1] + text_width * 1.2,
                    geometric_median[2] + text_size * 1.2,
                    col = "white",
                    border = "black"
                )

                graphics::text(
                    geometric_median[1] + text_width * 0.6,
                    geometric_median[2] + text_size * 0.6,
                    unique_val,
                    cex = text_size,
                    col = "black"
                )
            }
        }
    })
}

color_c_plot <- function(embedding,
                         color_info,
                         pt_size = 1,
                         pch = ".",
                         color_palette = "viridis::viridis",
                         text_size = 1) {
    colors <- paletteer::paletteer_c(
        palette = color_palette,
        n = 50
    )

    return({
        plot(
            embedding[, 1],
            embedding[, 2],
            pch = pch,
            col = colors[color_info],
            cex = pt_size,
            panel.first = {
                graphics::axis(1, tck = 1, lty = 2, col = "gray")
                graphics::axis(2, tck = 1, lty = 2, col = "gray")
            },
            xlab = "UMAP_1",
            ylab = "UMAP_2"
        )
    })
}



fill_ggplot <- function(embedding, fill_info) {
    ggplot2::ggplot(
        data.frame(embedding),
        ggplot2::aes(
            x = .data$UMAP_1,
            y = .data$UMAP_2,
            fill = fill_info
        )
    ) +
        ggplot2::geom_point() +
        ggplot2::theme_bw()
}

mask_ggplot <- function(embedding_df, mask) {
    ggplot2::ggplot() +
        ggplot2::geom_point(
            data = embedding_df[!mask, ],
            mapping = ggplot2::aes(x = .data$UMAP_1, .data$UMAP_2),
            color = "#C3C3C3"
        ) +
        ggplot2::geom_point(
            data = embedding_df[mask, ],
            mapping = ggplot2::aes(x = .data$UMAP_1, .data$UMAP_2),
            color = "red"
        ) +
        ggplot2::theme_bw()
}

voting_scheme_ggplot <- function(embedding,
                                 n_expressed_genes,
                                 n_genes_above_threshold) {
    umap_df <- data.frame(embedding)
    colnames(umap_df) <- paste("UMAP", 1:2, sep = "_")
    # n_expressed_genes <- rep(0, ncol(pkg_env$expression_matrix))
    # expression_threshold <- as.numeric(expression_threshold)
    # for (gene in selected_genes) {
    #   n_expressed_genes <- n_expressed_genes + (pkg_env$expression_matrix[gene, ] > expression_threshold)
    # }

    mask <- (n_expressed_genes >= n_genes_above_threshold)

    mask_ggplot(umap_df, mask)
}

metadata_ggplot <- function(embedding,
                            metadata_name) {
    metadata_value <- pkg_env$metadata[, metadata_name]
    if (is.factor(metadata_value)) {
        return(color_ggplot(embedding, metadata_value) +
            ggplot2::guides(color = ggplot2::guide_legend(title = metadata_name)))
    }

    if (!is.numeric(metadata_value)) {
        metadata_value <- as.numeric(factor(metadata_value))
    }

    color_ggplot(embedding, metadata_value) +
        ggplot2::scale_color_continuous(metadata_name)
}

download_plot_modal <- function(session, output, placeholder, download_id) {
    ns <- session$ns
    shiny::showModal(
        shiny::modalDialog(
            shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
            shiny::textInput(ns("filename"), "File name:", placeholder, width = "80%"),
            shiny::numericInput(ns("width"), "Width (in):", 7, 3, 100, 0.1),
            shiny::numericInput(ns("height"), "Height (in):", 7, 3, 100, 0.1),
            footer = shiny::tagList(
                shiny::modalButton("Cancel"),
                shiny::downloadButton(outputId = ns(download_id))
            ),
            easyClose = TRUE
        ),
        session = session
    )
}

download_plot_handler <- function(session, filename, to_save_plot, width, height) {
    extension <- stringr::str_split(filename, "\\.")[[1]]
    if (length(extension) == 1) {
        filename <- glue::glue("{filename}.pdf")
    }

    shiny::downloadHandler(
        filename = function() {
            filename
        },
        content = function(file) {
            on.exit(shiny::removeModal(session))
            ggplot2::ggsave(file, to_save_plot, width = width, height = height)
        }
    )
}

identify_cluster_markers <- function(
    cells1_name,
    cells2_name,
    expression_matrix,
    mean_function_name = "logNormalize",
    features = NULL) {
    foldchange_result <- Seurat::FoldChange(
        object = expression_matrix,
        cells.1 = cells1_name,
        cells.2 = cells2_name,
        features = features,
        mean.fxn = mean_functions[[mean_function_name]],
        fc.name = ifelse(mean_function_name == "rowMeans", "avg_diff", "avg_log2FC")
    )

    markers_result <- Seurat::FindMarkers(
        object = pkg_env$expression_matrix,
        slot = "data",
        cells.1 = cells1_name,
        cells.2 = cells2_name,
        fc.results = foldchange_result
    )

    return(markers_result)
}


create_seurat_object <- function(
    features,
    feature_type,
    clustering_algorithm,
    nclusters,
    ndims,
    calculate_dim_reduction = TRUE) {
    nfeatures <- as.character(length(features))

    so <- SeuratObject::CreateSeuratObject(
        counts = pkg_env$expression_matrix,
        meta.data = pkg_env$metadata
    )

    so <- Seurat::NormalizeData(so, verbose = FALSE)
    so <- Seurat::FindVariableFeatures(so, verbose = FALSE)
    so <- Seurat::ScaleData(so, verbose = FALSE)
    if (calculate_dim_reduction) {
        so <- Seurat::RunPCA(so, npcs = ndims, verbose = FALSE)
        so <- Seurat::RunUMAP(so, reduction = "pca", dims = seq_len(ndims), verbose = FALSE)
    }

    for (nclust in nclusters) {
        so[[glue::glue("stable_{nclust}_clusters")]] <- pkg_env$stab_obj[[feature_type]][[nfeatures]]$clustering_importance$split_by_k[[clustering_algorithm]][[nclust]]$partitions[[1]]$mb
        so[[glue::glue("stable_{nclust}_ecc")]] <- pkg_env$stab_obj[[feature_type]][[nfeatures]]$clustering_importance$split_by_k[[clustering_algorithm]][[nclust]]$ecc
    }

    Seurat::Idents(so) <- glue::glue("stable_{nclusters[1]}_clusters")

    return(so)
}

render_plot_by_height <- function(id, session) {
    session$output[[paste0(id, "_generator")]] <- shiny::renderUI({
        used_height <- floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2))
        shiny::plotOutput(
            outputId = session$ns(id),
            height = paste0(used_height, "px")
        )
    })
}

#' Calculate markers - Shiny
#'
#' @description Performs the Wilcoxon rank sum test to identify differentially
#' expressed genes between two groups of cells in the shiny context. The
#' method can be also used outside the shiny context, as long as the expression
#' matrix is stored in a h5 file.
#'
#' @param cells1 A vector of cell indices for the first group of cells.
#' @param cells2 A vector of cell indices for the second group of cells.
#' @param logfc_threshold The minimum absolute log fold change to consider a
#' gene as differentially expressed. Defaults to `0`, meaning all genes are
#' taken into considereation.
#' @param min_pct_threshold The minimum fraction of cells expressing a gene
#' form each cell population to consider the gene as differentially expressed.
#' Increasing the value will speed up the function. Defaults to `0.1`.
#' @param average_expression_threshold The minimum average expression that a
#' gene should have in order to be considered as differentially expressed.
#' @param average_expression_group1_threshold The minimum average expression
#' that a gene should have in the first group of cells to be considered as
#' differentially expressed. Defaults to `0`.
#' @param min_diff_pct_threshold The minimum difference in the fraction of cells
#' expressing a gene between the two cell populations to consider the gene as
#' differentially expressed. Defaults to `-Inf`.
#' @param used_slot Parameter that provides additional information about the
#' expression matrix, whether it was scaled or not. The value of this parameter
#' impacts the calculation of the fold change. If `data`, the function will
#' calculates the fold change as the fraction between the log value of the
#' average of the expression raised to exponential for the two cell groups. If
#' `scale.data`, the function will calculate the fold change as the fraction
#' between the average of the expression values for the two cell groups.
#' Other options will default to calculating the fold change as the fraction
#' between the log value of the average of the expression values for the two
#' cell groups. Defaults to `data`.
#' @param norm_method The normalization method used to normalize the expression
#' matrix. The value of this parameter impacts the calculation of the average
#' expression of the genes when `used_slot = "data"`. If `LogNormalize`, the
#' log fold change will be calculated as described for the `used_slot`
#' parameter. Otherwise, the log fold change will be calculated as the fraction
#' between the log value of the average of the expression values for the two
#' cell groups. Defaults to `SCT`.
#' @param expression_h5_path The path to the h5 file containing the expression
#' matrix. The h5 file should contain the following fields: `expression_matrix`,
#' `rank_matrix`, `average_expression`, `genes`. The file path
#' defaults to `expression.h5`.
#' @param pseudocount_use The pseudocount to add to the expression values when
#' calculating the average expression of the genes, to avoid the 0 value for
#' the denominator. Defaults to `1`.
#' @param base The base of the logharithm. Defaults to `2`.
#' @param verbose Whether to print messages about the progress of the function.
#' Defaults to TRUE.
#' @param check_difference Whether to perform set difference between the two
#' cells. Defaults to TRUE.

#' @return A data frame containing the following columns:
#' - `gene`: The gene name.
#' - `avg_log2FC`: The average log fold change between the two cell groups.
#' - `p_val`: The p-value of the Wilcoxon rank sum test.
#' - `p_val_adj`: The adjusted p-value of the Wilcoxon rank sum test.
#' - `pct.1`: The fraction of cells expressing the gene in the first cell group.
#' - `pct.2`: The fraction of cells expressing the gene in the second cell group.
#' - `avg_expr_group1`: The average expression of the gene in the first cell group.
#' - `avg_expr`: The average expression of the gene.
#'
#' @export
calculate_markers_shiny <- function(cells1,
                                    cells2,
                                    logfc_threshold = 0,
                                    min_pct_threshold = 0.1,
                                    average_expression_threshold = 0,
                                    average_expression_group1_threshold = 0,
                                    min_diff_pct_threshold = -Inf,
                                    used_slot = "data",
                                    norm_method = "SCT",
                                    expression_h5_path = "expression.h5",
                                    pseudocount_use = 1,
                                    base = 2,
                                    verbose = TRUE,
                                    check_difference = TRUE) {
    if (check_difference) {
        cells2 <- setdiff(cells2, cells1)
    }

    if (length(cells2) == 0 || length(cells1) == 0) {
        warning("One of the cell groups is empty!")
        return()
    }

    average_expressions <- rhdf5::h5read(expression_h5_path, "average_expression")
    nfiltered_genes <- sum(average_expressions >= average_expression_threshold)
    chunk_size <- as.integer(rhdf5::h5read(expression_h5_path, "chunk_size"))
    genes <- rhdf5::h5read(expression_h5_path, "genes")[seq_len(nfiltered_genes)]

    average_expressions <- average_expressions[seq_len(nfiltered_genes)]
    names(average_expressions) <- genes

    nchunks <- ceiling(nfiltered_genes / chunk_size)
    if (verbose) {
        print(glue::glue("[{Sys.time()}] Start DEG analysis"))
    }
    df_list <- lapply(
        seq_len(nchunks),
        function(i) {
            # print(i)
            if (i == nchunks) {
                index_list <- seq(from = chunk_size * (i - 1) + 1, by = 1, to = nfiltered_genes)
            } else {
                index_list <- seq(from = chunk_size * (i - 1) + 1, by = 1, length.out = chunk_size)
            }
            calculate_markers(
                expression_matrix = rhdf5::h5read(
                    expression_h5_path,
                    "expression_matrix",
                    index = list(index_list, NULL)
                ),
                rank_matrix = rhdf5::h5read(
                    expression_h5_path,
                    "rank_matrix",
                    index = list(index_list, NULL)
                ),
                cells1 = cells1,
                cells2 = cells2,
                logfc_threshold = logfc_threshold,
                min_pct_threshold = min_pct_threshold,
                min_diff_pct_threshold = min_diff_pct_threshold,
                avg_expr_threshold_group1 = average_expression_group1_threshold,
                feature_names = genes[index_list],
                used_slot = used_slot,
                norm_method = norm_method,
                pseudocount_use = pseudocount_use,
                base = base,
                check_cells_set_diff = FALSE,
                adjust_pvals = FALSE
            )
        }
    )
    df_list <- dplyr::bind_rows(df_list)
    df_list$p_val_adj <- stats::p.adjust(
        p = df_list$p_val,
        method = "bonferroni",
        n = nfiltered_genes
    )
    df_list$avg_expr <- average_expressions[df_list$gene]

    if (verbose) {
        print(glue::glue("[{Sys.time()}] DEG analysis - Done"))
    }

    df_list
}

#' Calculate markers
#'
#' @description Performs the Wilcoxon rank sum test to identify differentially
#' expressed genes between two groups of cells.
#'
#' @param expression_matrix A matrix of gene expression values having genes
#' in rows and cells in columns.
#' @param cells1 A vector of cell indices for the first group of cells.
#' @param cells2 A vector of cell indices for the second group of cells.
#' @param logfc_threshold The minimum absolute log fold change to consider a
#' gene as differentially expressed. Defaults to `0`, meaning all genes are
#' taken into considereation.
#' @param min_pct_threshold The minimum fraction of cells expressing a gene
#' form each cell population to consider the gene as differentially expressed.
#' Increasing the value will speed up the function. Defaults to `0.1`.
#' @param avg_expr_threshold_group1 The minimum average expression that a gene
#' should have in the first group of cells to be considered as differentially
#' expressed. Defaults to `0`.
#' @param min_diff_pct_threshold The minimum difference in the fraction of cells
#' expressing a gene between the two cell populations to consider the gene as
#' differentially expressed. Defaults to `-Inf`.
#' @param rank_matrix A matrix where the cells are ranked based on their
#' expression levels with respect to each gene. Defaults to `NULL`, in which
#' case the function will calculate the rank matrix. We recommend calculating
#' the rank matrix beforehand and passing it to the function to speed up the
#' computation.
#' @param feature_names A vector of gene names. Defaults to `NULL`, in which
#' case the function will use the row names of the expression matrix as gene
#' names.
#' @param used_slot Parameter that provides additional information about the
#' expression matrix, whether it was scaled or not. The value of this parameter
#' impacts the calculation of the fold change. If `data`, the function will
#' calculates the fold change as the fraction between the log value of the
#' average of the expression raised to exponential for the two cell groups. If
#' `scale.data`, the function will calculate the fold change as the fraction
#' between the average of the expression values for the two cell groups.
#' Other options will default to calculating the fold change as the fraction
#' between the log value of the average of the expression values for the two
#' cell groups. Defaults to `data`.
#' @param norm_method The normalization method used to normalize the expression
#' matrix. The value of this parameter impacts the calculation of the average
#' expression of the genes when `used_slot = "data"`. If `LogNormalize`, the
#' log fold change will be calculated as described for the `used_slot`
#' parameter. Otherwise, the log fold change will be calculated as the fraction
#' between the log value of the average of the expression values for the two
#' cell groups. Defaults to `SCT`.
#' @param pseudocount_use The pseudocount to add to the expression values when
#' calculating the average expression of the genes, to avoid the 0 value for
#' the denominator. Defaults to `1`.
#' @param base The base of the logharithm. Defaults to `2`.
#' @param adjust_pvals A logical value indicating whether to adjust the p-values
#' for multiple testing using the Bonferonni method. Defaults to `TRUE`.
#' @param check_cells_set_diff A logical value indicating whether to check if
#' thw two cell groups are disjoint or not. Defaults to `TRUE`.
#'
#' @return A data frame containing the following columns:
#' - `gene`: The gene name.
#' - `avg_log2FC`: The average log fold change between the two cell groups.
#' - `p_val`: The p-value of the Wilcoxon rank sum test.
#' - `p_val_adj`: The adjusted p-value of the Wilcoxon rank sum test.
#' - `pct.1`: The fraction of cells expressing the gene in the first cell group.
#' - `pct.2`: The fraction of cells expressing the gene in the second cell group.
#' - `avg_expr_group1`: The average expression of the gene in the first cell group.
#'
#' @export
#' @examples
#' set.seed(2024)
#' # create an artificial expression matrix
#' expr_matrix <- matrix(
#'     c(runif(100 * 50), runif(100 * 50, min = 3, max = 4)),
#'     ncol = 200, byrow = FALSE
#' )
#' colnames(expr_matrix) <- as.character(1:200)
#' rownames(expr_matrix) <- paste("feature", 1:50)
#'
#' calculate_markers(
#'     expression_matrix = expr_matrix,
#'     cells1 = 101:200,
#'     cells2 = 1:100
#' )
#' # TODO should be rewritten such that you don't create new matrix objects inside
#' # just
calculate_markers <- function(expression_matrix,
                              cells1,
                              cells2,
                              logfc_threshold = 0,
                              min_pct_threshold = 0.1,
                              avg_expr_threshold_group1 = 0,
                              min_diff_pct_threshold = -Inf,
                              rank_matrix = NULL,
                              feature_names = NULL,
                              used_slot = "data",
                              norm_method = "SCT",
                              pseudocount_use = 1,
                              base = 2,
                              adjust_pvals = TRUE,
                              check_cells_set_diff = TRUE) {
    rowMeans_function <- base::rowMeans
    if (inherits(expression_matrix, "dgCMatrix")) {
        rowMeans_function <- Matrix::rowMeans
    }

    # as of Seurat 5
    deprecated_default_mean_fxn <- function(x) {
        return(log(x = rowMeans_function(x = x) + pseudocount_use, base = base))
    }

    sct_data_mean_function <- function(x) {
       return(log(x = (rowSums(x = expm1(x = x)) + pseudocount_use)/NCOL(x), base = base))
    }

    counts_mean_function <- function(x) {
        return(log(x = (rowSums(x = x) + pseudocount_use)/NCOL(x), base = base))
    }

    log1pdata_mean_function <- function(x) {
        return(log(x = (rowSums(x = expm1(x = x)) + pseudocount_use)/NCOL(x), base = base))
    }

    if (is.null(feature_names)) {
        if (is.null(rownames(expression_matrix))) {
            rownames(expression_matrix) <- paste0("gene_", seq_len(nrow(expression_matrix)))
        }

        feature_names <- rownames(expression_matrix)
    } else {
        rownames(expression_matrix) <- feature_names
    }

    if (check_cells_set_diff) {
        cells2 <- setdiff(cells2, cells1)
    }

    if (length(cells2) == 0 || length(cells1) == 0) {
        return()
    }

    indices <- c(cells1, cells2)
    n <- length(feature_names)

    expression_matrix <- expression_matrix[, indices, drop = FALSE]

    base_text <- ifelse(test = base == exp(1), yes = "", no = base)
    fc_name <- ifelse(
        test = used_slot == "scale.data",
        yes = "avg_diff",
        no = paste0("avg_log", base_text, "FC")
    )
    mean_fxn <- switch(EXPR = used_slot,
        counts = counts_mean_function,
        data = switch(EXPR = norm_method,
            SCT = sct_data_mean_function,
            LogNormalize = log1pdata_mean_function,
            NULL = log1pdata_mean_function,
            counts_mean_function
        ),
        scale.data = rowMeans_function,
        counts_mean_function
    )

    # TODO a bit redundant; this is also done in FoldChange; check if you can combine the two steps
    avg_expr_group1 <- rowMeans_function(expression_matrix[, seq_along(cells1), drop = FALSE])
    # TODO calculate the average of the second group as well; should we perform the filter here as well?
    if (avg_expr_threshold_group1 > 0) {
        mask <- which(avg_expr_group1 >= avg_expr_threshold_group1)

        if (length(mask) == 0) {
            return(data.frame())
        }

        expression_matrix <- expression_matrix[mask, , drop = FALSE]
        feature_names <- feature_names[mask]
        avg_expr_group1 <- avg_expr_group1[mask]
    }

    # TODO create own calculation of FoldChange
    fc_results <- Seurat::FoldChange(
        object = expression_matrix,
        cells.1 = seq_along(cells1),
        cells.2 = seq(from = length(cells1) + 1, to = length(indices), by = 1),
        mean.fxn = mean_fxn,
        fc.name = fc_name
    )

    max_pct <- pmax(fc_results$pct.1, fc_results$pct.2)
    min_pct <- pmin(fc_results$pct.1, fc_results$pct.2)
    pct_diff <- max_pct - min_pct

    mask <- which(max_pct >= min_pct_threshold &
        pct_diff >= min_diff_pct_threshold &
        abs(fc_results[, 1]) >= logfc_threshold)

    if (length(mask) == 0) {
        return(data.frame())
    }

    expression_matrix <- expression_matrix[mask, , drop = FALSE]
    fc_results <- fc_results[mask, ]
    feature_names <- feature_names[mask]
    avg_expr_group1 <- avg_expr_group1[mask]

    if (is.null(rank_matrix)) {
        rank_matrix <- matrix(nrow = nrow(expression_matrix), ncol = ncol(expression_matrix))

        for (i in seq_len(nrow(expression_matrix))) {
            rank_matrix[i, ] <- rank(expression_matrix[i, ], ties.method = "min")
        }
    } else {
        rank_matrix <- rank_matrix[mask, indices, drop = FALSE]
    }

    if (length(mask) == 1) {
        # expression_matrix <- matrix(expression_matrix, nrow = 1)
        rank_matrix <- matrix(rank_matrix, nrow = 1)
    }

    fc_results$gene <- feature_names
    fc_results$p_val <- wilcox_test(rank_matrix, length(cells1), max(rank_matrix))

    if (adjust_pvals) {
        fc_results$p_val_adj <- stats::p.adjust(
            p = fc_results$p_val,
            method = "bonferroni",
            n = n
        )

        fc_results <- fc_results[, c(4, 1, 2, 3, 5, 6)]
        fc_results$avg_expr_group1 <- avg_expr_group1

        return(fc_results)
    }

    fc_results <- fc_results[, c(4, 1, 2, 3, 5)]
    fc_results$avg_expr_group1 <- avg_expr_group1

    return(fc_results)
}

# function copied from https://cran.r-project.org/web/packages/TeachingDemos/index.html
shadowtext <- function(x, y = NULL, labels, col = "white", bg = "black",
                       theta = seq(pi / 32, 2 * pi, length.out = 64), r = 0.1, cex = 1, ...) {
    xy <- grDevices::xy.coords(x, y)
    fx <- graphics::grconvertX(xy$x, to = "nfc")
    fy <- graphics::grconvertY(xy$y, to = "nfc")
    fxo <- r * graphics::strwidth("A", units = "figure", cex = cex)
    fyo <- r * graphics::strheight("A", units = "figure", cex = cex)

    for (i in theta) {
        graphics::text(graphics::grconvertX(fx + cos(i) * fxo, from = "nfc"),
            graphics::grconvertY(fy + sin(i) * fyo, from = "nfc"),
            labels,
            cex = cex, col = bg, ...
        )
    }
    graphics::text(xy$x, xy$y, labels, cex = cex, col = col, ...)
}

gear_overall <- function(ns, id) {
    shiny::tagList(
        shiny::sliderInput(
            inputId = ns(paste0(id, "_boxplot_width")),
            label = "Boxplot width",
            min = 0.00, max = 1.00, value = 0.50
        ),
        shiny::sliderInput(
            inputId = ns(paste0(id, "_inter_distance")),
            label = "Space between groups",
            min = 1, max = 15, value = 1, step = 1
        ),
        shiny::sliderInput(
            inputId = ns(paste0(id, "_intra_distance")),
            label = "Space inside groups",
            min = 1, max = 15, value = 1, step = 1
        ),
        shiny::sliderInput(
            inputId = ns(paste0(id, "_text_size")),
            label = "Text size",
            min = 0.10, max = 10.00, value = 1.00, step = 0.1
        )
    )
}

gear_umaps <- function(ns, id, discrete = TRUE, default_order = "original") {
    shinyWidgets::dropdownButton(
        shiny::tagList(
            if (discrete) {
                shiny::sliderInput(
                    inputId = ns(paste0(id, "_text_size")),
                    label = "Text size",
                    min = 0.10, max = 10.00, value = 1.00, step = 0.1
                )
            },
            shiny::sliderInput(
                inputId = ns(paste0(id, "_axis_size")),
                label = "Axis labels size",
                min = 0.10, max = 10.00, value = 1.00, step = 0.1
            ),
            shiny::sliderInput(
                inputId = ns(paste0(id, "_legend_size")),
                label = "Legend size",
                min = 0.10, max = 10.00, value = 1.00, step = 0.1
            ),
            shiny::sliderInput(
                inputId = ns(paste0(id, "_pt_size")),
                label = "Point size",
                min = 0.05, max = 5.00, value = 0.30, step = 0.05
            ),
            shinyWidgets::radioGroupButtons(
                inputId = ns(paste0(id, "_pt_type")),
                label = "Point type",
                choices = c("Circle", "Pixel")
            ),
            shinyWidgets::radioGroupButtons(
                inputId = ns(paste0(id, "_pt_order")),
                label = "Point ordering",
                choices = c("original", "highest", "lowest"),
                selected = default_order
            ),
            if (discrete) {
                shinyWidgets::prettySwitch(
                    inputId = ns(paste0(id, "_labels")),
                    label = "Show labels",
                    status = "success",
                    fill = TRUE
                )
            }
        ),
        circle = TRUE,
        status = "success",
        size = "sm",
        icon = shiny::icon("cog")
    )
}

gear_download <- function(ns, id, label = "") {
    shinyWidgets::dropdownButton(
        label = "",
        icon = shiny::icon("download"),
        status = "success",
        size = "sm",
        shiny::em(paste0("Note: Use one of the following extensions: ", paste(names(filetypes), collapse = ", "))),
        shiny::textInput(ns(paste0("filename_", id)), "File name:", width = "80%", value = label),
        shiny::numericInput(ns(paste0("width_", id)), "Width (in):", 7, 3, 100, 0.1),
        shiny::numericInput(ns(paste0("height_", id)), "Height (in):", 7, 3, 100, 0.1),
        shiny::selectInput(ns(paste0("filetype_", id)), "Filetype", choices = names(filetypes), selected = names(filetypes)[1], width = "80%"),
        shinyWidgets::prettyRadioButtons(ns(paste0("raster_", id)), "Rasterize", choices = c("Yes", "No"), selected = "Yes"),
        shiny::downloadButton(ns(paste0("download_", id)), label = "Download Plot")
    )
}

jaccard_index <- function(a, b) {
    if (length(a) == 0 && length(b) == 0) {
        return(1)
    } else {
        return(length(intersect(a, b)) / length(union(a, b)))
    }
}

update_gears_width <- function() {
    shiny::observe({
        win_dims <- pkg_env$dimension()
        shiny::req(win_dims)

        shinyjs::runjs(paste0(
            "$('.dropdown-menu .shiny-split-layout').css('width', '", win_dims[1] / 2.2, "px');"
        ))
    })
}

gene_name_transformation <- function(gene_name) {
    gene_name <- toupper(gene_name)
    gene_name <- gsub(" ", "", gene_name)
    gene_name <- gsub("-", "", gene_name)
    gene_name <- gsub("_", "", gene_name)
    gene_name <- gsub("\\.", "", gene_name)

    return(gene_name)
}

# split_vector_by_metadata <- function(vec, metadata) {
#     if (is.factor(metadata)) {
#         unique_values <- levels(metadata)
#     } else {
#         unique_values <- unique(metadata)
#     }

#     split_list <- lapply(unique_values, function(i) { rep(0, length(vec))})
#     iter_list <- rep(1, length(unique_values))
#     names(split_list) <- unique_values
#     names(iter_list) <- unique_values

#     for (i in)





# }
