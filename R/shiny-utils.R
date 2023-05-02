
# expression_matrix <- NA
# metadata <- NA
# genes <- NA
# stab_obj <- NA
# feature_types <- NA

# pkg_env <- as.environment("namsespace:ClustAsssess")#new.env(parent = emptyenv())
pkg_env <- .GlobalEnv#new.env(parent = baseenv())


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
                                      xlabel = NULL,
                                      ylabel = NULL,
                                      plot_title = "",
                                      text_size = 1) {

  
  unique_groups <- unique(dataframe[ , x_column])
  n_groups <- length(unique_groups)
  groups_colours <- rhdf5::h5read("stability.h5", paste0("colors/", n_groups))

  if (is.null(xlabel)) {
    xlabel <- x_column
  }

  if (is.null(ylabel)) {
    ylabel <- y_column
  }
  
  boxplot(
    formula = as.formula(paste0(y_column, " ~ ", x_column)),
    data = dataframe,
    col = groups_colours,
    cex.axis = text_size,
    cex.label = text_size,
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

  groups_colours <- rhdf5::h5read("stability.h5", paste0("colors/", n_groups))
  groups_values <- match(groups_values, unique_groups_values)
  at_values <- rep(0, length(groups_values))
  text_coords <- rep(0, length(unique_x_values))
  abline_coords <- rep(0, length(unique_x_values)-1)
  
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
      text_coords[count_diff] <- mean(at_values[start_index:(i-1)])
      start_index <- i
    }
  }

  # print(at_values)
  # print(abline_coords)
  # print(groups_values)
  
  text_coords[count_diff+1] <- mean(at_values[start_index:n_boxplots])

  if (display_legend) {
    predicted_width <- strwidth(unique_groups_values, units = "inches", cex = text_size)
    # predicted_height <- strheight(unique_groups_values[1], units = "inches", cex = text_size)
    space_width <- strwidth(" ", units = "inches", cex = text_size)
    number_columns <- min(
      max(
        plt_width %/% (5 * space_width + max(predicted_width)),
        1),
      length(unique_groups_values)
    )
    number_rows <- ceiling(length(unique_groups_values) / number_columns)

    text_height <- strheight(paste(
          rep("TEXT", number_rows + 1),
          collapse = "\n"
          ),
          units = "inches",
          cex = text_size)  + 0.1


    layout(
      matrix(c(1,2), nrow = 2),
      heights = c(
        lcm((plt_height - text_height) * 2.54),
        lcm(text_height * 2.54)
      )
    )
  }

  boxplot(
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
  abline(v = abline_coords, lty = "dashed", col = "grey")
  title(xlab = "k", ylab = "ecc", cex.lab = text_size)
  axis(side = 1, at = text_coords, labels = unique_x_values, las = 2, cex.axis = text_size)

  if (!display_legend) {
    return()
  }

  par(mar = c(0, 0, 0, 0))
  plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend(
    "topleft",
    legend = unique_groups_values,
    col = groups_colours,
    pch = 15,
    cex = text_size,
    pt.cex = text_size * 2,
    bty = 'n',
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

color_ggplot <- function(embedding, color_info, pt_size = 1) {
  df <- data.frame(embedding)
  colnames(df) <- paste("UMAP", 1:2, sep = "_")

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$UMAP_1,
      y = .data$UMAP_2,
    )
  ) +
    ggplot2::geom_point(ggplot2::aes(color = color_info), size = pt_size) +
    ggplot2::theme_bw()
}

metadata_plot <- function(embedding,
                          metadata_name,
                          plt_height,
                          plt_width,
                          groups_highlight = NULL,
                          pt_size = 1,
                          pch = '.',
                          text_size = 1,
                          axis_size = 1,
                          display_legend = FALSE,
                          predicted_height = NULL,
                          labels = FALSE) {


  # print(pkg_env$metadata[[metadata_name]])
  # print(embedding)

  color_plot2(
    embedding = embedding,
    color_info = pkg_env$metadata[[metadata_name]],
    color_values = pkg_env$metadata_colors[[metadata_name]],
    unique_values = pkg_env$metadata_unique[[metadata_name]],
    plt_height = plt_height,
    plt_width = plt_width,
    groups_highlight = groups_highlight,
    pt_size = pt_size,
    pch = pch,
    text_size = text_size,
    axis_size = axis_size,
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
    par(mai = c(0.1, 0, 0.1, 0))
    unique_values <- c(min(color_info), max(color_info))
  } else {
    par(mar = c(0, 0, 0, 0))
    # calculate space needed for the legend
    predicted_width <- strwidth(c(" ", unique_values), units = "inches", cex = text_size)
    space_width <- predicted_width[1]
    predicted_width <- predicted_width[2:length(predicted_width)]
    number_columns <- min(
      max(
        plt_width %/% (5 * space_width + max(predicted_width)),
        1),
      length(unique_values)
    )
  }

  if (is.null(color_values)) {
    color_values <- viridis::viridis(50)
  }

  if (is.function(color_values)) {
    ncolors <- ifelse(is_continuous, 50, length(unique_values))
    color_values <- color_values(ncolors)
  }

  if (is_continuous) {
    legend_image <- as.raster(matrix(color_values, nrow=1))
    plot(c(0, 1) , c(-1, 1), type = "n", axes = FALSE, bty = "n", ylab = "", xlab = "")
    text(y=-0.5, x = seq(0,1,l=5), labels = round(seq(from = unique_values[1], to = unique_values[2], length.out = 5), digits = 3), cex = text_size)
    rasterImage(legend_image, 0, 0, 1,1)
    
  } else {
    plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
    legend(
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
  # par(old_par)

  # dev.off()
}


only_legend_metadata_plot <- function(metadata_name,
                                      groups = NULL,
                                      text_size = 1,
                                      plt_width) {

  unique_values <- pkg_env$metadata_unique[[metadata_name]]
  color_values <- pkg_env$metadata_colors[[metadata_name]]
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
                        groups_highlight = NULL,
                        unique_values = NULL,
                        pt_size = 1,
                        pch = ".",
                        text_size = 1,
                        axis_size = 1,
                        display_legend = FALSE,
                        predicted_height = NULL,
                        labels = FALSE) {

  if (is.null(color_values)) {
    # TODO treat this case
  }


  # xlim <- c(min(embedding[, 1]), max(embedding[, 1]))
  # ylim <- c(min(embedding[, 2]), max(embedding[, 2]))
  
  # convert pixels to inches
  plt_height <- plt_height / ppi 
  plt_width <- plt_width / ppi
  
  is_continuous <- is.null(unique_values) 
  
  if (display_legend) {
    if(is_continuous) {
      unique_values <- c(min(color_info), max(color_info))
      number_rows <- 2

      if (is.null(predicted_height)) {
        predicted_height <- strheight("1", units = "inches", cex = text_size)
      }
    } else {
        # calculate space needed for the legend
      predicted_width <- strwidth(unique_values, units = "inches", cex = text_size)
      space_width <- strwidth(" ", units = "inches", cex = text_size)
      number_columns <- min(
        max(
          plt_width %/% (5 * space_width + max(predicted_width)),
          1),
        length(unique_values)
      )


      if (is.null(predicted_height)) {
        number_rows <- ceiling(length(unique_values) / number_columns)
        predicted_height <- strheight(paste(
              rep("TEXT", number_rows + 1),
              collapse = "\n"
              ),
              units = "inches",
              cex = text_size)
      }

    }

  }

  if (is.null(color_values)) {
    color_values <- viridis::viridis(50)
  }

  if (is.function(color_values)) {
    ncolors <- ifelse(is_continuous, 50, length(unique_values))
    color_values <- color_values(ncolors)
  }

  if (is_continuous) {
    color_info <- cut(color_info, breaks = 50)
    groups_included <- NULL
  } else {
    groups_included <- seq_along(unique_values)

    if (!is.null(groups_highlight)) {
      groups_included <- which(unique_values %in% groups_highlight)
      color_values[-groups_included] <- "lightgray"

      if (length(groups_included) != length(unique_values)) {
        mask <- (color_info %in% groups_highlight)
        embedding <- rbind(
          embedding[!mask, ],
          embedding[mask, ])
        color_info <- c(color_info[!mask], color_info[mask])
      }
    }
  }

  if (display_legend) {
    layout(
      matrix(c(1,2), nrow = 2),
      heights = c(
        lcm(plt_height*2.54),
        lcm(predicted_height *2.54)
      )
    )
  }

  if (is.logical(color_info)) {
    embedding <- rbind(
      embedding[!color_info, ],
      embedding[color_info, ]
    )

    ntrue <- sum(color_info)
    colrs <- c(
      rep(color_values[1], length(color_info) - ntrue),
      rep(color_values[2], ntrue)
    )
    # colrs[color_info] <- color_values[2]
  } else {
    colrs <- color_values[color_info]
  }

  plot(
    embedding,
    pch = pch,
    col = colrs,
    cex = pt_size,
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    cex.axis = axis_size,
    cex.lab = axis_size
  )

  if (!is_continuous && labels) {
      for (unique_val in unique_values[groups_included]) {
        mask <- (color_info == unique_val)

        # TODO precalculate the medians in the write_object funcrion
        geometric_median <- Gmedian::Gmedian(embedding[mask, ])

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
    par(mai = c(0.1, 0, 0.1, 0))
    legend_image <- as.raster(matrix(color_values, nrow=1))
    plot(c(0,1),c(-1,1), type = 'n', axes = F, bty='n',ylab='',xlab='')
    text(y=-0.5, x = seq(0,1,l=5), labels = round(seq(from = unique_values[1], to = unique_values[2], length.out = 5), digits = 3), cex = text_size)
    rasterImage(legend_image, 0, 0, 1,1)
    
  } else {
    
    par(mar = c(0, 0, 0, 0))
    plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend(
      "topleft",
      legend = unique_values,
      col = color_values,
      pch = 15,
      cex = text_size,
      pt.cex = text_size * 2,
      bty = 'n',
      ncol = number_columns,
      # text.width = NA,
      xpd = TRUE
    )
  }

  # dev.off()
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
    if(is.factor(color_info)) {
      unique_values <- levels(color_info)
    } else {
      unique_values <- unique(color_info)
    }
    ncolors <- length(unique_values)
    paletteer_function <- paletteer::paletteer_d
  } else {
    color_info <- cut(color_info, breaks = ncolors)
  }
  colors <- paletteer_function(palette = color_palette,
                                   n = ncolors)

  return({
    plot(
      embedding[, 1],
      embedding[, 2],
      pch = pch,
      col = colors[color_info],
      cex = pt_size,
      panel.first = {
        axis(1, tck = 1, lty = 2, col = "gray")
        axis(2, tck = 1, lty = 2, col = "gray")
      },
      xlab = "UMAP_1",
      ylab = "UMAP_2"
    )

    if (labels) {
      for (unique_val in unique_values) {
        geometric_median <- Gmedian::Gmedian(embedding[color_info == unique_val, ])

        text_width <- strwidth(unique_val, cex = text_size)
        rect(
          geometric_median[1],
          geometric_median[2],
          geometric_median[1] + text_width * 1.2,
          geometric_median[2] + text_size * 1.2,
          col = "white",
          border = "black"
        )

        text(
          geometric_median[1] + text_width * 0.6,
          geometric_median[2] + text_size * 0.6,
          unique_val,
          cex = text_size,
          col = "black"
        )
      }

      # title(xlab = "UMAP_1", ylab = "UMAP_2")
    }
    # legend(x = "bottom",
    #    inset = c(0, -0.2), # You will need to fine-tune the second
    #                        # value depending on the windows size
    #    legend = unique_values,
    #    ncol = 10,
    #    pch = 20,
    #    box.lwd = 0,
    #    col = as.character(colors),
    #    pt.cex = pt_size / 2,
    #    xpd = TRUE
    #    )
  })
}

color_c_plot <- function(embedding,
                       color_info,
                       pt_size = 1,
                       pch = ".",
                       color_palette = "viridis::viridis",
                       text_size = 1) {

  colors <- paletteer::paletteer_c(palette = color_palette,
                                   n = 50)
  
  return({
    plot(
      embedding[, 1],
      embedding[, 2],
      pch = pch,
      col = colors[color_info],
      cex = pt_size,
      panel.first = {
        axis(1, tck = 1, lty = 2, col = "gray")
        axis(2, tck = 1, lty = 2, col = "gray")
      },
      xlab = "UMAP_1",
      ylab = "UMAP_2"
    )


      # title(xlab = "UMAP_1", ylab = "UMAP_2")
    
    # legend(x = "bottom",
    #    inset = c(0, -0.2), # You will need to fine-tune the second
    #                        # value depending on the windows size
    #    legend = unique_values,
    #    ncol = 10,
    #    pch = 20,
    #    box.lwd = 0,
    #    col = as.character(colors),
    #    pt.cex = pt_size / 2,
    #    xpd = TRUE
    #    )
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
  features = NULL
) {
  print("a")
  print(mean_functions[[mean_function_name]])
  print(cells1_name)
  print(cells2_name)
  foldchange_result <- Seurat::FoldChange(
    object = expression_matrix,
    cells.1 = cells1_name,
    cells.2 = cells2_name,
    features = features,
    mean.fxn = mean_functions[[mean_function_name]],
    fc.name = ifelse(mean_function_name == "rowMeans", "avg_diff", "avg_log2FC")
  )
  print("b")

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
  calculate_dim_reduction = TRUE
) {
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

  for(nclust in nclusters) {
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

#' Calculate markers
#'
#' @description to be completed
#'
#' @export
calculate_markers <- function(expression_matrix,
                              cells1,
                              cells2,
                              logfc_threshold = 0,
                              min_pct_threshold = 0.1,
                              min_diff_pct_threshold = -Inf,
                              rank_matrix = NULL,
                              feature_names = NULL,
                              pseudocount_use = 1,
                              base = 2) {
  
  print(Sys.time())
  default_mean_fxn <- function(x) {
    return(log(x = rowMeans(x = x) + pseudocount_use, base = base))
  }

  if (is.null(feature_names)) {
    feature_names <- rownames(expression_matrix)
  }

  cells2 <- setdiff(cells2, cells1)
  
  indices <- c(cells1, cells2)
  n <- length(feature_names)

  expression_matrix <- expression_matrix[ , indices]


  used_slot <- "data"
  norm_method <- ""

  base_text <- ifelse(test = base == exp(1), yes = "", no = base)
  fc_name <- ifelse(
    test = used_slot == "scale.data", 
    yes = "avg_diff",
    no = paste0("avg_log", base_text, "FC")
  )
  mean_fxn <- switch(EXPR = used_slot,
    data = switch(EXPR = norm_method,
      LogNormalize = function(x) {
        return(log(x = rowMeans(x = expm1(x = x)) + pseudocount_use, base = base))
      },
      default_mean_fxn),
    scale.data = rowMeans,
    default_mean_fxn
  )
  
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
                abs(fc_results[ , 1]) >= logfc_threshold)

  if (length(mask) == 0) {
    return(data.frame())
  }

  expression_matrix <- expression_matrix[mask, ]
  fc_results <- fc_results[mask, ]
  feature_names <- feature_names[mask]

  # return(expression_matrix)

  if (!is.matrix(rank_matrix)) {
    print(Sys.time())
    rank_matrix <- matrix(nrow = nrow(expression_matrix), ncol = ncol(expression_matrix))

    for (i in seq_len(nrow(expression_matrix))) {
      rank_matrix[i, ] <- rank(expression_matrix[i, ], ties.method = "min")
    }
    print(Sys.time())
  } else {
    rank_matrix <- rank_matrix[mask, indices]
  }

  fc_results$p_val <- wilcox_test(rank_matrix, length(cells1), max(rank_matrix))
  print(Sys.time())
  fc_results$gene <- feature_names
  fc_results$p_val_adj <- p.adjust(
    p = fc_results$p_val,
    method = "bonferroni",
    n = n 
  )

  fc_results <- fc_results[, c(5, 1, 2, 3, 4, 6)]

  return(fc_results)
}

# function copied from https://cran.r-project.org/web/packages/TeachingDemos/index.html
shadowtext <- function(x, y = NULL, labels, col = 'white', bg = 'black',
	theta= seq(pi / 32, 2 * pi, length.out = 64), r = 0.1, cex=1, ... ) {

	xy <- xy.coords(x,y)
	fx <- grconvertX(xy$x, to='nfc')
	fy <- grconvertY(xy$y, to='nfc')
	fxo <- r*strwidth('A', units='figure', cex=cex)
	fyo <- r*strheight('A', units='figure', cex=cex)
	
	for (i in theta) {
	  text(grconvertX(fx + cos(i)*fxo, from="nfc"),
	       grconvertY(fy + sin(i)*fyo, from="nfc"),
	       labels, cex=cex, col=bg, ...)
	}
	text(xy$x, xy$y, labels, cex=cex, col=col, ... ) 
}