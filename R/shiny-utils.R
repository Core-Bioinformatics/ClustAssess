
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

color_ggplot <- function(embedding, color_info) {
  df <- data.frame(embedding)
  colnames(df) <- paste("UMAP", 1:2, sep = "_")

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$UMAP_1,
      y = .data$UMAP_2,
    )
  ) +
    ggplot2::geom_point(ggplot2::aes(color = color_info)) +
    ggplot2::theme_bw()
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

  # print(placeholder)

  # output[[download_id]] <- download_plot_handler(
  #   session,
  #   "test.pdf",
  #   empty_ggplot(),
  #   7,
  #   7
  # )

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
  mean_function_name = "logNormalize",
  features = NULL
) {
  print("a")
  print(mean_functions[[mean_function_name]])
  print(cells1_name)
  print(cells2_name)
  foldchange_result <- Seurat::FoldChange(
    object = pkg_env$expression_matrix,
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