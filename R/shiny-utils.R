
# expression_matrix <- NA
# metadata <- NA
# genes <- NA
# stab_obj <- NA
# feature_types <- NA

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

empty_ggplot <- function() {
  ggplot2::ggplot() +
    ggplot2::theme_void()
}

expression_ggplot <- function(embedding, expression) {
  ggplot2::ggplot(
    data.frame(embedding),
    ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2, color = expression)
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::scale_color_gradientn("", colors = cList[[1]]) +
    ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 15)) +
    ggplot2::theme(legend.position = "bottom")
}

color_ggplot <- function(embedding, color_info) {
  df <- data.frame(embedding)

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

voting_scheme_ggplot <- function(embedding,
                                 selected_genes,
                                 expression_threshold) {
  umap_df <- data.frame(embedding)
  n_expressed_genes <- rep(0, ncol(expression_matrix))
  expression_threshold <- as.numeric(expression_threshold)
  for (gene in selected_genes) {
    n_expressed_genes <- n_expressed_genes + (expression_matrix[gene, ] > expression_threshold)
  }

  mask <- (n_expressed_genes >= length(selected_genes))

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = umap_df[!mask, ],
      mapping = ggplot2::aes(x = .data$UMAP_1, .data$UMAP_2),
      color = "#C3C3C3"
    ) +
    ggplot2::geom_point(
      data = umap_df[mask, ],
      mapping = ggplot2::aes(x = .data$UMAP_1, .data$UMAP_2),
      color = "red"
    ) +
    ggplot2::theme_bw()
}

metadata_ggplot <- function(embedding,
                            metadata_name) {
  metadata_value <- metadata[, metadata_name]
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

download_modal <- function(session, placeholder, download_id) {
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

download_handler <- function(session, filename, to_save_plot, width, height) {
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
